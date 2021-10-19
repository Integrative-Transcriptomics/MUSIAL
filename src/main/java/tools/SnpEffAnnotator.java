package tools;

import exceptions.MusialIOException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.Objects;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import main.Musial;
import me.tongfei.progressbar.ProgressBar;
import org.apache.commons.io.FilenameUtils;
import runnables.SnpEffAnnotatorRunnable;
import utility.CL;
import utility.IO;
import utility.Logging;

/**
 * Runs local installation of SnpEff via command line.
 * <p>
 * Used to run SnpEff .vcf file annotations.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class SnpEffAnnotator {

  /**
   * Conducts SnpEff annotation on each sample.
   * <p>
   * A version of SnpEff is contained within the JAR locally. The snpEff.jar and snpEff.config files are copied into
   * a temporary directory in order to run. Next, the configuration file is adjusted and the reference sequence and
   * annotation files are copied into a temporary directory (as it is required for SnpEff in order to execute). After
   * generating a database with the specified reference files, SnpEff is run on each single sample .vcf file.
   *
   * @param arguments Instance of {@link ArgumentsParser} containing the command line options specified by the user.
   * @throws IOException          If any copy procedure or writing to a file fails.
   * @throws MusialIOException    If any generation of output directories or files fails.
   * @throws InterruptedException If a single thread was interrupted.
   */
  public static void runSnpEff(ArgumentsParser arguments) throws IOException, MusialIOException,
      InterruptedException {
    ProgressBar progress = Musial.buildProgress();
    progress.maxHint(4 + arguments.getSampleInput().size());
    try (progress) {
      // 1. Copy snpEff from internal installation to output directory.
      String snpEffPath = arguments.getOutputDirectory().toString() + "/snpEff/";
      File snpEffFile = new File(snpEffPath);
      if (snpEffFile.mkdirs()) {
        InputStream snpEffConfig = SnpEffAnnotator.class.getResourceAsStream("/snpEff/snpEff.config");
        Files.copy(
            Objects.requireNonNull(snpEffConfig),
            Path.of(snpEffPath + "snpEff.config"),
            StandardCopyOption.REPLACE_EXISTING
        );
        snpEffConfig.close();
        InputStream snpEffJar = SnpEffAnnotator.class.getResourceAsStream("/snpEff/snpEff.jar");
        Files.copy(
            Objects.requireNonNull(snpEffJar),
            Path.of(snpEffPath + "snpEff.jar"),
            StandardCopyOption.REPLACE_EXISTING
        );
        snpEffJar.close();
      } else {
        throw new MusialIOException("Unable to create SnpEff output directory:\t" + snpEffFile.getAbsolutePath());
      }
      progress.step();
      // 2. Copy reference .fasta and .gff into snpEff/data directory.
      String referenceName = FilenameUtils.removeExtension(arguments.getReferenceInput().getName());
      String snpEffDataPath = snpEffPath + "data/" + referenceName + "/";
      File snpEffDataFile = new File(snpEffDataPath);
      if (snpEffDataFile.mkdirs()) {
        Files.copy(
            arguments.getAnnotationInput().toPath(),
            Path.of(snpEffDataPath + "genes.gff"),
            StandardCopyOption.REPLACE_EXISTING
        );
      } else {
        throw new MusialIOException("Unable to create SnpEff output directory:\t" + snpEffDataFile.getAbsolutePath());
      }
      snpEffDataPath = snpEffPath + "data/genomes/";
      snpEffDataFile = new File(snpEffDataPath);
      if (snpEffDataFile.mkdirs()) {
        Files.copy(
            arguments.getReferenceInput().toPath(),
            Path.of(snpEffDataPath + referenceName + ".fa"),
            StandardCopyOption.REPLACE_EXISTING
        );
      } else {
        throw new MusialIOException("Unable to create SnpEff output directory:\t" + snpEffDataFile.getAbsolutePath());
      }
      progress.step();
      // 3. Add reference .fasta and .gff information to snpEff.config.
      Files.write(
          Path.of(snpEffPath + "snpEff.config"),
          ("\n# " + referenceName + " genome").getBytes(StandardCharsets.UTF_8),
          StandardOpenOption.APPEND
      );
      Files.write(
          Path.of(snpEffPath + "snpEff.config"),
          ("\n" + referenceName + ".genome : " + referenceName).getBytes(StandardCharsets.UTF_8),
          StandardOpenOption.APPEND
      );
      progress.step();
      // 4. Generate database with reference genome information.
      String[] createSNPEffDatabase = {"java", "-jar", "snpEff.jar", "build", "-gff3", "-v", referenceName};
      CL.runCommand(createSNPEffDatabase, snpEffPath + "/" + referenceName + "_snpEffDatabase.log", "", snpEffPath);
      progress.step();
      // 5. Run snpEff annotation on each input .vcf file (multi-threaded).
      ExecutorService executor = Executors.newFixedThreadPool(arguments.getNumThreads());
      String snpEffResultsPath = snpEffPath + "results/";
      File snpEffResultsFile = new File(snpEffResultsPath);
      if (snpEffResultsFile.mkdirs()) {
        for (int i = 0; i < arguments.getSampleInput().size(); i++) {
          executor.execute(new SnpEffAnnotatorRunnable(
              snpEffPath,
              snpEffResultsPath,
              arguments.getSampleNames().get(i),
              i,
              referenceName,
              arguments.getSampleInput().get(i),
              arguments,
              progress
          ));
        }
        executor.shutdown();
        //noinspection ResultOfMethodCallIgnored
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
      } else {
        throw new MusialIOException(
            "Unable to create SnpEff output directory:\t" + snpEffResultsFile.getAbsolutePath());
      }
      // 6. Delete temporarily generated files.
      IO.deleteDirectory(new File(snpEffPath + "data/"));
      //noinspection ResultOfMethodCallIgnored
      new File(snpEffPath + "snpEff.config").delete();
      //noinspection ResultOfMethodCallIgnored
      new File(snpEffPath + "snpEff.jar").delete();
      //noinspection ResultOfMethodCallIgnored
      new File(snpEffPath + referenceName + "_snpEffDatabase.log").delete();
      progress.setExtraMessage(Logging.getDoneMessage());
    }
  }

}
