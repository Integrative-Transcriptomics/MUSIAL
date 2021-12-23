package components;

import exceptions.MusialIOException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
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
   * Static {@link String} used as column headers for the SnpEff summary output file.
   */
  public static final String snpEffSummaryHeader =
      "Position\tReference\tAllele\tSamples\tAnnotation" +
          "\tAnnotation_Impact\tGene_Name\tGene_ID\tFeature_Type\tFeature_ID\tTranscript_BioType\tRank\tHGVS.c\tHGVS" +
          ".p\tcDNA.pos/cDNA.length\tCDS.pos/CDS.length\tAA.pos/AA.length\tDistance\tERRORS/WARNINGS/INFO" +
          IO.LINE_SEPARATOR;
  /**
   * Static {@link String} used to access the SnpEff annotations in an .vcf file.
   */
  public static final String snpEffAnnotationAttributeKey = "ANN";

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
    String snpEffPath;
    String referenceName;
    String snpEffDataPath;
    File snpEffDataFile;
    ProgressBar progress1 = Musial.buildProgress();
    progress1.maxHint(1);
    progress1.stepTo(0);
    Logging.logStatus("Generating SnpEff working directory.");
    try (progress1) {
      // 1. Copy snpEff from internal installation to output directory.
      snpEffPath = arguments.getOutputDirectory().toString() + "/SnpEff/";
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
        throw new MusialIOException("Unable to create SnpEff working directory:\t" + snpEffFile.getAbsolutePath());
      }
      // 2. Copy reference .fasta and .gff into snpEff/data directory.
      referenceName = FilenameUtils.removeExtension(arguments.getReferenceFile().getName());
      snpEffDataPath = snpEffPath + "data/" + referenceName + "/";
      snpEffDataFile = new File(snpEffDataPath);
      if (snpEffDataFile.mkdirs()) {
        Files.copy(
            arguments.getAnnotationInput().toPath(),
            Path.of(snpEffDataPath + "genes.gff"),
            StandardCopyOption.REPLACE_EXISTING
        );
      } else {
        throw new MusialIOException("Unable to create SnpEff working directory:\t" + snpEffDataFile.getAbsolutePath());
      }
      snpEffDataPath = snpEffPath + "data/genomes/";
      snpEffDataFile = new File(snpEffDataPath);
      if (snpEffDataFile.mkdirs()) {
        Files.copy(
            arguments.getReferenceFile().toPath(),
            Path.of(snpEffDataPath + referenceName + ".fa"),
            StandardCopyOption.REPLACE_EXISTING
        );
      } else {
        throw new MusialIOException("Unable to create SnpEff working directory:\t" + snpEffDataFile.getAbsolutePath());
      }
      progress1.step();
      progress1.setExtraMessage(Logging.getDoneMessage());
    }
    ProgressBar progress2 = Musial.buildProgress();
    progress2.maxHint(1);
    progress2.stepTo(0);
    Logging.logStatus("Generating SnpEff configuration and database.");
    try (progress2) {
      // 3. Add reference .fasta and .gff information to snpEff.config.
      Files.writeString(
          Path.of(snpEffPath + "snpEff.config"),
          "\n# " + referenceName + " genome",
          StandardOpenOption.APPEND
      );
      Files.writeString(
          Path.of(snpEffPath + "snpEff.config"),
          "\n" + referenceName + ".genome : " + referenceName,
          StandardOpenOption.APPEND
      );
      // 4. Generate database with reference genome information.
      String[] createSNPEffDatabase = {"java", "-jar", "snpEff.jar", "build", "-gff3", "-v", referenceName};
      CL.runCommand(createSNPEffDatabase, snpEffPath + "/" + referenceName + "_snpEffDatabase.log", "", snpEffPath);
      progress2.step();
      progress2.setExtraMessage(Logging.getDoneMessage());
    }
    ProgressBar progress3 = Musial.buildProgress();
    progress3.maxHint(arguments.getSampleInput().size());
    progress3.stepTo(0);
    Logging.logStatus("Running SnpEff on sample variant call files (" + arguments.getSampleInput().size() + ").");
    try (progress3) {
      // 5. Run snpEff annotation on each input .vcf file (multi-threaded).
      ExecutorService executor = Executors.newFixedThreadPool(arguments.getNumThreads());
      String snpEffResultsPath = snpEffPath + "AnnotatedVCFs/";
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
              progress3
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
      progress3.setExtraMessage(Logging.getDoneMessage());
    }
  }

}
