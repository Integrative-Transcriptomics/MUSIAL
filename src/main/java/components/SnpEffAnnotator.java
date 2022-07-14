package components;

import exceptions.MusialIOException;
import main.Musial;
import utility.CL;
import utility.IO;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.Objects;

/**
 * Runs a local installation of SnpEff via command line.
 * <p>
 * Used to run SnpEff .vcf file annotations.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class SnpEffAnnotator {

  @SuppressWarnings("DuplicateExpressions")
  public static void runSnpEff(File targetDir, File targetVcf, File referenceFasta, File referenceGff,
                               String referenceChromosome)
      throws IOException, MusialIOException {
    // 1. Copy snpEff from JAR to target directory.
    if (targetDir.isDirectory()) {
      InputStream snpEffConfig = SnpEffAnnotator.class.getResourceAsStream("/snpEff/snpEff.config");
      Files.copy(
          Objects.requireNonNull(snpEffConfig),
          Path.of(targetDir + "/snpEff.config"),
          StandardCopyOption.REPLACE_EXISTING
      );
      snpEffConfig.close();
      InputStream snpEffJar = SnpEffAnnotator.class.getResourceAsStream("/snpEff/snpEff.jar");
      Files.copy(
          Objects.requireNonNull(snpEffJar),
          Path.of(targetDir + "/snpEff.jar"),
          StandardCopyOption.REPLACE_EXISTING
      );
      snpEffJar.close();
    } else {
      throw new IOException("Failed to copy SnpEff.jar to " + targetDir.getAbsolutePath());
    }
    // 2. Copy reference .gff and .fasta into snpEff/data directory.
    IO.generateDirectory(new File(targetDir.toPath() + "/data/"));
    IO.generateDirectory(new File(targetDir.toPath() + "/data/" + referenceChromosome));
    Files.copy(
        referenceGff.toPath(),
        Path.of(targetDir.toPath() + "/data/" + referenceChromosome + "/genes.gff"),
        StandardCopyOption.REPLACE_EXISTING
    );
    IO.generateDirectory(new File(targetDir.toPath() + "/data/genomes/"));
    Files.copy(
        referenceFasta.toPath(),
        Path.of(targetDir.toPath() + "/data/genomes/" + referenceChromosome + ".fa"),
        StandardCopyOption.REPLACE_EXISTING
    );
    // 3. Add reference .fasta and .gff information to snpEff.config.
    Files.writeString(
        Path.of(targetDir + "/snpEff.config"),
        "\n# " + referenceChromosome + " genome",
        StandardOpenOption.APPEND
    );
    Files.writeString(
        Path.of(targetDir + "/snpEff.config"),
        "\n" + referenceChromosome + ".genome : " + referenceChromosome,
        StandardOpenOption.APPEND
    );

    // 4. Generate database with reference genome information.
    String[] createSNPEffDatabase = {"java", "-jar", "snpEff.jar", "build", "-gff3", "-v", referenceChromosome};
    CL.runCommand(createSNPEffDatabase, targetDir + "/" + "SnpEffDatabase.log", "",
        targetDir.getAbsolutePath());

    // 5. Run snpEff annotation on each input .vcf file (multi-threaded).
    String[] runSNPEff = {
        "java",
        "-jar",
        "snpEff.jar",
        "-noLog",
        "-t",
        String.valueOf(Musial.THREADS),
        "-v",
        "-no-downstream",
        "-no-upstream",
        "-no-intron",
        "-no-intergenic",
        "-no",
        "INTRAGENIC",
        /* TODO: Handle generation of statistics using SnpEff.
          "-s", new File(targetDir + "/" + "SnpEff.html").getAbsolutePath(),
         */
        referenceChromosome,
        targetVcf.getAbsolutePath()
    };
    CL.runCommand(runSNPEff,
        new File(targetDir + "/" + "SnpEff.log").getAbsolutePath(),
        new File(targetDir + "/" + "SnpEff.vcf").getAbsolutePath(),
        targetDir.getAbsolutePath());

  }

  public static HashMap<String, String> convertAnnotation(String snpEffAnnotationString) {
    HashMap<String, String> annotationsMap = new HashMap<>();
    String[] snpEffAnnotations = snpEffAnnotationString.split(",", -1);
    String[] snpEffAnnotationFields;
    String featureName;
    String effect;
    String impact;
    String HGVSp;
    for (String snpEffAnnotation : snpEffAnnotations) {
      snpEffAnnotationFields = snpEffAnnotation.split("\\|", -1);
      featureName = snpEffAnnotationFields[3];
      effect = snpEffAnnotationFields[1];
      impact = snpEffAnnotationFields[2];
      HGVSp = snpEffAnnotationFields[10];
      annotationsMap.put(featureName + "_EFF", effect + "|" + impact + "|" + HGVSp);
    }
    return annotationsMap;
  }

}
