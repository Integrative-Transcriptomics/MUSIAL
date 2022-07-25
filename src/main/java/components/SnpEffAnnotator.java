package components;

import exceptions.MusialIOException;
import main.Musial;

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
 * @version 2.1
 * @since 2.0
 */
public final class SnpEffAnnotator {

    /**
     * Initializes all necessary configuration files and runs SnpEff on a single .vcf file.
     *
     * @param targetDir           {@link File} instance pointing to the directory at which SnpEff shall be run.
     * @param targetVcf           {@link File} instance pointing to the vcf file to annotate.
     * @param referenceFasta      {@link File} instance pointing to the fasta file containing the reference genome.
     * @param referenceGff        {@link File} instance pointing to the gff file containing the reference genome annotation.
     * @param referenceChromosome {@link String} storing the reference chromosome name.
     * @throws IOException       Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialIOException Thrown if any input or output file is missing or unable to being generated (caused by any MUSIAL method).
     */
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
        CL.runCommand(createSNPEffDatabase, targetDir + "/" + "snpEffDatabase.log", "",
                targetDir.getAbsolutePath());

        // 5. Run snpEff annotation on each input .vcf file (multi-threaded).
        String[] runSNPEff = {
                "java",
                "-jar",
                "snpEff.jar",
                "-v",
                "-no-downstream",
                "-no-upstream",
                "-no-intron",
                "-no-intergenic",
                "-s",
                new File(targetDir + "/" + "snpEff.html").getAbsolutePath(),
                referenceChromosome,
                targetVcf.getAbsolutePath()
        };
        CL.runCommand(runSNPEff,
                new File(targetDir + "/" + "snpEff.log").getAbsolutePath(),
                new File(targetDir + "/" + "annotated_variants.vcf").getAbsolutePath(),
                targetDir.getAbsolutePath()
        );
    }

    /**
     * Converts a SnpEff annotation string into a {@link HashMap} object.
     *
     * @param snpEffAnnotationString {@link String} containing the content of a SnpEff annotation.
     * @return {@link HashMap} containing a key/value pair for each SnpEff annotation of the passed {@link String}.
     */
    @SuppressWarnings("unused")
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
