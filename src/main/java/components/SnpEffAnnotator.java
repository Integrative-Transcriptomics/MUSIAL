package components;

import datastructure.FastaContainer;
import exceptions.MusialException;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Objects;
import java.util.stream.Collectors;

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
     * @param targetDir      {@link File} instance pointing to the directory at which SnpEff shall be run.
     * @param targetVcf      {@link File} instance pointing to the vcf file to annotate.
     * @param referenceFasta {@link File} instance pointing to the fasta file containing the reference genome.
     * @param referenceGff   {@link File} instance pointing to the gff file containing the reference genome annotation.
     * @throws IOException     Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialException Thrown if any input or output file is missing or unable to being generated (caused by any MUSIAL method).
     */
    @SuppressWarnings({"DuplicateExpressions", "unused"})
    public static void runSnpEff(File targetDir, File targetVcf, File referenceFasta, File referenceGff)
            throws IOException, MusialException {
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
        IO.generateDirectory(new File(targetDir.toPath() + "/data/reference/"));
        Files.copy(
                referenceGff.toPath(),
                Path.of(targetDir.toPath() + "/data/reference/genes.gff"),
                StandardCopyOption.REPLACE_EXISTING
        );
        IO.generateDirectory(new File(targetDir.toPath() + "/data/genomes/"));
        Files.copy(
                referenceFasta.toPath(),
                Path.of(targetDir.toPath() + "/data/genomes/reference.fa"),
                StandardCopyOption.REPLACE_EXISTING
        );
        // 3. Add reference .fasta and .gff information to snpEff.config.
        Files.writeString(
                Path.of(targetDir + "/snpEff.config"),
                "\n# reference genome",
                StandardOpenOption.APPEND
        );
        HashSet<FastaContainer> referenceFastaEntries = IO.readFastaToSet(referenceFasta);
        Files.writeString(
                Path.of(targetDir + "/snpEff.config"),
                "\nreference.genome : reference",
                StandardOpenOption.APPEND
        );

        // 4. Generate database with reference genome information.
        String[] createSNPEffDatabase = {"java", "-jar", "snpEff.jar", "build", "-gff3", "-v", "reference"};
        CL.runCommand(createSNPEffDatabase, targetDir + "/" + "snpEffDatabase.log", "",
                targetDir.getAbsolutePath());

        // 5. Run snpEff annotation on each input .vcf file (multi-threaded).
        String[] runSNPEff = {
                "java",
                "-jar",
                "snpEff.jar",
                "-v",
                "-no-downstream", // Do not show DOWNSTREAM changes
                "-no-intergenic", // Do not show INTERGENIC changes
                "-no-intron", // Do not show INTRON changes
                "-no-upstream", // Do not show UPSTREAM changes
                "-s",
                new File(targetDir + "/" + "snpEff.html").getAbsolutePath(),
                "reference",
                targetVcf.getAbsolutePath()
        };
        CL.runCommand(runSNPEff,
                new File(targetDir + "/" + "snpEff.log").getAbsolutePath(),
                new File(targetDir + "/" + "annotated_variants.vcf").getAbsolutePath(),
                targetDir.getAbsolutePath()
        );
    }

}
