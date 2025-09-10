package op;

import exceptions.MusialException;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Tuple;
import model.Storage;
import model.Variant;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.RandomStringUtils;
import org.apache.commons.lang3.tuple.Triple;
import util.Constants;
import util.IO;
import util.Logging;
import util.OS;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Objects;
import java.util.stream.Collectors;

/**
 * The {@link VariantAnnotator} class provides functionality to annotate variants.
 * <p>
 * Currently, it supports annotation using SnpEff, a popular tool for annotating and predicting the effects of genetic variants.
 * <p>
 * This class is responsible for managing the workflow of variant annotation, including preparing input data, executing SnpEff commands, and
 * processing the output to update the storage with annotated variants. It handles temporary file management and error logging to ensure a
 * smooth annotation process.
 */
public class VariantAnnotator {

    /**
     * {@link Storage} object for managing data.
     */
    private final Storage storage;

    /**
     * Constructs a new instance of the {@link VariantAnnotator} class.
     * <p>
     * This constructor initializes the {@link VariantAnnotator} with the provided {@link Storage} object.
     *
     * @param storage The {@link Storage} object used to manage genomic data.
     */
    public VariantAnnotator(Storage storage) {
        this.storage = storage;
    }

    /**
     * Runs variant annotation using SnpEff on the {@link Storage}.
     * <p>
     * This method performs the following steps:
     * <ul>
     *   <li>Checks if the storage contains novel variants to annotate. If not, throws an {@link IllegalArgumentException}.</li>
     *   <li>Creates a temporary directory for SnpEff files and configurations.</li>
     *   <li>Writes the novel variants from the storage to a temporary VCF file.</li>
     *   <li>Copies the SnpEff configuration and JAR files to the temporary directory.</li>
     *   <li>Writes the reference genome and features files to the appropriate locations in the temporary directory.</li>
     *   <li>Updates the SnpEff configuration file with reference genome information.</li>
     *   <li>Builds the SnpEff database using the reference genome information.</li>
     *   <li>Runs the SnpEff annotation on the temporary VCF file.</li>
     *   <li>Processes the annotation results and updates the storage with the annotated data.</li>
     *   <li>Handles errors during the SnpEff build and annotation processes, logging them and saving error logs to the output
     *   directory.</li>
     *   <li>Cleans up the temporary directory after the analysis is complete.</li>
     * </ul>
     *
     * @param workingDirectory The working directory where temporary files will be created.
     * @throws MusialException If the SnpEff annotation process fails.
     * @throws IOException     If an error occurs while reading or writing files.
     */
    public void runSnpEff(Path workingDirectory) throws MusialException, IOException {
        // Generate a temporary directory for snpEff.
        String prefix = "%s-%s".formatted("snpEff", RandomStringUtils.randomAlphanumeric(6));
        Path temp = Files.createTempDirectory(prefix);
        try {
            // Create map to store variant pointers and write storage variants to temporary VCF file.
            Collection<Triple<String, Integer, String>> novelVariants = storage.getNovelVariants();
            ArrayList<Tuple<Triple<String, Integer, String>, Variant>> variants = new ArrayList<>(novelVariants.size());
            for (Triple<String, Integer, String> variant : novelVariants) {
                variants.add(new Tuple<>(variant, storage.getContig(variant.getLeft()).getVariant(variant.getMiddle(),
                        variant.getRight())));
            }
            IO.writeFile(Path.of(temp + "/variants" + FileExtensions.VCF), StorageIO.toVCF(variants));

            // Write reference .gff and .fasta to temp. target directory.
            IO.writeFile(Path.of(temp + "/data/reference/genes.gff"), StorageIO.toGFF3(storage));
            IO.writeFile(Path.of(temp + "/data/genomes/reference.fa"), StorageIO.toFASTA(storage));

            // Copy snpEff config and JAR to temp. target directory.
            Path snpEffConfigPath = Path.of(temp + "/snpEff.config");
            IO.copyResourceToFile("/snpEff/snpEff.config", snpEffConfigPath);
            IO.copyResourceToFile("/snpEff/snpEff.jar", Path.of(temp + "/snpEff.jar"));

            // Add reference .fasta and .gff information to snpEff.config.
            String codonTableConfig = storage.getContigs()
                    .stream()
                    .map(contig -> "reference.genome.%s : Bacterial_and_Plant_Plastid".formatted(contig._id))
                    .collect(Collectors.joining("\n"));
            Files.writeString(snpEffConfigPath, "\n# reference genome\nreference.genome : reference\n%s".formatted(codonTableConfig),
                    StandardOpenOption.APPEND);

            // Generate database with reference genome information.
            String[] cmdSnpEffBuild = {"java", "-jar", "snpEff.jar", "build", "-gff3", "-noLog", "-nodownload", "-maxErrorRate", "0.0",
                    "-noCheckCds", "-noCheckProtein", "reference"};
            OS.runCommand(cmdSnpEffBuild, temp + "/snpEff.build.err", temp + "/snpEff.build.log",
                    temp.toString());

            // Run snpEff annotation on variants file.
            String[] cmdSnpEffAnn = {"java", "-jar", "snpEff.jar", "eff", "-noLog", "-noStats", "-nodownload", "-noShiftHgvs",
                    "-noHgvs", "-ud", "0", "reference", temp + "/variants" + FileExtensions.VCF};
            OS.runCommand(cmdSnpEffAnn, temp + "/snpEff.ann.err", temp + "/annotation" + FileExtensions.VCF,
                    temp.toString());

            // Transfer annotation results to storage.
            try (BufferedReader br = new BufferedReader(new FileReader(temp + "/annotation" + FileExtensions.VCF))) {
                String line = br.readLine();
                int index = 0;
                while (Objects.nonNull(line)) {
                    if (!line.startsWith(Constants.SIGN)) {
                        String[] annotationFields = line.split("\t");
                        if (!annotationFields[7].equals(".")) {
                            annotationFields = annotationFields[7].replace("ANN=", "").split(Constants.COMMA)[0].split("\\|");
                            for (int i = 0; i < annotationFields.length; i++) {
                                if (i == 1 || i == 2 || i == 5 || i == 7 || i == 12 || i == 13) {
                                    variants.get(index).b.addAttributeIfAbsent(
                                            Constants.SNP_EFF_PREFIX + Constants.SNP_EFF_KEYS.get(i),
                                            i == 1 ? annotationFields[i].replaceAll("&", Constants.COMMA) : annotationFields[i]
                                    );
                                } else if (i == 6) {
                                    variants.get(index).b.addAttributeIfAbsent(
                                            Constants.SNP_EFF_PREFIX + Constants.SNP_EFF_KEYS.get(i),
                                            annotationFields[i].split("-")[1]
                                    );
                                }
                            }
                        }
                        index++;
                    }
                    line = br.readLine();
                }
            } catch (FileNotFoundException e) {
                throw new MusialException(String.format("Failed to read SnpEff annotation. %s", e.getMessage()));
            }
        } finally {
            File buildErrorFile = new File(temp + "/snpEff.build.err");
            if (buildErrorFile.exists() && buildErrorFile.length() != 0) {
                Logging.logSevere("SnpEff `build` has raised an error or warning; a copy of the log file is in the output directory -" +
                        " the annotations may be incorrect.");
                FileUtils.copyFile(buildErrorFile, new File(workingDirectory.toAbsolutePath()
                        + "/musial_snpeff_build_%s.error".formatted(Logging.getDate())));
            }
            File annErrorFile = new File(temp + "/snpEff.ann.err");
            if (annErrorFile.exists() && annErrorFile.length() != 0) {
                Logging.logSevere("SnpEff `ann` has raised an error or warning; a copy of the log file is in the output directory - " +
                        "the annotations may be incorrect.");
                FileUtils.copyFile(annErrorFile, new File(workingDirectory.toAbsolutePath()
                        + "/musial_snpeff_ann_%s.error".formatted(Logging.getDate())));
            }
            // Clean up temporary directory.
            FileUtils.deleteDirectory(temp.toFile());
        }
    }

}
