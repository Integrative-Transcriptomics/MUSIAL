package utility;

import datastructure.MusialStorage;
import datastructure.VariantInformation;
import exceptions.MusialException;
import main.Constants;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

/**
 * Implements static methods to run SnpEff locally via command line.
 *
 * @author Simon Hackl
 * @version 2.3
 * @since 2.0
 */
public final class SnpEffAnnotator {

    /**
     * The target directory for SnpEff annotation.
     */
    private final File targetDirectory;

    /**
     * The target .vcf format file to annotate.
     */
    private final File targetVcf;

    /**
     * The reference sequences of the target .vcf file in .fasta format.
     */
    private final File referenceSequences;

    /**
     * The reference sequence features of the target .vcf file in .gff format.
     */
    private final File referenceFeatures;

    /**
     * The MusialStorage instance to transfer annotations to.
     */
    private final MusialStorage musialStorage;

    /**
     * Prefix for exceptions thrown by the SnpEff Annotation process.
     */
    @SuppressWarnings("FieldCanBeLocal")
    private final String EXCEPTION_PREFIX = "(SnpEff Annotation)";

    /**
     * Constructor of {@link SnpEffAnnotator} class.
     *
     * @param targetDirectory    The working directory of SnpEff.
     * @param targetVcf          The .vcf format file to annotate.
     * @param referenceSequences The reference sequence of the target .vcf file in .fasta format.
     * @param referenceFeatures  The reference sequence features of the target .vcf file in .gff format.
     * @param musialStorage      The {@link MusialStorage} instance to transfer annotations to.
     */
    public SnpEffAnnotator(File targetDirectory, File targetVcf, File referenceSequences, File referenceFeatures, MusialStorage musialStorage) {
        this.targetDirectory = targetDirectory;
        this.targetVcf = targetVcf;
        this.referenceSequences = referenceSequences;
        this.referenceFeatures = referenceFeatures;
        this.musialStorage = musialStorage;
    }

    /**
     * Runs the SnpEff variant call annotation process.
     *
     * @throws MusialException If an error occurs during the annotation process.
     * @throws IOException     If an I/O error occurs.
     */
    @SuppressWarnings("ResultOfMethodCallIgnored")
    public void run() throws MusialException, IOException {
        boolean directoryCreated = this.targetDirectory.mkdirs();
        if (!directoryCreated) {
            throw new MusialException(EXCEPTION_PREFIX + " Failed to generate temporary directory at " + this.targetDirectory.getAbsolutePath());
        }
        try {
            // Write musial storage variants to temporary .vcf file.
            IO.writeFile(
                    this.targetVcf,
                    musialStorage.toVcf()
            );
            // Copy snpEff from JAR to target directory.
            try {
                InputStream snpEffConfig = SnpEffAnnotator.class.getResourceAsStream("/snpEff/snpEff.config");
                Files.copy(
                        Objects.requireNonNull(snpEffConfig),
                        Path.of(this.targetDirectory + "/snpEff.config"),
                        StandardCopyOption.REPLACE_EXISTING
                );
                snpEffConfig.close();
                InputStream snpEffJar = SnpEffAnnotator.class.getResourceAsStream("/snpEff/snpEff.jar");
                Files.copy(
                        Objects.requireNonNull(snpEffJar),
                        Path.of(this.targetDirectory + "/snpEff.jar"),
                        StandardCopyOption.REPLACE_EXISTING
                );
                snpEffJar.close();
            } catch (IOException e) {
                throw new MusialException(EXCEPTION_PREFIX + " Failed to extract SnpEff to " + this.targetDirectory.getAbsolutePath());
            }
            // Copy reference .gff and .fasta into /data directory.
            try {
                new File(this.targetDirectory.toPath() + "/data/reference/").mkdirs();
                Files.copy(
                        this.referenceFeatures.toPath(),
                        Path.of(this.targetDirectory.toPath() + "/data/reference/genes.gff"),
                        StandardCopyOption.REPLACE_EXISTING
                );
                new File(this.targetDirectory.toPath() + "/data/genomes/").mkdirs();
                Files.copy(
                        this.referenceSequences.toPath(),
                        Path.of(this.targetDirectory.toPath() + "/data/genomes/reference.fa"),
                        StandardCopyOption.REPLACE_EXISTING
                );
            } catch (IOException e) {
                throw new MusialException(EXCEPTION_PREFIX + " Failed to extract reference .fasta and/or .gff.");
            }
            // Add reference .fasta and .gff information to snpEff.config.
            try {
                Files.writeString(
                        Path.of(this.targetDirectory + "/snpEff.config"),
                        "\n# reference genome",
                        StandardOpenOption.APPEND
                );
                Files.writeString(
                        Path.of(this.targetDirectory + "/snpEff.config"),
                        "\nreference.genome : reference",
                        StandardOpenOption.APPEND
                );
            } catch (IOException e) {
                throw new MusialException(EXCEPTION_PREFIX + " Failed to add reference  " + this.targetDirectory.getAbsolutePath());
            }
            // Generate database with reference genome information.
            String[] createSNPEffDatabase = {"java", "-jar", "snpEff.jar", "build",
                    "-gff3",
                    "-v", "reference"
            };
            CL.runCommand(createSNPEffDatabase,
                    this.targetDirectory + "/" + "snpEffDatabase.error",
                    this.targetDirectory + "/" + "snpEffDatabase.log",
                    this.targetDirectory.getAbsolutePath()
            );
            // Run snpEff annotation on each input .vcf file (multi-threaded).
            String[] runSNPEff = {"java", "-jar", "snpEff.jar",
                    "-v",
                    "-no-downstream", // Do not show DOWNSTREAM changes
                    "-no-intergenic", // Do not show INTERGENIC changes
                    "-no-intron", // Do not show INTRON changes
                    "-no-upstream", // Do not show UPSTREAM changes
                    "reference",
                    targetVcf.getAbsolutePath()
            };
            CL.runCommand(runSNPEff,
                    new File(this.targetDirectory + "/" + "snpEff.log").getAbsolutePath(),
                    new File(this.targetDirectory + "/" + "annotated_variants.vcf").getAbsolutePath(),
                    this.targetDirectory.getAbsolutePath()
            );
            // Transfer annotation results to musial storage.
            transferAnnotation();
        } finally {
            if (this.targetDirectory.isDirectory()) IO.deleteDirectory(this.targetDirectory);
        }
    }

    /**
     * Transfers all annotations to the MusialStorage instance.
     */
    private void transferAnnotation() {
        List<String> variantCallAnnotations;
        try {
            variantCallAnnotations = IO.readLinesFromFile(this.targetDirectory + "/" + "annotated_variants.vcf").stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());
            String pos;
            String ref;
            StringBuilder alt;
            String[] annotationFields;
            for (String variantCallAnnotation : variantCallAnnotations) {
                annotationFields = variantCallAnnotation.split("\t");
                pos = annotationFields[1];
                ref = annotationFields[3];
                alt = new StringBuilder(annotationFields[4]);
                if (ref.length() != 1 && alt.length() != 1)
                    continue; // TODO: Ambiguous calls not retrievable.
                if (ref.length() > alt.length() && alt.length() == 1)
                    alt.append(Constants.DELETION_OR_GAP_STRING.repeat(ref.length() - alt.length()));
                if (annotationFields[7].equals(".")) {
                    continue;
                } else {
                    annotationFields = annotationFields[7].replace("ANN=", "").split("\\|");
                }
                Iterator<String> affectedFeatureNames = musialStorage.getFeatureNameIterator(Integer.parseInt(pos));
                String affectedFeatureName;
                VariantInformation variantInformation;
                while (affectedFeatureNames.hasNext()) {
                    affectedFeatureName = affectedFeatureNames.next();
                    for (int i = 1; i < annotationFields.length; i++) {
                        if (i >= 12) {
                            break;
                        }
                        variantInformation = musialStorage.getFeature(affectedFeatureName).getNucleotideVariant(Integer.parseInt(pos), alt.toString());
                        variantInformation.addInfo(
                                Constants.SNP_EFF_INFO_PREFIX + Constants.SNP_EFF_INFO_KEYS.get(i - 1),
                                annotationFields[i]
                        );
                    }
                }
                alt.setLength(0);
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }
}
