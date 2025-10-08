package task;

import cli.CLISequence;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import model.Contig;
import model.Feature;
import model.Sample;
import model.Storage;
import op.AminoacidSequenceGenerator;
import op.NucleotideSequenceGenerator;
import op.SequenceGenerator;
import op.StorageFactory;
import util.Constants;
import util.Logging;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * The {@code ExecutorSequence} class is responsible for executing the {@code sequence} task to profile samples.
 */
public class ExecutorSequence {

    /**
     * Command-line interface for the sequence task.
     */
    private final CLISequence cli;

    /**
     * Storage instance to get genomic data.
     */
    private final Storage storage;

    /**
     * Set of sample identifiers to be processed. If empty, all samples in the storage will be processed.
     */
    private Set<String> sampleIdentifiers;

    /**
     * Set of features to be processed.
     */
    private final Set<Feature> features;

    /**
     * Set of contigs to be processed.
     */
    private final Set<Contig> contigs;

    /**
     * Map of contig names to their specified regions (start and end positions). Contigs do not have to be in this map to be processed.
     */
    private final Map<String, Tuple<Integer, Integer>> contigRegions;

    /**
     * Constructs an {@code ExecutorSequence} instance with the specified command-line interface.
     *
     * @param cli The command-line interface for the sequence task.
     * @throws IOException     If an I/O error occurs while loading the storage.
     * @throws MusialException If a validation error occurs with the input parameters.
     */
    public ExecutorSequence(CLISequence cli) throws IOException, MusialException {
        this.cli = cli;
        int numberOfLoci = this.cli.loci.size();
        this.features = new HashSet<>(numberOfLoci);
        this.contigs = new HashSet<>(numberOfLoci);
        this.contigRegions = new HashMap<>(numberOfLoci);
        Logging.logInfo("Load storage.");
        this.storage = StorageFactory.fromPath(cli.input);
        Logging.logDone("");
        Logging.logInfo("Validate parameters.");
        validate();
        Logging.logDone("");
    }

    /**
     * Executes the sequence generation task based on the specified content type.
     * <p>
     * This method determines the type of sequence data to generate (amino acid or nucleotide) based on the {@code cli.content} value. It
     * then processes the features and, if applicable, the contigs to generate the required sequence data. The sequence data is written to
     * output files as specified in the command-line arguments.
     *
     * @throws IOException     If an I/O error occurs during sequence generation or file writing.
     * @throws MusialException If an error occurs during sequence generation or validation.
     */
    public void run() throws IOException, MusialException {
        Logging.logInfo("Generate sequence data.");
        if (cli.content.equals(CLISequence.Content.AMINOACID)) {
            processFeatures(true);
        } else {
            processFeatures(false);
            processContigs();
        }
        Logging.logDone("");
    }

    /**
     * Processes features and generates sequence data for each feature.
     * <p>
     * This method iterates over all features specified in the {@code features} set. For each feature, it retrieves the associated contig
     * and initializes a {@link SequenceGenerator} based on the content type (amino acid or nucleotide). The sequence data for each sample
     * is written to a file using a {@code BufferedWriter}. If the merge option is enabled, merged sequences are written; otherwise,
     * individual sequences are written.
     *
     * @param isAminoAcid A boolean indicating whether the content type is amino acid (`true`) or nucleotide (`false`).
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processFeatures(boolean isAminoAcid) throws IOException, MusialException {
        for (Feature feature : features) {
            // Open a writer for the output file corresponding to the feature
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(cli.outputGenerator.apply(feature.name), cli.append))) {
                // Retrieve the contig associated with the feature
                Contig contig = storage.getContig(feature.contig);

                // Create a sequence generator for the feature
                SequenceGenerator generator = createFeatureSequenceGenerator(isAminoAcid, contig, feature);

                // Write sequences based on the merge option
                if (cli.merge) {
                    writeMergedSequences(writer, generator, feature, isAminoAcid);
                } else {
                    writeIndividualSequences(writer, generator, feature);
                }
            }
        }
    }

    /**
     * Processes contigs and generates sequence data for each contig.
     * <p>
     * This method iterates over all contigs specified in the {@code contigs} set. For each contig, it determines the locus identifier based
     * on whether the contig has a specified region in the {@code contigRegions} map. It then initializes a
     * {@link NucleotideSequenceGenerator} for the contig, either for the entire contig or for the specified region. The sequence data for
     * each sample is written to a file using a {@code BufferedWriter}.
     *
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processContigs() throws IOException, MusialException {
        for (Contig contig : contigs) {
            // Determine the locus identifier based on whether the contig has a specified region
            String locusId = contigRegions.containsKey(contig._id)
                    ? "%sg%d_%d".formatted(contig._id, contigRegions.get(contig._id).a, contigRegions.get(contig._id).b)
                    : contig._id;

            // Initialize the sequence generator for the contig or the specified region
            SequenceGenerator generator = contigRegions.containsKey(contig._id)
                    ? new NucleotideSequenceGenerator(storage, contig, contigRegions.get(contig._id).a,
                    contigRegions.get(contig._id).b, !cli.variable, cli.align, sampleIdentifiers)
                    : new NucleotideSequenceGenerator(storage, contig, !cli.variable, cli.align, sampleIdentifiers);

            // Write the sequence data for each sample to the output file
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(cli.outputGenerator.apply(locusId), cli.append))) {
                for (Sample sample : storage.getSamples()) {
                    writer.write(">lcl|%s|%s%n%s%n".formatted(sample._id, locusId, generator.getSequence(sample._id)));
                }
            }
        }
    }

    /**
     * Creates a sequence generator based on the specified content type (amino acid or nucleotide).
     * <p>
     * This method initializes and returns an appropriate sequence generator for the given feature and contig. If the content type is amino
     * acid, an {@link AminoacidSequenceGenerator} is created. Otherwise, a {@link NucleotideSequenceGenerator} is created. The generator is
     * configured based on the command-line options provided in the {@code cli} object.
     *
     * @param isAminoAcid A boolean indicating whether the content type is amino acid (`true`) or nucleotide (`false`).
     * @param contig      The contig associated with the feature for which the generator is being created.
     * @param feature     The feature for which the generator is being created.
     * @return A {@link SequenceGenerator} instance configured for the specified content type, contig, and feature.
     * @throws IOException     If an I/O error occurs during generator initialization.
     * @throws MusialException If an error occurs during generator creation or validation.
     */
    private SequenceGenerator createFeatureSequenceGenerator(boolean isAminoAcid, Contig contig, Feature feature) throws IOException,
            MusialException {
        return isAminoAcid
                ? new AminoacidSequenceGenerator(storage, contig, feature, !cli.variable, cli.align, sampleIdentifiers)
                : new NucleotideSequenceGenerator(storage, contig, feature, !cli.variable, cli.align, sampleIdentifiers);
    }

    /**
     * Writes merged sequences for a given feature to the specified writer.
     * <p>
     * This method iterates over all samples in the storage and generates sequence data for each sample. It ensures that duplicate sequences
     * are not written by maintaining a set of observed identifiers. The sequence data includes the identifier, feature name, and the
     * generated sequence. The identifier is determined based on whether the content is amino acid or nucleotide.
     *
     * @param writer      The `BufferedWriter` used to write the sequence data to a file.
     * @param generator   The `SequenceGenerator` used to generate sequences for the feature.
     * @param feature     The feature for which sequences are being written.
     * @param isAminoAcid A boolean indicating whether the content is amino acid (`true`) or nucleotide (`false`).
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void writeMergedSequences(BufferedWriter writer, SequenceGenerator generator, Feature feature, boolean isAminoAcid) throws IOException, MusialException {
        Set<String> observed = new HashSet<>();
        for (Sample sample : storage.getSamples()) {
            String id = isAminoAcid
                    ? getProteoformId(feature, sample)
                    : sample.getRelatedAllele(feature._id);
            if (observed.add(id)) {
                writer.write(">lcl|%s|%s%n%s%n".formatted(id, feature.name, generator.getSequence(sample._id)));
            }
        }
    }

    /**
     * Writes individual sequences per sample for a given feature with the specified writer.
     * <p>
     * This method iterates over all samples in the storage and writes the sequence data for each sample to the provided `BufferedWriter`.
     * The sequence data is generated using the specified `SequenceGenerator` and includes the sample ID, feature name, and the generated
     * sequence.
     *
     * @param writer    The `BufferedWriter` used to write the sequence data to a file.
     * @param generator The `SequenceGenerator` used to generate sequences for the feature.
     * @param feature   The feature for which sequences are being written.
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void writeIndividualSequences(BufferedWriter writer, SequenceGenerator generator, Feature feature) throws IOException,
            MusialException {
        for (Sample sample : storage.getSamples()) {
            writer.write(">lcl|%s|%s%n%s%n".formatted(sample._id, feature.name, generator.getSequence(sample._id)));
        }
    }

    /**
     * Retrieves the proteoform identifier for a given feature and sample.
     * <p>
     * This method determines the proteoform ID based on the allele associated with the sample for the specified feature. If the allele is
     * the reference allele, the proteoform ID is set to a constant representing a synonymous change. Otherwise, the proteoform ID is
     * derived from the related proteoform of the allele.
     *
     * @param feature The feature for which the proteoform ID is being retrieved.
     * @param sample  The sample associated with the feature.
     * @return The proteoform identifier as a string.
     */
    private String getProteoformId(Feature feature, Sample sample) {
        String alleleId = sample.getRelatedAllele(feature._id);
        return alleleId.equals(Constants.REFERENCE)
                ? Constants.SYNONYMOUS
                : feature.getAllele(alleleId).getRelatedProteoform();
    }

    /**
     * Validates the input parameters for the sequence task.
     * <p>
     * This method performs a series of validation checks to ensure that the input parameters provided in the command-line arguments are
     * correct and consistent. It validates the following:
     * <ul>
     *   <li>Samples: Ensures that the specified samples exist in the storage.</li>
     *   <li>Loci: Ensures that the specified loci are valid and properly formatted.</li>
     *   <li>Configuration: Ensures that the overall configuration settings are valid.</li>
     * </ul>
     * If any of these validations fail, a {@link MusialException} is thrown with an appropriate error message.
     *
     * @throws MusialException If any of the validation checks fail.
     */
    private void validate() throws MusialException {
        validateSamples();
        validateLoci();
        validateConfiguration();
    }

    /**
     * Validates the configuration settings provided in the command-line arguments.
     * <p>
     * This method performs several checks to ensure the configuration is valid:
     * <ul>
     *   <li>Ensures that loci without sequence information can only proceed if the variable option is enabled.</li>
     *   <li>Validates that the merge option requires the storage to include typing information.</li>
     *   <li>Checks that feature-based sequence generation requires typing information in the storage.</li>
     *   <li>Ensures that amino acid sequence generation has at least one feature specified.</li>
     * </ul>
     * If any of these conditions are not met, a {@link MusialException} is thrown with an appropriate error message.
     *
     * @throws MusialException If the configuration is invalid based on the checks performed.
     */
    private void validateConfiguration() throws MusialException {
        // Check if any contig is missing sequence information and the variable option is not enabled
        if (!cli.variable && contigs.stream().anyMatch(contig -> !contig.hasSequence())) {
            throw new MusialException("Loci without sequence information require the -v/--variable option to proceed. Please enable this " +
                    "option and try again.");
        }

        // Validate the merge option and typing information
        if (cli.merge && storage.parameters.skipTyping()) {
            throw new MusialException("The -m/--merge option requires the storage to include typing information. Please recreate the " +
                    "storage with typing information and try again.");
        }

        // Ensure merge option is only used with feature-based sequence generation.
        if (cli.merge && !contigs.isEmpty()) {
            throw new MusialException("The -m/--merge option is only applicable for feature-based sequence generation. Please remove " +
                    "contigs from the loci and try again.");
        }

        // Validate feature-based sequence generation and typing information
        if (!features.isEmpty() && storage.parameters.skipTyping()) {
            throw new MusialException("Sequence generation of features requires the storage to include typing information. Please " +
                    "recreate the storage with typing information and try again.");
        }

        // Ensure amino acid sequence generation has at least one feature specified
        if (features.isEmpty() && cli.content.equals(CLISequence.Content.AMINOACID)) {
            throw new MusialException("Amino acid sequence generation requires at least one feature to be specified.");
        }
    }

    /**
     * Validates the samples provided in the command-line arguments.
     * <p>
     * This method filters the samples specified in the command-line input to include only those that exist in the storage. If no valid
     * samples are found, an exception is thrown. Additionally, it logs a warning if some specified samples do not exist in the storage and
     * are removed.
     * <p>
     * If no samples are specified in the command-line input, all samples from the storage are used; this is indicated by an empty set.
     *
     * @throws MusialException If none of the specified samples exist in the storage.
     */
    private void validateSamples() throws MusialException {
        // Filter samples based on their existence in the storage or initialize with empty set.
        this.sampleIdentifiers = cli.samples.isEmpty() ? cli.samples :
                cli.samples.stream().filter(storage::hasSample).collect(Collectors.toSet());

        // Throw an exception if no valid samples exist.
        if (!cli.samples.isEmpty() && this.sampleIdentifiers.isEmpty()) {
            throw new MusialException("None of the specified samples exists in the storage.");
        }

        // Log a warning if some samples were removed
        int removedSamples = cli.samples.size() - this.sampleIdentifiers.size();
        if (removedSamples > 0) {
            Logging.logWarning("%d of the specified samples do not exist in the storage and were removed.".formatted(removedSamples));
        }
    }

    /**
     * Validates the loci provided in the command-line arguments.
     * <p>
     * This method iterates through each locus specified in the command-line input and validates it using the {@code validateLocus} method.
     * After validation, it ensures that at least one valid locus (feature or contig) has been identified. If no valid loci are found, an
     * exception is thrown.
     *
     * @throws MusialException If no valid loci are specified or if an error occurs during validation.
     */
    private void validateLoci() throws MusialException {
        // Validate each locus in the command-line input
        for (String locus : cli.loci) {
            validateLocus(locus);
        }

        // Ensure at least one valid locus is specified
        if (features.isEmpty() && contigs.isEmpty()) {
            throw new MusialException("No valid loci were specified.");
        }
    }

    /**
     * Validates a given locus and categorizes it as a feature, contig, or contig region.
     * <p>
     * This method checks if the provided locus corresponds to a feature, a contig, or a contig region. If the locus is a feature, it is
     * added to the {@code features} set. If it is a contig, it is added to the {@code contigs} set. If it is a contig region (in the format
     * CHROM:START-END), it is validated and added to the {@code contigRegions} map and {@code contigs} set. If the locus does not match any
     * of these categories, an exception is thrown.
     *
     * @param locus The locus to validate, which can be a feature name, contig name, or contig region.
     * @throws MusialException If the locus is invalid, does not exist, or is not in the expected format.
     */
    private void validateLocus(String locus) throws MusialException {
        // Check if the locus is a feature
        Feature feature = storage.hasFeature(locus)
                ? storage.getFeature(locus)
                : storage.getFeatures().stream().filter(f -> f.name.equals(locus)).findFirst().orElse(null);
        if (Objects.nonNull(feature)) {
            if (cli.content.equals(CLISequence.Content.AMINOACID) && !feature.isCoding()) {
                throw new MusialException("Feature '%s' is not a coding feature and cannot be used for amino acid sequence generation.".formatted(locus));
            }
            features.add(feature);
            return;
        }

        // Split the locus into parts for contig validation
        String[] parts = locus.split(":");
        if (parts.length == 1 && storage.hasContig(parts[0])) {
            if (cli.content.equals(CLISequence.Content.AMINOACID)) {
                throw new MusialException(("Contigs ('%s') cannot be used for amino acid sequence generation, please specify only " +
                        "coding features").formatted(locus));
            }
            contigs.add(storage.getContig(parts[0]));
            return;
        }

        if (parts.length == 2) {
            String contigIdentifier = parts[0];
            if (storage.hasContig(contigIdentifier)) {
                try {
                    String[] regionParts = parts[1].split("-");
                    if (regionParts.length != 2) {
                        throw new MusialException("Locus '%s' is not valid, expected CHROM:START-END.".formatted(locus));
                    }
                    int start = Integer.parseInt(regionParts[0]);
                    int end = Integer.parseInt(regionParts[1]);
                    if (start < 1 || end < 1 || start > end) {
                        throw new MusialException("Locus '%s' is not valid, expected CHROM:START-END with 1 < START < END.".formatted(locus));
                    }
                    if (cli.content.equals(CLISequence.Content.AMINOACID)) {
                        throw new MusialException(("Contigs ('%s') cannot be used for amino acid sequence generation, please specify only" +
                                " " +
                                "coding features").formatted(locus));
                    }
                    contigs.add(storage.getContig(contigIdentifier));
                    contigRegions.put(contigIdentifier, new Tuple<>(start, end));
                    return;
                } catch (NumberFormatException e) {
                    throw new MusialException("Locus '%s' is not valid, expected CHROM:START-END.".formatted(locus));
                }
            }
        }

        throw new MusialException("Locus '%s' is not valid, it is neither a feature nor a contig.".formatted(locus));
    }


}
