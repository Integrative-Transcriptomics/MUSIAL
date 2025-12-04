package task;

import cli.CLISequence;
import com.google.common.base.Splitter;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import model.Contig;
import model.Feature;
import model.Storage;
import op.AminoAcidSequenceGenerator;
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
    private Set<Feature> features;

    /**
     * Set of contigs to be processed.
     */
    private Set<Contig> contigs;

    /**
     * Map of contig names to their specified regions (start and end positions). Contigs do not have to be in this map to be processed.
     */
    private Map<String, Tuple<Integer, Integer>> contigRegions;

    /**
     * Constructs an {@code ExecutorSequence} instance with the specified command-line interface.
     *
     * @param cli The command-line interface for the sequence task.
     * @throws IOException     If an I/O error occurs while loading the storage.
     * @throws MusialException If a validation error occurs with the input parameters.
     */
    public ExecutorSequence(CLISequence cli) throws IOException, MusialException {
        this.cli = cli;
        Logging.logInfo("Load storage.");
        this.storage = StorageFactory.fromPath(cli.input);
        Logging.logDone("");
        Logging.logInfo("Validate parameters.");
        validateLoci();
        validateConfiguration();
        validateSamples();
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
     * This method determines how to process the features based on the `split` parameter provided in the command-line interface. The `split`
     * parameter specifies the output file structure:
     * <ul>
     *   <li><b>FEATURE</b>: One file per feature/contig containing all sample sequences.</li>
     *   <li><b>SAMPLE</b>: One file per sample containing all feature sequences.</li>
     *   <li><b>NONE</b>: All sequences are written into a single file.</li>
     *   <li><b>BOTH</b>: Each sequence is written into its own file.</li>
     * </ul>
     * Depending on the `split` mode, the appropriate processing method is invoked.
     *
     * @param asAminoAcid A boolean indicating whether the content type is amino acid (`true`) or nucleotide (`false`).
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processFeatures(boolean asAminoAcid) throws IOException, MusialException {
        switch (cli.split) {
            case FEATURE -> processFeaturesByLocus(asAminoAcid); // Process features into one file per feature.
            case SAMPLE -> processFeaturesBySample(asAminoAcid);  // Process features into one file per sample.
            case NONE -> processFeaturesInOneFile(asAminoAcid);   // Process all features into a single file.
            case BOTH -> processFeaturesIndividually(asAminoAcid); // Process each sequence into its own file.
        }
    }

    /**
     * Creates a sequence generator based on the specified content type (amino acid or nucleotide).
     * <p>
     * This method initializes and returns an appropriate sequence generator for the given feature and contig. If the content type is amino
     * acid, an {@link AminoAcidSequenceGenerator} is created. Otherwise, a {@link NucleotideSequenceGenerator} is created. The generator is
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
                ? new AminoAcidSequenceGenerator(storage, contig, feature, !cli.variable, cli.align, sampleIdentifiers)
                : new NucleotideSequenceGenerator(storage, contig, feature, !cli.variable, cli.align, sampleIdentifiers);
    }

    /**
     * Writes per sample feature sequences into files organized by features.
     * <p>
     * Each file will contain sequences of all samples for a specific feature. The sequence content is handled automatically by the
     * {@link SequenceGenerator} based on the specified content type. Merging of sequences is managed by the
     * {@link #processSamplesWithWriter} method.
     *
     * @param isAminoAcid A boolean indicating whether the content type is amino acid (`true`) or nucleotide (`false`).
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processFeaturesByLocus(boolean isAminoAcid) throws IOException, MusialException {
        for (Feature feature : features) {
            // Retrieve the contig associated with the feature
            Contig contig = storage.getContig(feature.contig);
            // Create a sequence generator for the feature
            SequenceGenerator generator = createFeatureSequenceGenerator(isAminoAcid, contig, feature);
            // Open a writer for the feature-specific output file
            try (BufferedWriter writer =
                         new BufferedWriter(new FileWriter(cli.outputGenerator.apply(generator.getName(true)), false))) {
                // Process and write sequences for all samples associated with the feature
                processSamplesWithWriter(isAminoAcid, writer, feature, generator);
            }
        }
    }

    /**
     * Writes per sample feature sequences into files organized by samples.
     * <p>
     * Each file will contain sequences of all features for a specific sample. The sequence content is handled automatically by the
     * {@link SequenceGenerator} based on the specified content type. Merging of sequences is not implemented, as it is not allowed by the
     * {@code -m/--merge} option when using sample-based splitting, see {@link CLISequence}.
     *
     * @param isAminoAcid A boolean indicating whether the content type is amino acid (`true`) or nucleotide (`false`).
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processFeaturesBySample(boolean isAminoAcid) throws IOException, MusialException {
        for (Feature feature : features) {
            // Retrieve the contig associated with the feature
            Contig contig = storage.getContig(feature.contig);
            // Create a sequence generator for the feature
            SequenceGenerator generator = createFeatureSequenceGenerator(isAminoAcid, contig, feature);
            for (String sampleIdentifier : sampleIdentifiers) {
                // Open a writer for the sample-specific output file
                try (BufferedWriter writer =
                             new BufferedWriter(new FileWriter(cli.outputGenerator.apply("%s-sequences".formatted(sampleIdentifier)),
                                     true))) {
                    // Write the sequence for the sample in FASTA format
                    writer.write(">lcl|%s|%s%n%s%n".formatted(sampleIdentifier, generator.getName(false),
                            formatSequence(generator.getSequence(sampleIdentifier))));
                }
            }
        }
    }

    /**
     * Writes per sample feature sequences into a single output file.
     * <p>
     * The sequence content is handled automatically by the {@link SequenceGenerator} based on the specified content type. Merging of
     * sequences is managed by the {@link #processSamplesWithWriter} method.
     *
     * @param isAminoAcid A boolean indicating whether the content type is amino acid (`true`) or nucleotide (`false`).
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processFeaturesInOneFile(boolean isAminoAcid) throws IOException, MusialException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(cli.outputGenerator.apply("musial-sequences"), true))) {
            for (Feature feature : features) {
                Contig contig = storage.getContig(feature.contig);
                SequenceGenerator generator = createFeatureSequenceGenerator(isAminoAcid, contig, feature);
                processSamplesWithWriter(isAminoAcid, writer, feature, generator);
            }
        }
    }

    /**
     * Writes per sample feature sequences into individual output files.
     * <p>
     * The sequence content is handled automatically by the {@link SequenceGenerator} based on the specified content type. Merging of
     * sequences is managed by the {@link #processSamplesWithWriter} method.
     *
     * @param isAminoAcid A boolean indicating whether the content type is amino acid (`true`) or nucleotide (`false`).
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processFeaturesIndividually(boolean isAminoAcid) throws IOException, MusialException {
        for (Feature feature : features) {
            // Retrieve the contig associated with the feature
            Contig contig = storage.getContig(feature.contig);
            // Create a sequence generator for the feature
            SequenceGenerator generator = createFeatureSequenceGenerator(isAminoAcid, contig, feature);

            if (cli.merge) {
                // Set to track observed merged identifiers to avoid duplicate entries
                Set<String> observed = new HashSet<>(feature.getAlleleCount());
                for (String sampleIdentifier : sampleIdentifiers) {
                    // Determine the merged identifier based on the content type
                    String mergedIdentifier = isAminoAcid
                            ? getProteoformId(feature, sampleIdentifier)
                            : getAlleleIdentifier(feature, sampleIdentifier);
                    try (BufferedWriter writer =
                                 new BufferedWriter(new FileWriter(cli.outputGenerator.apply("%s-%s".formatted(mergedIdentifier,
                                         generator.getName(true))), false))) {
                        // Write the sequence if the merged identifier has not been observed
                        if (observed.add(mergedIdentifier)) {
                            writer.write(">lcl|%s|%s%n%s%n".formatted(mergedIdentifier, generator.getName(false),
                                    formatSequence(generator.getSequence(sampleIdentifier))));
                        }
                    }
                }
            } else {
                for (String sampleIdentifier : sampleIdentifiers) {
                    try (BufferedWriter writer =
                                 new BufferedWriter(new FileWriter(cli.outputGenerator.apply("%s-%s".formatted(sampleIdentifier,
                                         generator.getName(true))), false))) {
                        // Write the sequence for each sample
                        writer.write(">lcl|%s|%s%n%s%n".formatted(sampleIdentifier, generator.getName(false),
                                formatSequence(generator.getSequence(sampleIdentifier))));
                    }
                }
            }
        }
    }

    /**
     * Processes sample sequences and writes them to the provided writer.
     * <p>
     * This method handles the generation and writing of sequences for each sample based on the `merge` option:
     * <ul>
     *   <li>If `merge` is enabled, sequences are grouped by a merged identifier (proteoform ID for amino acid content or allele ID for
     *   nucleotide content).</li>
     *   <li>If `merge` is disabled, sequences are written individually for each sample.</li>
     * </ul>
     *
     * @param isAminoAcid A boolean indicating whether the content type is amino acid (`true`) or nucleotide (`false`).
     * @param writer      The {@link BufferedWriter} used to write the sequences to the output file.
     * @param feature     The feature being processed, used to retrieve allele or proteoform information.
     * @param generator   The {@link SequenceGenerator} responsible for generating sequences for the feature.
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processSamplesWithWriter(boolean isAminoAcid, BufferedWriter writer, Feature feature, SequenceGenerator generator) throws IOException, MusialException {
        if (cli.merge) {
            // Set to track observed merged identifiers to avoid duplicate entries
            Set<String> observed = new HashSet<>(feature.getAlleleCount());
            for (String sampleIdentifier : sampleIdentifiers) {
                // Determine the merged identifier based on the content type
                String mergedIdentifier = isAminoAcid
                        ? getProteoformId(feature, sampleIdentifier)
                        : getAlleleIdentifier(feature, sampleIdentifier);
                // Write the sequence if the merged identifier has not been observed
                if (observed.add(mergedIdentifier)) {
                    writer.write(">lcl|%s|%s%n%s%n".formatted(mergedIdentifier, generator.getName(false),
                            formatSequence(generator.getSequence(sampleIdentifier))));
                }
            }
        } else {
            // Write sequences individually for each sample
            for (String sampleIdentifier : sampleIdentifiers) {
                writer.write(">lcl|%s|%s%n%s%n".formatted(sampleIdentifier, generator.getName(false),
                        formatSequence(generator.getSequence(sampleIdentifier))));
            }
        }
    }

    /**
     * Processes contigs and generates sequence data for each contig.
     * <p>
     * This method determines how to process the contigs based on the `split` parameter provided in the command-line interface. The `split`
     * parameter specifies the output file structure:
     * <ul>
     *   <li><b>FEATURE</b>: One file per contig containing all sample sequences.</li>
     *   <li><b>SAMPLE</b>: One file per sample containing all contig sequences.</li>
     *   <li><b>NONE</b>: All sequences are written into a single file.</li>
     *   <li><b>BOTH</b>: Each sequence is written into its own file.</li>
     * </ul>
     * Depending on the `split` mode, the appropriate processing method is invoked.
     *
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processContigs() throws IOException, MusialException {
        switch (cli.split) {
            case FEATURE -> processContigsByLocus(); // Process contigs into one file per contig.
            case SAMPLE -> processContigsBySample();  // Process contigs into one file per sample.
            case NONE -> processContigsInOneFile();   // Process all contigs into a single file.
            case BOTH -> processContigsIndividually(); // Process each sequence into its own file.
        }
    }

    /**
     * Creates a sequence generator for a given contig.
     * <p>
     * This method initializes and returns a {@link NucleotideSequenceGenerator} for the specified contig. If the contig is associated with
     * a specific region (start and end positions) in the {@code contigRegions} map, the generator is configured to process only that
     * region. Otherwise, the generator processes the entire contig.
     * <p>
     * The generator is further configured based on the command-line options provided in the {@code cli} object.
     *
     * @param contig The contig for which the sequence generator is being created.
     * @return A {@link NucleotideSequenceGenerator} instance configured for the specified contig.
     * @throws IOException     If an I/O error occurs during generator initialization.
     * @throws MusialException If an error occurs during generator creation or validation.
     */
    private SequenceGenerator createContigSequenceGenerator(Contig contig) throws IOException, MusialException {
        return contigRegions.containsKey(contig._id)
                ? new NucleotideSequenceGenerator(storage, contig, contigRegions.get(contig._id).a, contigRegions.get(contig._id).b,
                !cli.variable, cli.align, sampleIdentifiers)
                : new NucleotideSequenceGenerator(storage, contig, !cli.variable, cli.align, sampleIdentifiers);
    }

    /**
     * Writes per sample contig sequences into files organized by contigs.
     * <p>
     * Each file will contain sequences for all samples for a specific contig.
     *
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processContigsByLocus() throws IOException, MusialException {
        for (Contig contig : contigs) {
            SequenceGenerator generator = createContigSequenceGenerator(contig);
            int bufferSize = estimateBufferSize(generator);
            try (BufferedWriter writer =
                         new BufferedWriter(new FileWriter(cli.outputGenerator.apply(generator.getName(true)), false),
                                 bufferSize)) {
                for (String sampleIdentifier : sampleIdentifiers) {
                    writer.write(">lcl|%s|%s%n%s%n".formatted(sampleIdentifier, generator.getName(false),
                            formatSequence(generator.getSequence(sampleIdentifier))));
                }
            }
        }
    }

    /**
     * Writes per sample contig sequences into files organized by samples.
     * <p>
     * Each file will contain sequences for all contigs for a specific sample.
     *
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processContigsBySample() throws IOException, MusialException {
        for (Contig contig : contigs) {
            SequenceGenerator generator = createContigSequenceGenerator(contig);
            int bufferSize = estimateBufferSize(generator);
            for (String sampleIdentifier : sampleIdentifiers) {
                try (BufferedWriter writer =
                             new BufferedWriter(new FileWriter(cli.outputGenerator.apply("%s-sequences".formatted(sampleIdentifier)),
                                     true), bufferSize)) {
                    writer.write(">lcl|%s|%s%n%s%n".formatted(sampleIdentifier, generator.getName(false),
                            formatSequence(generator.getSequence(sampleIdentifier))));
                }
            }
        }
    }

    /**
     * Writes per sample contig sequences into a single output file.
     *
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processContigsInOneFile() throws IOException, MusialException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(cli.outputGenerator.apply("musial-sequences"), true), 1048576)) {
            for (Contig contig : contigs) {
                SequenceGenerator generator = createContigSequenceGenerator(contig);
                for (String sampleIdentifier : sampleIdentifiers) {
                    writer.write(">lcl|%s|%s%n%s%n".formatted(sampleIdentifier, generator.getName(false),
                            formatSequence(generator.getSequence(sampleIdentifier))));
                }
            }
        }
    }

    /**
     * Writes per sample contig sequences into individual output files.
     *
     * @throws IOException     If an I/O error occurs while writing to the file.
     * @throws MusialException If an error occurs during sequence generation.
     */
    private void processContigsIndividually() throws IOException, MusialException {
        for (Contig contig : contigs) {
            SequenceGenerator generator = createContigSequenceGenerator(contig);
            int bufferSize = estimateBufferSize(generator);
            for (String sampleIdentifier : sampleIdentifiers) {
                try (BufferedWriter writer =
                             new BufferedWriter(new FileWriter(cli.outputGenerator.apply("%s-%s".formatted(sampleIdentifier,
                                     generator.getName(true))), false), bufferSize)) {
                    writer.write(">lcl|%s|%s%n%s%n".formatted(sampleIdentifier, generator.getName(false),
                            formatSequence(generator.getSequence(sampleIdentifier))));
                }
            }
        }
    }

    /**
     * Retrieves the allele identifier for a given feature and sample.
     * <p>
     * This method fetches the sample from the storage using the provided sample identifier and retrieves the allele related to the
     * specified feature.
     *
     * @param feature          The feature for which the allele identifier is being retrieved.
     * @param sampleIdentifier The identifier of the sample associated with the feature.
     * @return A {@link String} representing the allele identifier related to the feature for the given sample.
     */
    private String getAlleleIdentifier(Feature feature, String sampleIdentifier) {
        return storage.getSample(sampleIdentifier).getRelatedAllele(feature._id);
    }

    /**
     * Retrieves the proteoform identifier for a given feature and sample.
     * <p>
     * This method determines the proteoform ID based on the allele associated with the sample for the specified feature.
     * <ul>
     *   <li>If the allele is the reference allele, the proteoform ID is set to a constant representing a synonymous change.</li>
     *   <li>If the allele is not the reference allele, the proteoform ID is derived from the related proteoform of the allele.</li>
     * </ul>
     *
     * @param feature          The feature for which the proteoform ID is being retrieved.
     * @param sampleIdentifier The identifier of the sample associated with the feature.
     * @return The proteoform identifier as a string.
     */
    private String getProteoformId(Feature feature, String sampleIdentifier) {
        // Retrieve the allele identifier for the given feature and sample
        String alleleIdentifier = getAlleleIdentifier(feature, sampleIdentifier);

        // Return the proteoform ID based on whether the allele is the reference allele
        return alleleIdentifier.equals(Constants.REFERENCE)
                ? Constants.SYNONYMOUS
                : feature.getAllele(alleleIdentifier).getRelatedProteoform();
    }

    private int estimateBufferSize(SequenceGenerator generator) {
        int opt = (int) Math.pow(2, Math.ceil(Math.log(generator.getSize()) / Math.log(2)));
        return Math.clamp(opt, 8192, 1048576);
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
        if (this.cli.loci.isEmpty()) {
            // If no loci are specified, load all features or contigs based on the content type.
            loadLoci();
        } else {
            // Initialize sets and map with an estimated size based on the number of loci provided.
            int numberOfLoci = this.cli.loci.size();
            this.features = new HashSet<>(numberOfLoci);
            this.contigs = new HashSet<>(numberOfLoci);
            this.contigRegions = new HashMap<>(numberOfLoci);

            // Validate each locus in the command-line input
            for (String locus : cli.loci) {
                validateLocus(locus);
            }

            // Ensure at least one valid locus is specified
            if (features.isEmpty() && contigs.isEmpty()) {
                throw new MusialException("No valid loci were specified.");
            }
        }
    }

    /**
     * Loads loci (features or contigs) based on the content type specified in the command-line interface.
     * <p>
     * This method initializes the `features` and `contigs` sets based on the content type:
     * <ul>
     *   <li>If the content type is amino acid, it filters and loads only coding features.</li>
     *   <li>If the content type is nucleotide, it loads all features.</li>
     * </ul>
     * If no features are found:
     * <ul>
     *   <li>For amino acid content, an exception is thrown indicating no coding features are available.</li>
     *   <li>For nucleotide content, all contigs are loaded instead.</li>
     * </ul>
     * The `contigRegions` map is initialized as empty in any case.
     *
     * @throws MusialException If no coding features are found for amino acid content.
     */
    private void loadLoci() throws MusialException {
        // Check if the content type is amino acid
        if (this.cli.content.equals(CLISequence.Content.AMINOACID)) {
            // Load only coding features for amino acid content
            this.features = this.storage.getFeatures().stream().filter(Feature::isCoding).collect(Collectors.toSet());
        } else {
            // Load all features for nucleotide content
            this.features = new HashSet<>(this.storage.getFeatures());
        }

        // If no features are found, handle based on content type
        if (this.features.isEmpty()) {
            if (this.cli.content.equals(CLISequence.Content.AMINOACID)) {
                // Throw an exception if no coding features are found for amino acid content
                throw new MusialException("The storage does not contain any coding features required for amino acid sequence generation.");
            } else {
                // Load all contigs for nucleotide content
                this.contigs = new HashSet<>(this.storage.getContigs());
                this.contigRegions = Collections.emptyMap();
            }
        } else {
            // Initialize contigs and contigRegions as empty if features are loaded
            this.contigs = Collections.emptySet();
            this.contigRegions = Collections.emptyMap();
        }
    }

    /**
     * Formats a given sequence string into lines of a specified length.
     * <p>
     * This method splits the input sequence into lines of 80 characters each, making it suitable for formats that require line breaks, such
     * as FASTA. The formatted sequence is returned as a single string with newline characters separating the lines.
     *
     * @param sequence The input sequence string to be formatted.
     * @return A formatted sequence string with lines of 80 characters each.
     */
    private String formatSequence(String sequence) {
        return String.join("\n", Splitter.fixedLength(80).splitToList(sequence));
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
        if (cli.samples.isEmpty()) {
            // If no samples are specified, use all samples from the storage.
            this.sampleIdentifiers = storage.getSamples().stream().map(s -> s._id).collect(Collectors.toSet());
        } else {
            // Filter samples based on their existence in the storage.
            this.sampleIdentifiers = cli.samples.stream().filter(storage::hasSample).collect(Collectors.toSet());
        }

        // Throw an exception if no valid samples exist.
        if (this.sampleIdentifiers.isEmpty()) {
            throw new MusialException("None of the specified samples exists in the storage.");
        }

        // Log a warning if some samples were removed
        int removedSamples = cli.samples.size() - this.sampleIdentifiers.size();
        if (removedSamples > 0) {
            Logging.logWarning("%d of the specified samples do not exist in the storage and were removed.".formatted(removedSamples));
        }
    }

}
