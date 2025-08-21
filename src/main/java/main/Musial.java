package main;

import com.google.common.base.Splitter;
import model.*;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import org.apache.commons.io.FileUtils;
import utility.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.logging.Level;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Main class of MUSIAL (MUlti Sample varIant AnaLysis).
 * <p>
 * MUSIAL is a Java command-line tool designed to analyze and summarize single nucleotide variants (SNVs) and insertions/deletions (indels)
 * across multiple prokaryotic samples. The software aggregates and analyzes variant calls from multiple samples of a prokaryotic species
 * and provides an interface to generate comprehensive statistics and alignments at the genome, gene and protein level. MUSIAL enables a
 * comprehensive assessment of variability within a species at the genome, gene and protein level, providing insights into, for example,
 * conserved and variable regions, diversity at the gene level and common proteoforms among samples.
 */
public final class Musial {

    /**
     * Name of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String name = "";

    /**
     * Version of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String version = "";

    /**
     * Author contact of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String contact = "";

    /**
     * License information of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String license = "";

    /**
     * Specifies the task to execute.
     */
    public static Task task;

    /**
     * Output directory used for generated MUSIAL storage files.
     */
    public static File outputDirectory;

    /**
     * File extension used for generated MUSIAL storage files.
     * <p>
     * This constant specifies the default file extension for storage files
     * created or used by the MUSIAL application. It is marked as `transient`
     * to indicate that it should not be serialized as part of the class state.
     * <p>
     * The extension should either be `.json` or `.json.gz` depending on the
     * compression method used. For production, use `.json.gz` for compressed files.
     */
    public static final String outputExtension = ".json.gz";

    /**
     * Start time of the program.
     */
    public static long startTime;

    /**
     * {@link Enum} specifying MUSIAL tasks.
     */
    public enum Task {
        /**
         * Task to build a MUSIAL storage file.
         */
        BUILD,
        /**
         * Task to add (sample) data from VCF files to a MUSIAL storage file.
         */
        EXPAND,
        /**
         * Task to generate tables of various content from a MUSIAL storage file.
         */
        VIEW,
        /**
         * Task to export sequence data in FASTA format from a MUSIAL storage file.
         */
        SEQUENCE,
        /**
         * Task is undefined.
         */
        UNDEFINED
    }

    /**
     * The main entry point of the MUSIAL application.
     * <p>
     * This method initializes the program, determines the task to execute based on the provided arguments,
     * and executes the corresponding functionality. It handles errors gracefully and logs relevant information.
     *
     * @param args Command-line arguments specifying the task and its parameters.
     *             <ul>
     *                 <li><b>args[0]</b>: The task to execute (e.g., BUILD, UPDATE, VIEW, SEQUENCE).</li>
     *                 <li>Additional arguments are parsed by the {@link CLI} class.</li>
     *             </ul>
     */
    public static void main(String[] args) {
        try {
            // Initialize the logging system.
            Logging.init(Level.CONFIG);

            // Record the start time of the program.
            startTime = System.currentTimeMillis();

            // Load metadata such as software name, version, and contact information.
            loadMetadata();

            // Check if any arguments were provided; if not, display usage information and exit.
            if (args.length == 0) {
                System.out.printf("No arguments were specified. Call `java -jar %s-%s.jar [-h|--help]` for more information.%n",
                        Musial.name, Musial.version);
                System.exit(0);
            }

            // Parse the first argument to determine the task to execute.
            try {
                task = Task.valueOf(args[0].toUpperCase());
            } catch (IllegalArgumentException e) {
                // If the task is invalid, set it to UNDEFINED.
                task = Task.UNDEFINED;
            }

            // Parse additional arguments using the CLI utility.
            CLI.parse(args);

            // Execute the task based on the parsed value.
            switch (task) {
                case BUILD -> {
                    Logging.logInfo("Execute task \033[1mbuild\033[0m");
                    UpdateUtility.build();
                }
                case EXPAND -> {
                    Logging.logInfo("Execute task \033[1mexpand\033[0m");
                    UpdateUtility.expand();
                }
                case VIEW -> {
                    Logging.logInfo("Execute task \033[1mtable\033[0m");
                    ContentUtility.run();
                }
                case SEQUENCE -> {
                    Logging.logInfo("Execute task \033[1msequence\033[0m");
                    SequenceUtility.run();
                }
                // Exit the program, if the task is undefined.
                default -> System.exit(-2);
            }
        } catch (Exception e) {
            // Log the error message and stack trace, then exit with an error code.
            if (e.getClass().equals(MusialException.class))
                Logging.logExit("An internal error has occurred.");
            else
                Logging.logExit("An unexpected error has occurred.");
            e.printStackTrace();
            System.exit(-1);
        }
    }

    /**
     * Loads metadata, such as the software title and version from `/src/main/resources/info.properties` and prints the information to stdout.
     *
     * @throws IOException If any metadata can not be loaded.
     */
    private static void loadMetadata() throws IOException {
        Properties properties = new Properties();
        InputStream in = Musial.class.getResourceAsStream("/info.properties");
        properties.load(in);
        Musial.name = properties.getProperty("name");
        Musial.version = properties.getProperty("version");
        Musial.contact = properties.getProperty("contact");
        Musial.license = properties.getProperty("license");
        // Print information to stdout.
        Logging.printSoftwareInfo();
    }

    /**
     * Inner utility class for updating MUSIAL storage files.
     * <p>
     * Handles the {@code BUILD} and {@code EXPAND} tasks.
     */
    public static class UpdateUtility {

        /**
         * Updates the storage by processing variant calls, running annotations, inferring sequence types,
         * and computing statistics.
         *
         * @param storage The {@link Storage} instance to update.
         * @throws IOException     If an I/O error occurs during the update process.
         * @throws MusialException If a MUSIAL-specific error occurs.
         */
        private static void update(Storage storage) throws IOException, MusialException {
            Logging.logInfo("Process variant calls.");
            storage.updateVariants();

            // Check and run SnpEff annotation if applicable.
            if (storage.getSkipSnpEff()) {
                Logging.logInfo("Skip SnpEff analysis per user request.");
            } else if (storage.getHasMissingContigSequences()) {
                Logging.logWarning("Skip SnpEff annotation; no reference sequence for contigs available.");
            } else if (storage.getFeatures().isEmpty()) {
                Logging.logWarning("Skip SnpEff annotation; no features available.");
            } else if (storage.getFeatures().stream().allMatch(f -> f.type.equals("region"))) {
                Logging.logWarning("Skip SnpEff annotation; all features are of type region.");
            } else if (!storage.getHasNovelVariants()) {
                Logging.logWarning("Skip SnpEff annotation; no variants to annotate.");
            } else {
                Logging.logInfo("Run SnpEff annotation.");
                storage.annotateVariants();
            }

            // Infer sequence types if reference sequences are available.
            if (storage.getFeatures().isEmpty()) {
                Logging.logWarning("Skip sequence type inference; no features available.");
            } else {
                Logging.logInfo("Infer sequence types.");
                storage.updateSequenceTypes();
            }

            // Compute statistics for the storage.
            Logging.logInfo("Compute statistics.");
            storage.updateStatistics();
        }

        /**
         * Builds a MUSIAL storage file by initializing storage, processing variant calls, running annotations,
         * inferring sequence types, and computing statistics. The results are written to the specified output file.
         *
         * @throws MusialException If the output file is not specified or other MUSIAL-specific errors occur.
         * @throws IOException     If an I/O error occurs during file operations.
         */
        private static void build() throws MusialException, IOException {
            // Validate the output file parameter.
            String outputPath = (String) CLI.parameters.get("output");
            if (outputPath == null || outputPath.isBlank()) {
                throw new IOException("No valid output file or directory was specified.");
            }

            // Ensure the output path is a file, not a directory.
            File outputFile = new File(outputPath);
            if (outputFile.isDirectory()) {
                outputPath += File.separator + "musial_storage_%s%s".formatted(Logging.getDate(), outputExtension);
                outputFile = new File(outputPath);
            }

            // Set the output directory based on the parent directory of the output file.
            outputDirectory = outputFile.getParentFile();

            // Ensure the output directory exists by creating any necessary parent directories and is writable.
            FileUtils.createParentDirectories(outputDirectory);
            if (!outputDirectory.canWrite()) {
                throw new IOException("No write permission for output directory %s.".formatted(outputDirectory));
            }

            Logging.logInfo("Initialize storage.");
            Storage storage = Storage.Factory.fromCli();

            // Update the storage with variant calls, annotations, and statistics.
            update(storage);

            // Write the storage data to the specified output file.
            Logging.logInfo("Write storage to file: " + CLI.parameters.get("output"));
            Storage.Factory.serialize(storage, outputFile);

            // Log summary information about the storage and execution time.
            long processedGenotypes = storage.getProcessedGenotypesCount();
            float filteredGenotypes = storage.getFilteredGenotypesCount() / (float) processedGenotypes * 100;
            float ignoredGenotypes = storage.getIgnoredGenotypesCount() / (float) processedGenotypes * 100;
            Logging.logDone(
                    "Storage contains %d samples, %d features, %d variants. Processed %d genotypes (%.2f%% filtered, %.2f%% reference or excluded). Execution time: %.2f seconds."
                            .formatted(
                                    storage.getSamples().size(),
                                    storage.getFeatures().size(),
                                    storage.getVariantsCount(),
                                    processedGenotypes,
                                    filteredGenotypes,
                                    ignoredGenotypes,
                                    (System.currentTimeMillis() - startTime) / 1000.0
                            )
            );
        }

        /**
         * Expands an existing MUSIAL storage file by adding new sample data from variant call files,
         * updating annotations, and computing statistics. The results can be written to a new or existing file.
         *
         * @throws IOException     If an I/O error occurs during file operations.
         * @throws MusialException If a MUSIAL-specific error occurs.
         */
        private static void expand() throws IOException, MusialException {
            // Log and start the storage reading process.
            Logging.logInfo("Read storage.");
            File inputFile = new File((String) CLI.parameters.get("input"));
            Storage storage = Storage.Factory.deserialize(inputFile);

            // Record the original sample count and variant count for logging purposes.
            int originalSampleCount = storage.getSamples().size();
            long originalVariantsCount = storage.getVariantsCount();

            // Determine whether the updated storage should be written to a file.
            boolean write = (Boolean) CLI.parameters.get("write");

            // Retrieve the output mode/file path from the CLI parameters.
            String output = (String) CLI.parameters.get("output");
            File outputFile;

            // If the output file is set to "overwrite", use the input file path as the output file.
            if (output.equals("overwrite")) {
                outputFile = inputFile;
            } else {
                // Ensure the output path is a file, not a directory.
                outputFile = new File(output);
                if (outputFile.isDirectory()) {
                    outputFile = new File(outputFile.getAbsolutePath()
                            + File.separator
                            + "musial_storage_%s_%s".formatted(Logging.getDate(), outputExtension)
                    );
                }
            }

            // Set the output directory based on the parent directory of the output file.
            outputDirectory = outputFile.getParentFile();

            // Ensure the output directory exists by creating any necessary parent directories and is writable.
            FileUtils.createParentDirectories(outputDirectory);
            if (!outputDirectory.canWrite()) {
                throw new IOException("No write permission for output directory %s.".formatted(outputDirectory));
            }

            // Add sample information from the specified metadata file, if provided.
            String sampleInfoFile = (String) CLI.parameters.get("vcfMeta");
            if (sampleInfoFile != null) {
                Storage.Factory.cacheSampleInformation(storage, new File(sampleInfoFile));
            }

            // Add VCF files to the storage for processing.
            //noinspection unchecked
            Storage.Factory.updateVcfFiles(storage, (List<String>) CLI.parameters.get("vcfInput"));

            // Update the storage with new data, annotations, and statistics.
            update(storage);

            // Write the updated storage to the specified output file, if the write flag is enabled.
            if (write) {
                Logging.logInfo("Write storage to file: " + outputFile);
                Storage.Factory.serialize(storage, outputFile);
            }

            // Log summary information about the expanded storage and execution time.
            Logging.logDone(
                    "Storage %s with %d samples, %d variants. Processed %d genotypes. Execution time: %.2f seconds."
                            .formatted(
                                    write ? "updated" : "updatable", // Indicate whether the storage was expanded or just expandable.
                                    storage.getSamples().size() - originalSampleCount, // Number of new samples added.
                                    storage.getVariantsCount() - originalVariantsCount, // Number of new variants added.
                                    storage.getProcessedGenotypesCount(), // Total number of genotypes processed.
                                    (System.currentTimeMillis() - startTime) / 1000.0 // Total execution time in seconds.
                            )
            );
        }
    }

    /**
     * Inner utility class for extracting and displaying non-sequence data from a MUSIAL storage file.
     * <p>
     * Handles the {@code VIEW} and {@code PROFILE} tasks.
     */
    public static class ContentUtility {

        /**
         * Table structure for storing and displaying data in a tabular format.
         * <p>
         * This class provides functionality to manage rows and columns, add entries, and generate
         * a string representation of the table. It supports sorting of row identifiers using a
         * custom comparator and allows specifying default content for missing entries.
         */
        private static class Table {

            /**
             * A sorted set of unique identifiers for the rows in the table.
             */
            protected final NavigableSet<String> identifiers;

            /**
             * Comparator used for sorting the row identifiers.
             */
            private final Comparator<String> comparator;

            /**
             * Header for the identifier column.
             */
            protected final String identifierHeader;

            /**
             * List of column headers in the table.
             */
            protected final List<String> headers = new ArrayList<>();

            /**
             * Map storing the content of the table. Each key represents a column header,
             * and the value is a map of row identifiers to cell values.
             */
            protected final Map<String, Map<String, Object>> content = new HashMap<>();

            /**
             * Maximum capacity of the table.
             */
            private final int capacity;

            /**
             * Default content to display for missing entries.
             */
            private final String defaultContent;

            /**
             * Constructs a new `Table` instance.
             *
             * @param idHeader       The header for the identifier column.
             * @param capacity       The maximum capacity of the table.
             * @param comparator     The comparator used for sorting the identifiers.
             * @param defaultContent The default content for missing entries.
             */
            protected Table(String idHeader, int capacity, Comparator<String> comparator, String defaultContent) {
                this.identifiers = new TreeSet<>();
                this.comparator = comparator;
                this.identifierHeader = idHeader;
                this.capacity = capacity;
                this.defaultContent = defaultContent;
            }

            /**
             * Adds a new entry to the table.
             *
             * @param id    The identifier for the row.
             * @param items A list of key-value pairs representing the column header and its value.
             */
            protected void addContent(String id, List<Tuple<String, String>> items) {
                this.identifiers.add(id);
                Set<String> uniqueHeaders = new HashSet<>(headers);
                for (Tuple<String, String> item : items) {
                    if (uniqueHeaders.add(item.a)) {
                        headers.add(item.a);
                    }
                    content.computeIfAbsent(item.a, k -> new HashMap<>(capacity)).put(id, item.b);
                }
            }

            /**
             * Converts the table to a string representation.
             *
             * @return A string representation of the table, including headers and rows.
             */
            @Override
            public String toString() {
                StringBuilder sb = new StringBuilder();

                // Append the header row.
                sb.append(identifierHeader).append("\t").append(String.join("\t", headers)).append("\n");

                // Append each row of the table.
                identifiers.stream().sorted(comparator).forEach(id -> {
                    sb.append(id);
                    for (String header : headers) {
                        sb.append("\t").append(content.getOrDefault(header, Map.of()).getOrDefault(id, defaultContent));
                    }
                    sb.append("\n");
                });
                return sb.toString();
            }
        }

        /**
         * Set of supported content types for the table and profile tasks.
         * <p>
         * Defines the types of data that can be generated and displayed in
         * tabular format by the table and profile tasks of MUSIAL.
         */
        public static final Set<String> content = Set.of(
                "features", // Induces .run() to generate a table of features.
                "samples", // ... samples.
                "variants", // ... variants.
                "alleles", // ... alleles of one feature.
                "proteoforms", // ... proteoforms of one feature.
                "calls", // Induces .run() to generate a sample profile of variant calls.
                "types" // ... sequence types.
        );

        /**
         * Generates and displays a table based on the specified content type and filters.
         * <p>
         * This method reads the storage file, applies filters for features, samples, and positions,
         * and generates a table for the specified content type. The table can be displayed on the console
         * or written to a file.
         *
         * @throws IOException     If an I/O error occurs during file operations.
         * @throws MusialException If an error specific to MUSIAL occurs.
         */
        private static void run() throws IOException, MusialException {
            // Log and start the storage reading process.
            Logging.logInfo("Read storage.");
            File inputFile = new File((String) CLI.parameters.get("input"));
            Storage storage = Storage.Factory.deserialize(inputFile);

            // Retrieve and validate the content type to view.
            String content = ((String) CLI.parameters.get("content")).toLowerCase();
            if (!content.matches(String.join(Constants.pipe, ContentUtility.content))) {
                throw new MusialException("Content (-c) has to be one of %s, but %s was provided."
                        .formatted(String.join(", ", ContentUtility.content), content));
            }

            // Initialize sets to store filters for features, samples, and positions.
            Set<String> features = new HashSet<>(), samples = new HashSet<>(), positions = new HashSet<>();
            // Parse filter parameter to populate filter sets.
            //noinspection unchecked
            for (String value : (Set<String>) CLI.parameters.get("filter")) {
                if (storage.hasFeature(value)) features.add(value);
                if (storage.hasSample(value)) samples.add(value);
                if (value.matches("\\d+")) positions.add(value);
            }

            // sample.(gyrA.allele=[]10)

            // Retrieve and validate the output destination.
            String output = (String) CLI.parameters.get("output");
            File outputFile = null;
            if (!output.equals("stdout")) {
                // Ensure the output path is a file, not a directory.
                outputFile = new File(output);
                if (outputFile.isDirectory()) {
                    outputFile = new File(outputFile.getAbsolutePath()
                            + File.separator
                            + "musial_view_%s_%s.tsv".formatted(content, Logging.getDate())
                    );
                }

                // Set the output directory based on the parent directory of the output file.
                outputDirectory = outputFile.getParentFile();

                // Ensure the output directory exists by creating any necessary parent directories and is writable.
                FileUtils.createParentDirectories(outputDirectory);
                if (!outputDirectory.canWrite()) {
                    throw new IOException("No write permission for output directory %s.".formatted(outputDirectory));
                }
            }

            // Generate the table based on the specified content type.
            Logging.logInfo("Generate `%s` content.".formatted(content));
            Table table = switch (content) {
                case "features" -> featureTable(storage, features);
                case "samples" -> sampleTable(storage, samples);
                case "variants" -> variantTable(storage, positions, samples, features);
                case "calls" -> callMatrix(storage, samples, positions);
                default -> throw new MusialException("Unknown task `view` content %s.".formatted(content));
            };

            // Handle the case where no entries match the filters.
            if (table.identifiers.isEmpty()) {
                Logging.logWarning("No entries to view. Check your filter parameter.");
                // Log the completion of the task with the execution time.
                Logging.logDone("Execution time %.2f seconds.".formatted((System.currentTimeMillis() - startTime) / 1000.0));
            }
            // Output the table to the console if "stdout" is specified.
            else if ("stdout".equals(output)) {
                // Log the completion of the task with the execution time.
                Logging.logDone("Execution time %.2f seconds.".formatted((System.currentTimeMillis() - startTime) / 1000.0));
                System.out.println(table);
            }
            // Write the table to the specified file.
            else {
                IO.writeFile(outputFile.toPath(), table.toString());
                Logging.logInfo("Write results to %s.".formatted(output));
                // Log the completion of the task with the execution time.
                Logging.logDone("Execution time %.2f seconds.".formatted((System.currentTimeMillis() - startTime) / 1000.0));
            }
        }

        /**
         * Generates a table containing information about genomic features.
         * <p>
         * This method creates a table with rows representing features and columns representing
         * various attributes of each feature. The table can be filtered to include only specific
         * features based on the provided set of feature names.
         *
         * @param storage The {@link Storage} instance containing the features to be included in the table.
         * @param include A set of feature names to include in the table. If empty, all features are included.
         * @return A {@link Table} object containing the feature information.
         */
        private static Table featureTable(Storage storage, Set<String> include) {
            // Initialize the table with the header "name" and a comparator for sorting by feature start position.
            Table table = new Table("name", include.isEmpty() ? storage.getFeatures().size() : include.size(),
                    Comparator.comparingInt(i -> storage.getFeature(i).start), Constants.empty);

            // Stream through the features in the storage, filtering based on the include set.
            storage.getFeatures().stream()
                    .filter(feature -> include.isEmpty() || include.contains(feature.name))
                    .forEach(feature -> {
                        // Create a list of tuples representing the feature's attributes.
                        List<Tuple<String, String>> items = new ArrayList<>(List.of(
                                new Tuple<>("chromosome", feature.contig),
                                new Tuple<>("start", String.valueOf(feature.start)),
                                new Tuple<>("end", String.valueOf(feature.end)),
                                new Tuple<>("strand", String.valueOf(feature.strand)),
                                new Tuple<>("type", feature.type)
                        ));

                        // Add additional attributes of the feature to the list.
                        feature.getAttributes().forEach((key, value) -> items.add(new Tuple<>(key, value)));

                        // Add the feature's name and its attributes to the table.
                        table.addContent(feature.name, items);
                    });

            // Return the populated table.
            return table;
        }

        /**
         * Generates a table containing information about samples.
         * <p>
         * This method creates a table with rows representing samples and columns representing
         * various attributes of each sample. The table can be filtered to include only specific
         * samples based on the provided set of sample names.
         *
         * @param storage The {@link Storage} instance containing the samples to be included in the table.
         * @param include A set of sample names to include in the table. If empty, all samples are included.
         * @return A {@link Table} object containing the sample information.
         */
        private static Table sampleTable(Storage storage, Set<String> include) {
            // Initialize the table with the header "name" and a comparator for natural ordering of sample names.
            Table table = new Table("name", include.isEmpty() ? storage.getSamples().size() : include.size(),
                    Comparator.naturalOrder(), Constants.empty);

            // Stream through the samples in the storage, filtering based on the include set.
            storage.getSamples().stream()
                    .filter(sample -> include.isEmpty() || include.contains(sample.name))
                    .forEach(sample -> {
                        // Create a list of tuples representing the sample's attributes.
                        List<Tuple<String, String>> items = sample.getAttributes().entrySet().stream()
                                .map(entry -> new Tuple<>(entry.getKey(), entry.getValue()))
                                .toList();

                        // Add the sample's name and its attributes to the table.
                        table.addContent(sample.name, items);
                    });

            // Return the populated table.
            return table;
        }

        /**
         * Generates a table containing information about genomic variants.
         * <p>
         * This method creates a table with rows representing variants and columns representing
         * various attributes of each variant. The table can be filtered to include only specific
         * positions, samples, and features based on the provided sets.
         *
         * @param storage          The {@link Storage} instance containing the variants to be included in the table.
         * @param includePositions A set of positions to include in the table. If empty, all positions are included.
         * @param includeSamples   A set of sample names to include in the table. If empty, all samples are included.
         * @param includeFeatures  A set of feature names to include in the table. If empty, all features are included.
         * @return A {@link Table} object containing the variant information.
         */
        private static Table variantTable(Storage storage, Set<String> includePositions, Set<String> includeSamples, Set<String> includeFeatures) {
            // Initialize the table with the header and a comparator for sorting by position.
            Table table = new Table("contig\tpos\tref\talt", (int) storage.getVariantsCount(),
                    Comparator.comparingInt(i -> Integer.parseInt(i.split(Constants.tab)[1])), Constants.empty);

            // Iterate through each contig in the storage.
            storage.getContigs().forEach(contig ->
                    // Iterate through each variant in the contig.
                    contig.getVariants().forEach(variant -> {
                        // Check if the variant's position is included in the filter set.
                        if (includePositions.isEmpty() || includePositions.contains(String.valueOf(variant.a))) {
                            // Retrieve variant information for the current variant.
                            VariantInformation variantInfo = contig.getVariantInformation(variant.a, variant.b);

                            // Check if the variant is associated with any of the included features.
                            boolean hasFeature = includeFeatures.isEmpty() || includeFeatures.stream().anyMatch(variantInfo::hasOccurrence);

                            // Check if the variant is associated with any of the included samples.
                            boolean hasSample = includeSamples.isEmpty() || includeSamples.stream().anyMatch(
                                    sample -> variantInfo.hasOccurrence(Attributes.sampleOccurrence, sample));

                            // If the variant matches the feature and sample filters, add it to the table.
                            if (hasFeature && hasSample) {
                                // Create a list of tuples representing the variant's attributes.
                                List<Tuple<String, String>> items = new ArrayList<>(List.of(
                                        new Tuple<>("type", variantInfo.type.name())
                                ));

                                // Add additional attributes of the variant to the list.
                                variantInfo.getAttributes().forEach((key, value) -> items.add(new Tuple<>(key, value)));

                                // Add the occurrence information of the variant to the list.
                                items.add(new Tuple<>("samples", String.join(Constants.comma, variantInfo.getSampleOccurrence())));

                                // Add the variant's information to the table.
                                table.addContent(contig.name + "\t" + variant.a + "\t" + variantInfo.reference + "\t" + variant.b, items);
                            }
                        }
                    })
            );

            // Return the populated table.
            return table;
        }

        /**
         * Generates a table containing information about variant calls for samples.
         * <p>
         * This method creates a table with rows representing variant calls and columns representing
         * the contig, position, reference, and sample-specific call information. The table can be
         * filtered to include only specific samples and positions based on the provided sets.
         *
         * @param storage           The {@link Storage} instance containing the variant calls to be included in the table.
         * @param includedSamples   A set of sample names to include in the table. If empty, all samples are included.
         * @param includedPositions A set of positions to include in the table. If empty, all positions are included.
         * @return A {@link Table} object containing the variant call information.
         */
        private static Table callMatrix(Storage storage, Set<String> includedSamples, Set<String> includedPositions) {
            // Initialize the table with the header and a comparator for sorting by position.
            Table table = new Table("contig\tposition\treference", (int) storage.getVariantsCount(),
                    Comparator.comparingInt(s -> Integer.parseInt(s.split("\t")[1])), Constants.dot);

            // Stream through the samples in the storage, filtering based on the includedSamples set.
            storage.getSamples().stream()
                    .filter(sample -> includedSamples.isEmpty() || includedSamples.contains(sample.name))
                    .forEach(sample ->
                            // Stream through the contigs for each sample.
                            storage.getContigs().forEach(contig ->
                                    // Stream through the variant calls for each contig, filtering based on the includedPositions set.
                                    sample.getVariantCalls(contig.name).entrySet().stream()
                                            .filter(variantCall -> includedPositions.isEmpty() || includedPositions.contains(String.valueOf(variantCall.getKey())))
                                            .forEach(variantCall -> {
                                                // Create a list of tuples representing the variant call's attributes.
                                                List<Tuple<String, String>> items = List.of(new Tuple<>(sample.name, variantCall.getValue()));

                                                // Add the variant call's information to the table.
                                                table.addContent(
                                                        contig.name + "\t" + variantCall.getKey() + "\t" + Sample.getReferenceOfCall(variantCall.getValue()),
                                                        items
                                                );
                                            })
                            )
                    );

            // Return the populated table.
            return table;
        }
    }

    /**
     * Inner utility class for extracting sequence data from a MUSIAL storage file.
     * <p>
     * Handles the {@code SEQUENCE} task.
     */
    public static class SequenceUtility {

        /**
         * Executes the sequence export task for the MUSIAL application.
         * <p>
         * This method validates the output directory, reads the storage file, and exports
         * nucleotide or amino acid sequences for the specified features and samples based
         * on the provided parameters. The sequences are written to the specified output directory.
         *
         * @throws IOException     If an I/O error occurs during file operations.
         * @throws MusialException If a MUSIAL-specific error occurs, such as missing parameters or invalid paths.
         */
        private static void run() throws IOException, MusialException {
            // Log and start the storage reading process.
            Logging.logInfo("Read storage.");
            File inputFile = new File((String) CLI.parameters.get("input"));
            Storage storage = Storage.Factory.deserialize(inputFile);

            // Validate the output file parameter.
            String output = (String) CLI.parameters.get("output");
            if (output.equals("parent")) {
                outputDirectory = inputFile.getParentFile();
            } else {
                File outputFile = new File(output);
                if (!outputFile.isDirectory()) {
                    throw new IOException("Output path is not a directory.");
                }
                outputDirectory = outputFile;
                FileUtils.createParentDirectories(outputDirectory);
            }

            // Retrieve task parameters from the CLI.
            String content = (String) CLI.parameters.get("content");
            if (!content.matches("nt|aa")) {
                throw new MusialException("Content (-c) has to be one of nt or aa, but %s was provided.".formatted(content));
            }
            boolean nt = content.equals("nt"); // Determines if nucleotide sequences are exported.
            boolean merge = (Boolean) CLI.parameters.get("merge"); // If true, merges sequences for all samples.
            boolean strip = (Boolean) CLI.parameters.get("strip"); // If true, un-aligns sequences by removing gaps.
            boolean conserved = (Boolean) CLI.parameters.get("conserved"); // If true, includes conserved reference content.
            boolean reference = (Boolean) CLI.parameters.get("reference"); // If true, includes the reference sequence.

            // Retrieve the list of features and samples to process.
            //noinspection unchecked
            HashSet<String> featureNames = (HashSet<String>) CLI.parameters.get("features");
            //noinspection unchecked
            HashSet<String> sampleNames = (HashSet<String>) CLI.parameters.get("samples");

            // If no samples are specified, include all samples from the storage.
            if (sampleNames.isEmpty()) {
                sampleNames.addAll(storage.getSamples().stream().map(s -> s.name).collect(Collectors.toSet()));
            }

            // Log the start of the sequence export process.
            Logging.logInfo("Export sequences.");

            // Iterate through each feature and export its sequences.
            for (String featureName : featureNames) {

                // Retrieve the corresponding feature from storage.
                Feature feature = storage.getFeature(featureName.split(":")[0]);
                if (Objects.isNull(feature)) {
                    throw new MusialException("Feature %s not available in storage.".formatted(featureName));
                }

                Contig contig = storage.getContig(feature.contig);
                int from, to;

                // If the feature name contains a region specification (e.g., "feature:100..200"), parse and validate the region.
                if (featureName.contains(":") && featureName.matches("^.+:g.[0-9]+\\.\\.[0-9]+$")) {
                    String[] region = featureName.split(":g.")[1].split("\\.\\.");
                    from = Integer.parseInt(region[0]);
                    to = Integer.parseInt(region[1]);

                    if (from < feature.start || to > feature.end) {
                        throw new MusialException("Specified region %d..%d is out of bounds for feature %s (%d..%d)."
                                .formatted(from, to, feature.name, feature.start, feature.end));
                    }
                } else { // If no region is specified, use the full range of the feature.
                    from = feature.start;
                    to = feature.end;
                }

                // Export nucleotide or amino acid sequences based on the task parameters.
                if (nt) {
                    exportNtSequences(contig, feature, from, to, sampleNames, conserved, merge, strip, reference);
                } else {
                    exportAaSequences(contig, feature, from, to, sampleNames, conserved, merge, strip, reference);
                }
            }

            // Log the completion of the task with the execution time.
            Logging.logDone("Execution time %.2f seconds.".formatted(((float) (System.currentTimeMillis() - startTime) / 1000)));
        }

        /**
         * Exports nucleotide sequences for a given feature and contig to a FASTA file.
         * <p>
         * This method processes variants, resolves reference sequences, and generates
         * nucleotide sequences for alleles based on the provided parameters. The sequences
         * are written to a FASTA file in the specified output directory.
         *
         * @param contig         The contig containing the feature and its variants.
         * @param feature        The genomic feature for which sequences are exported.
         * @param from           The start position of the sequence to export.
         * @param to             The end position of the sequence to export.
         * @param sampleNames    A set of sample names to filter alleles for sequence generation.
         * @param writeConserved If true, generates sequences with conserved reference content.
         * @param merge          If true, merges sequences for all samples into a single output.
         * @param strip          If true, un-aligns sequences by removing gaps.
         * @param writeReference If true, includes the reference sequence in the output.
         * @throws IOException If an I/O error occurs during file writing.
         */
        private static void exportNtSequences(Contig contig, Feature feature, int from, int to, Set<String> sampleNames,
                                              boolean writeConserved, boolean merge, boolean strip, boolean writeReference) throws IOException {

            // Check if conserved sequences are requested but the contig lacks reference sequence information.
            if (writeConserved && !contig.hasSequence()) {
                Logging.logWarning("Skip feature %s as contig %s has no reference sequence information (incompatible with conserved export)."
                        .formatted(feature.name, contig.name));
                return;
            }

            // Collect allele UIDs that match with the provided sample names.
            final Set<String> alleleUids = feature.getAlleles().stream()
                    .filter(allele -> sampleNames.stream().anyMatch(allele::hasOccurrence))
                    .sorted(Comparator.comparing(SequenceType::getCount))
                    .map(allele -> allele.uid)
                    .collect(Collectors.toSet());
            if (alleleUids.isEmpty()) {
                Logging.logWarning("Skipping feature %s: Only sequence type is the reference sequence."
                        .formatted(feature.name));
                return;
            }

            // Retrieve variants associated with the selected alleles.
            ArrayList<Tuple<Integer, String>> variants = contig.getVariantsByAlleles(feature, alleleUids);

            // Filter variants to only include those within the specified position range.
            if (feature.start != from || feature.end != to) {
                variants = variants.stream()
                        .filter(variant -> variant.a >= from && variant.a <= to)
                        .collect(Collectors.toCollection(ArrayList::new));
            }

            // Retrieve the reference content if conserved sequences are requested.
            final char[] referenceContent;
            if (writeConserved) referenceContent = contig.getSubsequence(from, to).toCharArray();
            else referenceContent = null;

            // Store context (reference content and max. indel length) per position.
            HashMap<Integer, Tuple<String, Integer>> positionalContext = new HashMap<>();

            // Function to update the positional context from variant information.
            BiConsumer<Tuple<Integer, String>, Integer> updatePositionalContext = (context, insertionLength) ->
                    positionalContext.merge(context.a, new Tuple<>(context.b, insertionLength), (e1, e2) -> {
                        if (!Objects.equals(e2.a, Constants.empty) && !Objects.equals(e1.a, Constants.empty) && !Objects.equals(e1.a, e2.a)) {
                            Logging.logWarning("Reference content conflict at position %d (%s and %s).".formatted(context.a, e1.a, e2.a));
                        }
                        return new Tuple<>(e1.a.equals(Constants.empty) ? e2.a : e1.a, Math.max(e1.b, e2.b));
                    });

            // Process each variant to populate positional context.
            VariantInformation variantInformation;
            for (Tuple<Integer, String> variant : variants) {
                variantInformation = contig.getVariantInformation(variant.a, variant.b);
                char[] ref = variantInformation.reference.toCharArray();
                if (variantInformation.type.equals(VariantInformation.Type.SNV)) {
                    updatePositionalContext.accept(new Tuple<>(variant.a, variantInformation.reference), 0);
                } else if (variantInformation.type.equals(VariantInformation.Type.DELETION)) {
                    if (variant.b.charAt(0) != ref[0])
                        updatePositionalContext.accept(new Tuple<>(variant.a, String.valueOf(ref[0])), 0);
                    for (int i = 1; i < ref.length; i++) {
                        updatePositionalContext.accept(new Tuple<>(variant.a + i, String.valueOf(ref[i])), 0);
                    }
                } else if (variantInformation.type.equals(VariantInformation.Type.INSERTION)) {
                    int length = variant.b.length() - 1;
                    if (variant.b.charAt(0) != ref[0])
                        updatePositionalContext.accept(new Tuple<>(variant.a, String.valueOf(ref[0])), length);
                    else updatePositionalContext.accept(new Tuple<>(variant.a, Constants.empty), length);
                }
            }

            // String builder to store sequence content.
            StringBuilder content = new StringBuilder(writeConserved ? referenceContent.length : variants.size());

            // Function to resolve reference content for a position.
            Consumer<Integer> resolveReference = position -> {
                if (positionalContext.containsKey(position)) {
                    Tuple<String, Integer> context = positionalContext.get(position);
                    String referenceBase = context.a.isEmpty() && writeConserved ? String.valueOf(referenceContent[position - from]) : context.a;
                    content.append(SequenceOperations.padGaps(referenceBase, referenceBase.length() + context.b));
                } else if (writeConserved) {
                    content.append(referenceContent[position - from]);
                }
            };

            // Write sequences to FASTA file.
            String fileName = String.format("%s_%d_%d_%s_%s%s_nt.fasta",
                    feature.name,
                    from,
                    to,
                    writeConserved ? "conserved" : "variable",
                    merge ? "merged" : "samples",
                    strip ? "" : "_aligned"
            );
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputDirectory + File.separator + fileName, StandardCharsets.UTF_8))) {
                // Function to write header and sequence content to the file.
                BiConsumer<String, String> dump = (header, sequence) -> {
                    try {
                        writer.write("%s\n%s\n".formatted(header, String.join("\n", Splitter.fixedLength(80).split(sequence))));
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                };

                // Write the reference sequence if requested.
                if (writeReference) {
                    // Generate the reference sequence.
                    IntStream.rangeClosed(from, to).forEach(resolveReference::accept);

                    // Build the header.
                    StringBuilder header = new StringBuilder(">%s.%s".formatted(feature.name, Constants.reference));
                    long count = sampleNames.stream()
                            .filter(sampleName -> feature.getAlleles().stream().noneMatch(allele -> allele.hasOccurrence(sampleName)))
                            .count();
                    header.append(" N=%d/%d".formatted(count, sampleNames.size()));

                    String attributes = feature.attributesAsString(
                            Set.of(Constants.$Feature_children, Constants.$Attributable_frequencyDisrupted,
                                    Constants.$Attributable_frequencyReference, Constants.$Feature_numberOfAlleles,
                                    Constants.$Feature_numberOfProteoforms),
                            Constants.pipe
                    );
                    if (!attributes.isEmpty()) {
                        header.append(" ").append(attributes);
                    }

                    // Remove gaps if requested.
                    String sequence = strip ? content.toString().replace(Constants.gap, Constants.empty) : content.toString();

                    // Write the sequence to the file.
                    dump.accept(header.toString(), sequence);
                }

                // Write sequences for each allele.
                Feature.Allele allele;
                String alt;
                boolean positionConserved;
                int deletedPositions;
                Tuple<String, Integer> context;
                Tuple<String, Integer> nullContext = new Tuple<>(Constants.empty, 0);
                for (String alleleUid : alleleUids) {
                    allele = feature.getAllele(alleleUid);
                    content.setLength(0);
                    deletedPositions = 0;
                    for (int position = from; position <= to; position++) {
                        if (deletedPositions > 0) {
                            context = positionalContext.get(position);
                            content.append(SequenceOperations.padGaps(Constants.gap, 1 + context.b));
                            deletedPositions--;
                            if (allele.hasVariant(position))
                                Logging.logWarning("Conflict with variant %s at deleted position %d for allele %s of feature %s."
                                        .formatted(allele.getVariant(position), position, alleleUid, feature.name));
                        } else if (allele.hasVariant(position)) {
                            alt = allele.getVariant(position);
                            context = positionalContext.getOrDefault(position, nullContext);
                            positionConserved = context.a.isEmpty();
                            if (VariantInformation.isSubstitution(alt)) {
                                content.append(SequenceOperations.padGaps(alt, 1 + context.b));
                            } else if (VariantInformation.isDeletion(alt)) {
                                if (positionConserved && !writeConserved) {
                                    content.append(SequenceOperations.padGaps(Constants.empty, context.b));
                                } else {
                                    content.append(SequenceOperations.padGaps(alt.substring(0, 1), 1 + context.b));
                                }
                                deletedPositions += (alt.length() - 1);
                            } else if (VariantInformation.isInsertion(alt)) {
                                if (positionConserved && !writeConserved) {
                                    content.append(SequenceOperations.padGaps(alt.substring(1), context.b));
                                } else {
                                    content.append(SequenceOperations.padGaps(alt, 1 + context.b));
                                }
                            }
                        } else {
                            resolveReference.accept(position);
                        }
                    }

                    // Remove gaps if requested.
                    String sequence = strip ? content.toString().replace(Constants.gap, Constants.empty) : content.toString();

                    if (merge) {
                        // Build the header.
                        StringBuilder header = new StringBuilder(">%s".formatted(alleleUid));
                        long count = sampleNames.stream().filter(allele::hasOccurrence).count();
                        header.append(" N=%d/%d".formatted(count, sampleNames.size()));

                        String attributes = allele.attributesAsString(Constants.pipe);
                        if (!attributes.isEmpty()) {
                            header.append(" ").append(attributes);
                        }

                        // Write the sequence to the file.
                        dump.accept(header.toString(), sequence);
                    } else {
                        for (String sampleName : allele.getOccurrence()) {
                            if (sampleNames.contains(sampleName)) {
                                // Build the header.
                                StringBuilder header = new StringBuilder(">%s.%s".formatted(feature.name, sampleName));

                                String attributes = allele.attributesAsString(Constants.pipe);
                                if (!attributes.isEmpty()) {
                                    header.append(" ").append(attributes);
                                }

                                // Write the sequence to the file.
                                dump.accept(header.toString(), sequence);
                            }
                        }
                    }
                }
            }
        }

        /**
         * Exports amino acid sequences for a given feature and contig to a FASTA file.
         * <p>
         * This method processes variants, resolves reference sequences, and generates
         * amino acid sequences for proteoforms based on the provided parameters. The sequences
         * are written to a FASTA file in the specified output directory.
         *
         * @param contig      The contig containing the feature and its variants.
         * @param feature     The genomic feature for which sequences are exported.
         * @param from        The start position of the sequence to export.
         * @param to          The end position of the sequence to export.
         * @param sampleNames A set of sample names to filter proteoforms for sequence generation.
         * @param conserved   If true, generates sequences with conserved reference content.
         * @param merge       If true, merges sequences for all samples into a single output.
         * @param strip       If true, un-aligns sequences by removing gaps.
         * @param reference   If true, includes the reference sequence in the output.
         * @throws IOException     If an I/O error occurs during file writing.
         * @throws MusialException If a MUSIAL-specific error occurs.
         */
        private static void exportAaSequences(Contig contig, Feature feature, int from, int to, Set<String> sampleNames,
                                              boolean conserved, boolean merge, boolean strip, boolean reference) throws IOException, MusialException {
            // Check if the contig has reference sequence information; required for amino acid export.
            if (!contig.hasSequence()) {
                Logging.logWarning("Skip feature %s; contig %s has to have reference sequence information for amino acid export."
                        .formatted(feature.name, contig.name));
                return;
            }

            // Validate that the feature is coding.
            if (!feature.isCoding()) {
                throw new MusialException("Feature %s is not coding, cannot export amino acid sequences."
                        .formatted(feature.name));
            }

            // Collect proteoform UIDs that match the provided sample names, excluding synonymous proteoforms.
            final List<String> proteoformUids = feature.getAlleles().stream()
                    .filter(allele -> sampleNames.stream().anyMatch(allele::hasOccurrence))
                    .map(allele -> allele.getAttribute(Constants.$Allele_proteoform))
                    .collect(Collectors.toList());
            proteoformUids.remove(Constants.synonymous);
            if (proteoformUids.isEmpty()) {
                Logging.logWarning("Skip feature %s; all proteoforms are synonymous.".formatted(feature.name));
                return;
            }

            // Collect all variants associated with the selected proteoforms.
            Set<Tuple<Integer, String>> variantsSet = new HashSet<>();
            for (String proteoformUid : proteoformUids) {
                Feature.Proteoform proteoform = feature.getProteoform(proteoformUid);
                proteoform.getVariants().forEach((key, value) -> variantsSet.add(new Tuple<>(key, value)));
            }
            ArrayList<Tuple<Integer, String>> variants = new ArrayList<>(variantsSet);
            variants.sort(Comparator.comparingInt(i -> i.a));

            int fromRelative = ((from - feature.start + 1) + 2) / 3;
            int toRelative = (to - feature.start + 1) / 3;

            if (feature.start != from || feature.end != to) {
                // Check if the range is a multiple of 3, as amino acid sequences require this.
                if ((to - from + 1) % 3 != 0) {
                    throw new MusialException("Amino acid sequence export requires range to be a multiple of 3, but %d..%d is not."
                            .formatted(from, to));
                }

                // Filter variants to only include those within the specified range.
                variants = variants.stream()
                        .filter(variant -> variant.a >= fromRelative && variant.a <= toRelative)
                        .collect(Collectors.toCollection(ArrayList::new));
            }

            // Translate the reference nucleotide sequence to amino acid sequence.
            final char[] referenceContent = SequenceOperations
                    .translateSequence(contig.getSubsequence(from, to), feature.isReverse()).toCharArray();

            // Map to store positional context for variants.
            HashMap<Integer, Tuple<String, Integer>> positionalContext = new HashMap<>();

            // Function to update the positional context with variant information.
            BiConsumer<Tuple<Integer, String>, Integer> updatePositionalContext = (context, insertionLength) ->
                    positionalContext.merge(context.a, new Tuple<>(context.b, insertionLength), (e1, e2) -> {
                        if (!Objects.equals(e2.a, Constants.empty) && !Objects.equals(e1.a, Constants.empty) && !Objects.equals(e1.a, e2.a)) {
                            Logging.logWarning("Reference content conflict at position %d (%s and %s).".formatted(context.a, e1.a, e2.a));
                        }
                        return new Tuple<>(e1.a.equals(Constants.empty) ? e2.a : e1.a, Math.max(e1.b, e2.b));
                    });

            // Process each variant to populate the positional context.
            for (Tuple<Integer, String> variant : variants) {
                int referenceContentIndex = variant.a - fromRelative;
                if (VariantInformation.isSubstitution(variant.b)) {
                    updatePositionalContext.accept(new Tuple<>(variant.a, String.valueOf(referenceContent[referenceContentIndex])), 0);
                } else {
                    boolean mixed = variant.b.charAt(0) != referenceContent[referenceContentIndex];
                    if (VariantInformation.isDeletion(variant.b)) {
                        if (mixed) {
                            updatePositionalContext.accept(new Tuple<>(variant.a, String.valueOf(referenceContent[referenceContentIndex])), 0);
                        }
                        for (int i = 1; i < variant.b.length(); i++) {
                            updatePositionalContext.accept(new Tuple<>(variant.a + i, String.valueOf(referenceContent[referenceContentIndex + i])), 0);
                        }
                    } else if (VariantInformation.isInsertion(variant.b)) {
                        int length = variant.b.length() - 1;
                        if (mixed) {
                            updatePositionalContext.accept(new Tuple<>(variant.a, String.valueOf(referenceContent[referenceContentIndex])), length);
                        } else {
                            updatePositionalContext.accept(new Tuple<>(variant.a, Constants.empty), length);
                        }
                    }
                }
            }

            // StringBuilder to construct the sequence content.
            StringBuilder content = new StringBuilder(conserved ? referenceContent.length : variants.size());

            // Function to resolve reference content for a given position.
            Consumer<Integer> resolveReference = position -> {
                int referenceContentIndex = position - fromRelative;
                if (positionalContext.containsKey(position)) {
                    Tuple<String, Integer> context = positionalContext.get(position);
                    String referenceBase = context.a.isEmpty() && conserved ? String.valueOf(referenceContent[referenceContentIndex]) : context.a;
                    content.append(SequenceOperations.padGaps(referenceBase, referenceBase.length() + context.b));
                } else if (conserved) {
                    content.append(referenceContent[referenceContentIndex]);
                }
            };

            // Write the sequences to a FASTA file.
            String fileName = String.format("%s%s%s%s_%s.faa",
                    feature.name,
                    conserved ? "_conserved" : "_variant",
                    merge ? "_merged" : "_sample",
                    strip ? "" : "_aligned",
                    Logging.getDate()
            );
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputDirectory + File.separator + fileName, StandardCharsets.UTF_8))) {

                // Function to write a sequence to the file with a given header.
                Consumer<String> dump = header -> {
                    try {
                        String sequence = strip ? content.toString().replaceAll(Constants.gap, Constants.empty) : content.toString();
                        writer.write("%s\n%s\n".formatted(header, String.join("\n", Splitter.fixedLength(80).split(sequence))));
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                };

                // Write the reference sequence if requested.
                if (reference) {
                    IntStream.rangeClosed(fromRelative, toRelative).forEach(resolveReference::accept);
                    dump.accept(">reference [allelic_frequency=%s]".formatted(feature.getAttribute(Constants.$Attributable_frequencyReference)));
                }

                // Write sequences for each proteoform.
                for (String proteoformUid : proteoformUids) {
                    Feature.Proteoform proteoform = feature.getProteoform(proteoformUid);
                    content.setLength(0);
                    int deletedPositions = 0;
                    for (int position = fromRelative; position <= toRelative; position++) {
                        if (deletedPositions > 0) {
                            Tuple<String, Integer> context = positionalContext.get(position);
                            content.append(SequenceOperations.padGaps(Constants.gap, 1 + context.b));
                            deletedPositions--;
                            if (proteoform.hasVariant(position))
                                Logging.logWarning("Conflict with variant %s at deleted position %d for allele %s of feature %s."
                                        .formatted(proteoform.getVariant(position), position, proteoformUid, feature.name));
                        } else if (proteoform.hasVariant(position)) {
                            String alt = proteoform.getVariant(position);
                            Tuple<String, Integer> context = positionalContext.getOrDefault(position, new Tuple<>(Constants.empty, 0));
                            boolean positionConserved = context.a.isEmpty();
                            if (VariantInformation.isSubstitution(alt)) {
                                content.append(SequenceOperations.padGaps(alt, 1 + context.b));
                            } else if (VariantInformation.isDeletion(alt)) {
                                if (positionConserved && !conserved) {
                                    content.append(SequenceOperations.padGaps(Constants.empty, context.b));
                                } else {
                                    content.append(SequenceOperations.padGaps(alt.substring(0, 1), 1 + context.b));
                                }
                                deletedPositions += (alt.length() - 1);
                            } else if (VariantInformation.isInsertion(alt)) {
                                if (positionConserved && !conserved) {
                                    content.append(SequenceOperations.padGaps(alt.substring(1), context.b));
                                } else {
                                    content.append(SequenceOperations.padGaps(alt, 1 + context.b));
                                }
                            }
                        } else {
                            resolveReference.accept(position);
                        }
                    }
                    if (merge) {
                        dump.accept(proteoform.getFastaHeader(proteoform.getIdentifier()));
                    } else {
                        for (String alleleUid : proteoform.getOccurrence()) {
                            for (String sampleName : feature.getAllele(alleleUid).getOccurrence()) {
                                if (sampleNames.isEmpty() || sampleNames.contains(sampleName))
                                    dump.accept(">%s [proteoform=%s]".formatted(sampleName, proteoform.getIdentifier()));
                            }
                        }
                    }
                }
            }
        }

    }

}