package main;

import com.google.common.base.Splitter;
import datastructure.*;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import org.apache.commons.io.FileUtils;
import utility.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Main class of MUSIAL (MUlti Sample varIant AnaLysis), a tool to calculate SNV, gene, and whole genome alignments,
 * together with other relevant statistics based on .vcf files.
 */
public final class Musial {

    /**
     * Original system output stream.
     */
    public static final PrintStream originalSysOut = System.out;

    /**
     * Alternative empty output stream to ignore logging.
     */
    public static final PrintStream emptySysOut = new PrintStream(new OutputStream() {
        public void write(int b) {
        }
    });

    /**
     * Unique run ID of this program instance; generated from the current date and time.
     */
    public static final String runId = IO.md5Hash(Logging.getTimestamp());

    /**
     * Name of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String softwareName = "";

    /**
     * Version of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String softwareVersion = "";

    /**
     * Author contact of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String softwareContact = "";

    /**
     * License information of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String softwareLicense = "";

    /**
     * Specifies the task to execute.
     */
    public static Task task;

    /**
     * Specifies the output directory for the program.
     */
    public static File outputDirectory;

    /**
     * Start time of the program.
     */
    public static long startTime;

    /**
     * File extension used for MUSIAL storage files.
     * <p>
     * This constant specifies the default file extension for storage files
     * created or used by the MUSIAL application. It is marked as `transient`
     * to indicate that it should not be serialized as part of the class state.
     * <p>
     * The extension should either be `.json` or `.json.gz` depending on the
     * compression method used. For production, use `.json.gz` for compressed files.
     */
    public static final String storageExtension = ".json.gz";

    /**
     * {@link Enum} specifying tasks of MUSIAL to choose from.
     */
    public enum Task {
        /**
         * Task to build a MUSIAL storage file.
         */
        BUILD,
        /**
         * Task to add sample data from variant call files to a MUSIAL storage file.
         */
        EXPAND,
        /**
         * Task to view content of a MUSIAL storage file.
         */
        VIEW,
        /**
         * Task to export sequence data from a MUSIAL storage file.
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
            // Record the start time of the program.
            startTime = System.currentTimeMillis();

            // Load metadata such as software name, version, and contact information.
            loadMetadata();

            // Check if any arguments were provided; if not, display usage information and exit.
            if (args.length == 0) {
                System.out.printf("No arguments were specified. Call `java -jar %s-%s.jar [-h|--help]` for more information.%n",
                        Musial.softwareName, Musial.softwareVersion);
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
                    Logging.logInfo("Execute task \033[1;1mbuild\033[0m");
                    Update.build();
                }
                case EXPAND -> {
                    Logging.logInfo("Execute task \033[1;1mexpand\033[0m");
                    Update.expand();
                }
                case VIEW -> {
                    Logging.logInfo("Execute task \033[1;1mview\033[0m");
                    View.run();
                }
                case SEQUENCE -> {
                    Logging.logInfo("Execute task \033[1;1msequence\033[0m");
                    Sequence.run();
                }
                default -> System.exit(-2); // Exit with an error code if the task is undefined.
            }
        } catch (Exception e) {
            // Log the error message and stack trace, then exit with an error code.
            Logging.logError(e.getMessage());
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
        Musial.softwareName = properties.getProperty("name");
        Musial.softwareVersion = properties.getProperty("version");
        Musial.softwareContact = properties.getProperty("contact");
        Musial.softwareLicense = properties.getProperty("license");
        // Print information to stdout.
        Logging.printSoftwareInfo();
    }

    /**
     * Provides functionality to update (build or expand) a MUSIAL storage.
     */
    public static class Update {

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
            if (storage.skipSnpEff()) {
                Logging.logInfo("Skipping SnpEff analysis as per user request.");
            } else if (storage.hasMissingContigSequences()) {
                Logging.logWarning("Skip SnpEff annotation; no reference sequence provided.");
            } else if (storage.getFeatures().stream().allMatch(f -> f.type.equals("region"))) {
                Logging.logWarning("Skip SnpEff annotation; all features are of type region.");
            } else if (!storage.hasNovelVariants()) {
                Logging.logWarning("Skip SnpEff annotation; no variants to annotate.");
            } else {
                Logging.logInfo("Run SnpEff annotation.");
                storage.annotateVariants();
            }

            // Infer sequence types if reference sequences are available.
            Logging.logInfo("Infer sequence types.");
            storage.updateSequenceTypes();

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
                outputPath += File.separator + "musial_storage_" + Logging.getDate().hashCode() + storageExtension;
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
            Storage storage = Storage.Factory.fromCLI();

            // Update the storage with variant calls, annotations, and statistics.
            update(storage);

            // Write the storage data to the specified output file.
            Logging.logInfo("Write storage to file " + CLI.parameters.get("output"));
            Storage.Factory.serialize(storage, outputFile);

            // Log summary information about the storage and execution time.
            Logging.logDone(
                    "Storage contains %d samples, %d features, %d variants. Processed %d genotypes. Execution time: %.2f seconds."
                            .formatted(
                                    storage.getSamples().size(),
                                    storage.getFeatures().size(),
                                    storage.getVariantsCount(),
                                    storage.getProcessedGenotypes(),
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
                            + "musial_storage_"
                            + Logging.getTimestamp().hashCode()
                            + storageExtension
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
                Storage.Factory.setSampleInformation(storage, new File(sampleInfoFile));
            }

            // Add VCF files to the storage for processing.
            //noinspection unchecked
            Storage.Factory.setVcfFiles(storage, (List<String>) CLI.parameters.get("vcfInput"));

            // Update the storage with new data, annotations, and statistics.
            update(storage);

            // Write the updated storage to the specified output file, if the write flag is enabled.
            if (write) {
                Logging.logInfo("Write storage to file " + outputFile);
                Storage.Factory.serialize(storage, outputFile);
            }

            // Log summary information about the expanded storage and execution time.
            Logging.logDone(
                    "Storage %s with %d samples, %d variants. Processed %d genotypes. Execution time: %.2f seconds."
                            .formatted(
                                    write ? "expanded" : "expandable", // Indicate whether the storage was expanded or just expandable.
                                    storage.getSamples().size() - originalSampleCount, // Number of new samples added.
                                    storage.getVariantsCount() - originalVariantsCount, // Number of new variants added.
                                    storage.getProcessedGenotypes(), // Total number of genotypes processed.
                                    (System.currentTimeMillis() - startTime) / 1000.0 // Total execution time in seconds.
                            )
            );
        }
    }

    /**
     * Provides functionality to generate and display various types of tables based on the data stored in the `Storage` object.
     * <p>
     * Supports filtering and formatting the output for different content types such as features, alleles, samples, and variants,
     * as well as sequence types and variant calls per sample in a matrix format.
     */
    public static class View {

        /**
         * Represents a table structure for storing and displaying data in a tabular format.
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
         * A set of supported content types for the `view` task.
         * <p>
         * This set defines the types of data that can be generated and displayed
         * in tabular format by the `view` functionality of the MUSIAL application.
         * Each string in the set corresponds to a specific type of content:
         * <ul>
         *     <li><b>feature</b>: Represents genomic features.</li>
         *     <li><b>allele</b>: Represents alleles associated with features.</li>
         *     <li><b>sample</b>: Represents sample-specific data.</li>
         *     <li><b>type</b>: Represents sequence types for samples.</li>
         *     <li><b>variant</b>: Represents genomic variants.</li>
         *     <li><b>call</b>: Represents variant calls for samples.</li>
         * </ul>
         */
        public static final Set<String> content = Set.of(
                "feature",
                "allele",
                "sample",
                "type",
                "variant",
                "call"
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
            if (!content.matches(String.join(Constants.PIPE, Musial.View.content))) {
                throw new MusialException("Content (-c) has to be one of %s, but %s was provided."
                        .formatted(String.join(", ", Musial.View.content), content));
            }

            // Initialize sets to store filters for features, samples, and positions.
            Set<String> features = new HashSet<>(), samples = new HashSet<>(), positions = new HashSet<>();
            // Parse the `confine` parameter to populate the filter sets.
            //noinspection unchecked
            for (String value : (Set<String>) CLI.parameters.get("filter")) {
                if (storage.hasFeature(value)) features.add(value);
                if (storage.hasSample(value)) samples.add(value);
                if (value.matches("\\d+")) positions.add(value);
            }

            // Retrieve and validate the output destination.
            String output = (String) CLI.parameters.get("output");
            File outputFile = null;
            if (!output.equals("stdout")) {
                // Ensure the output path is a file, not a directory.
                outputFile = new File(output);
                if (outputFile.isDirectory()) {
                    outputFile = new File(outputFile.getAbsolutePath()
                            + File.separator
                            + "musial_view_%s_".formatted(content)
                            + Logging.getTimestamp().hashCode()
                            + ".tsv"
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
                case "feature" -> featureTable(storage, features);
                case "allele" -> alleleTable(storage, features, samples);
                case "sample" -> sampleTable(storage, samples);
                case "type" -> typeMatrix(storage, samples, features);
                case "variant" -> variantTable(storage, positions, samples, features);
                case "call" -> callMatrix(storage, samples, positions);
                default -> throw new MusialException("Unknown task `view` content %s.".formatted(content));
            };

            // Handle the case where no entries match the filters.
            if (table.identifiers.isEmpty()) {
                Logging.logWarning("No entries to view. Check your filter parameter.");
            }
            // Output the table to the console if "stdout" is specified.
            else if ("stdout".equals(output)) {
                System.out.println(table);
            }
            // Write the table to the specified file.
            else {
                IO.writeFile(outputFile.toPath(), table.toString());
                Logging.logInfo("Write results to %s.".formatted(output));
            }

            // Log the completion of the task with the execution time.
            Logging.logDone("Execution time %.2f seconds.".formatted((System.currentTimeMillis() - startTime) / 1000.0));
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
                    Comparator.comparingInt(i -> storage.getFeature(i).start), Constants.EMPTY);

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
         * Generates a table containing information about alleles and proteoforms for genomic features.
         * <p>
         * This method creates a table with rows representing alleles and proteoforms associated with features.
         * The table can be filtered to include only specific features and samples based on the provided sets.
         *
         * @param storage         The {@link Storage} instance containing the features and alleles to be included in the table.
         * @param includeFeatures A set of feature names to include in the table. If empty, all features are included.
         * @param includeSamples  A set of sample names to include in the table. If empty, all samples are included.
         * @return A {@link Table} object containing the allele and proteoform information.
         */
        private static Table alleleTable(Storage storage, Set<String> includeFeatures, Set<String> includeSamples) {
            // Initialize the table with the header "feature    type    uid" and a comparator for sorting by feature start position.
            Table table = new Table("feature\ttype\tuid", storage.getFeatures().size(),
                    Comparator.comparingInt(i -> storage.getFeature(i.split(Constants.TAB)[0]).start), Constants.EMPTY);

            // Stream through the features in the storage, filtering based on the includeFeatures set.
            storage.getFeatures().stream()
                    .filter(feature -> includeFeatures.isEmpty() || includeFeatures.contains(feature.name))
                    .forEach(feature -> {
                        // Set to track included alleles for the current feature.
                        Set<String> includeAlleles = new HashSet<>();

                        // Process each allele of the feature.
                        feature.getAlleles().forEach(allele -> {
                            // Check if the allele should be included based on the includeSamples set.
                            if (includeSamples.isEmpty() || includeSamples.stream().anyMatch(allele::hasOccurrence)) {
                                includeAlleles.add(allele._uid);

                                // Create a list of tuples representing the allele's attributes.
                                List<Tuple<String, String>> items = new ArrayList<>(allele.getAttributes().entrySet().stream()
                                        .map(entry -> new Tuple<>(entry.getKey(), entry.getValue()))
                                        .toList());

                                // Add the occurrence information of the allele to the list.
                                items.add(new Tuple<>("samples", allele.occurrenceAsString()));

                                // Add the allele's information to the table.
                                table.addContent("%s\tallele\t%s".formatted(feature.name, allele._uid), items);
                            }
                        });

                        // Process proteoforms if proteoform inference is not skipped and the feature is coding.
                        if (storage.runProteoformInference() && feature.isCoding()) {
                            feature.getProteoforms().forEach(proteoform -> {
                                // Check if the proteoform should be included based on the included alleles.
                                if (includeAlleles.stream().anyMatch(proteoform::hasOccurrence)) {
                                    // Create a list of tuples representing the proteoform's attributes.
                                    List<Tuple<String, String>> items = proteoform.getAttributes().entrySet().stream()
                                            .map(entry -> new Tuple<>(entry.getKey(), entry.getValue()))
                                            .toList();

                                    // Add the proteoform's information to the table.
                                    table.addContent("%s\tproteoform\t%s".formatted(feature.name, proteoform._uid), items);
                                }
                            });
                        }
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
                    Comparator.naturalOrder(), Constants.EMPTY);

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
         * Generates a table containing information about sequence types for samples and features.
         * <p>
         * This method creates a table with rows representing alleles and proteoforms associated with features.
         * The table can be filtered to include only specific samples and features based on the provided sets.
         *
         * @param storage         The {@link Storage} instance containing the samples and features to be included in the table.
         * @param includedSamples A set of sample names to include in the table. If empty, all samples are included.
         * @param includeFeatures A set of feature names to include in the table. If empty, all features are included.
         * @return A {@link Table} object containing the sequence type information.
         */
        private static Table typeMatrix(Storage storage, Set<String> includedSamples, Set<String> includeFeatures) {
            // Initialize the table with the header "name   type" and a comparator for sorting by feature start position.
            Table table = new Table("name\ttype", includeFeatures.isEmpty() ? storage.getFeatures().size() : includeFeatures.size(),
                    Comparator.comparingInt(i -> storage.getFeature(i.split(Constants.TAB)[0]).start), Constants.synonymous);

            // Stream through the samples in the storage, filtering based on the includedSamples set.
            storage.getSamples().stream()
                    .filter(sample -> includedSamples.isEmpty() || includedSamples.contains(sample.name))
                    .forEach(sample ->
                            // Stream through the alleles of the sample, filtering based on the includeFeatures set.
                            sample.getAlleles().stream()
                                    .filter(allele -> includeFeatures.isEmpty() || includeFeatures.contains(allele.getKey()))
                                    .forEach(allele -> {
                                        // Add the allele information to the table.
                                        List<Tuple<String, String>> items = List.of(new Tuple<>(sample.name, allele.getValue()));
                                        table.addContent(allele.getKey() + "\tallele", items);

                                        // Retrieve the feature associated with the allele.
                                        Feature feature = storage.getFeature(allele.getKey());

                                        // Add proteoform information if the feature is coding and proteoform inference is not skipped.
                                        if (feature.isCoding() && storage.runProteoformInference()) {
                                            items = List.of(new Tuple<>(sample.name, feature.getAllele(allele.getValue()).getAttribute(Constants.$Allele_proteoform)));
                                            table.addContent(allele.getKey() + "\tproteoform", items);
                                        }
                                    })
                    );

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
                    Comparator.comparingInt(i -> Integer.parseInt(i.split(Constants.TAB)[1])), Constants.EMPTY);

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
                                    sample -> variantInfo.hasOccurrence(Constants.$Attributable_samplesOccurrence, sample));

                            // If the variant matches the feature and sample filters, add it to the table.
                            if (hasFeature && hasSample) {
                                // Create a list of tuples representing the variant's attributes.
                                List<Tuple<String, String>> items = new ArrayList<>(List.of(
                                        new Tuple<>("type", variantInfo.type.name())
                                ));

                                // Add additional attributes of the variant to the list.
                                variantInfo.getAttributes().forEach((key, value) -> items.add(new Tuple<>(key, value)));

                                // Add the occurrence information of the variant to the list.
                                items.add(new Tuple<>("samples", String.join(Constants.COMMA, variantInfo.getSampleOccurrence())));

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
                    Comparator.comparingInt(s -> Integer.parseInt(s.split("\t")[1])), Constants.DOT);

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
     * Provides functionality for exporting nucleotide (nt) or amino acid (aa) sequences from genomic features and contigs
     * in the MUSIAL application.
     * <p>
     * It supports exporting sequences in FASTA format with options for alignment, merging, and including conserved reference content.
     */
    private static class Sequence {

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
                Feature feature = storage.getFeature(featureName); // Retrieve the feature by name.
                Contig contig = storage.getContig(feature.contig); // Retrieve the contig associated with the feature.

                // Export nucleotide or amino acid sequences based on the task parameters.
                if (nt) {
                    exportNtSequences(feature, contig, sampleNames, conserved, merge, strip, reference);
                } else {
                    exportAaSequences(feature, contig, sampleNames, conserved, merge, strip, reference);
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
         * @param feature     The genomic feature for which sequences are exported.
         * @param contig      The contig containing the feature and its variants.
         * @param sampleNames A set of sample names to filter alleles for sequence generation.
         * @param conserved   If true, generates sequences with conserved reference content.
         * @param merge       If true, merges sequences for all samples into a single output.
         * @param strip       If true, un-aligns sequences by removing gaps.
         * @param reference   If true, includes the reference sequence in the output.
         * @throws IOException If an I/O error occurs during file writing.
         */
        private static void exportNtSequences(Feature feature, Contig contig, Set<String> sampleNames,
                                              boolean conserved, boolean merge, boolean strip, boolean reference) throws IOException {

            // Check if conserved sequences are requested but the contig lacks reference sequence information.
            if (conserved && !contig.hasSequence()) {
                Logging.logWarning("Skip feature %s as contig %s has no reference sequence information (incompatible with conserved export)."
                        .formatted(feature.name, contig.name));
                return;
            }

            // Collect allele UIDs that match the provided sample names.
            final Set<String> alleleUids = feature.getAlleles().stream()
                    .filter(allele -> sampleNames.stream().anyMatch(allele::hasOccurrence))
                    .map(allele -> allele._uid)
                    .collect(Collectors.toSet());
            if (alleleUids.isEmpty()) {
                Logging.logWarning("Skipping feature %s: Only sequence type is the reference sequence."
                        .formatted(feature.name));
                return;
            }

            // Retrieve variants associated with the selected alleles.
            ArrayList<Tuple<Integer, String>> variants = contig.getVariantsByAlleles(feature, alleleUids);

            // Retrieve the reference content if conserved sequences are requested.
            final char[] referenceContent;
            if (conserved) referenceContent = contig.getSubsequence(feature.start, feature.end).toCharArray();
            else referenceContent = null;

            // Map to store positional context for variants, defined as the reference content and maximal insertion length per position.
            HashMap<Integer, Tuple<String, Integer>> positionalContext = new HashMap<>();

            // Function to update the positional context from variant information.
            BiConsumer<Tuple<Integer, String>, Integer> updatePositionalContext = (context, insertionLength) ->
                    positionalContext.merge(context.a, new Tuple<>(context.b, insertionLength), (e1, e2) -> {
                        if (!Objects.equals(e2.a, Constants.EMPTY) && !Objects.equals(e1.a, Constants.EMPTY) && !Objects.equals(e1.a, e2.a)) {
                            Logging.logWarning("Reference content conflict at position %d (%s and %s).".formatted(context.a, e1.a, e2.a));
                        }
                        return new Tuple<>(e1.a.equals(Constants.EMPTY) ? e2.a : e1.a, Math.max(e1.b, e2.b));
                    });

            // Process each variant to populate the positional context.
            VariantInformation info;
            for (Tuple<Integer, String> variant : variants) {
                info = contig.getVariantInformation(variant.a, variant.b);
                char[] ref = info.reference.toCharArray();
                if (info.type.equals(VariantInformation.Type.SNV)) {
                    updatePositionalContext.accept(new Tuple<>(variant.a, info.reference), 0);
                } else if (info.type.equals(VariantInformation.Type.DELETION)) {
                    if (variant.b.charAt(0) != ref[0])
                        updatePositionalContext.accept(new Tuple<>(variant.a, String.valueOf(ref[0])), 0);
                    for (int i = 1; i < ref.length; i++) {
                        updatePositionalContext.accept(new Tuple<>(variant.a + i, String.valueOf(ref[i])), 0);
                    }
                } else if (info.type.equals(VariantInformation.Type.INSERTION)) {
                    int length = variant.b.length() - 1;
                    if (variant.b.charAt(0) != ref[0])
                        updatePositionalContext.accept(new Tuple<>(variant.a, String.valueOf(ref[0])), length);
                    else updatePositionalContext.accept(new Tuple<>(variant.a, Constants.EMPTY), length);
                }
            }

            // StringBuilder to construct the sequence content.
            StringBuilder content = new StringBuilder(conserved ? feature.end - feature.start + 1 : variants.size());

            // Function to resolve reference content for a given position.
            Consumer<Integer> resolveReference = position -> {
                if (positionalContext.containsKey(position)) {
                    Tuple<String, Integer> context = positionalContext.get(position);
                    String referenceBase = context.a.isEmpty() && conserved ? String.valueOf(referenceContent[position - feature.start]) : context.a;
                    content.append(SequenceOperations.padGaps(referenceBase, referenceBase.length() + context.b));
                } else if (conserved) {
                    content.append(referenceContent[position - feature.start]);
                }
            };

            // Write the sequences to a FASTA file.
            String fileName = String.format("%s_c%d_m%d_x%d_r%d.fasta",
                    feature.name,
                    conserved ? 1 : 0,
                    merge ? 1 : 0,
                    strip ? 1 : 0,
                    reference ? 1 : 0
            );
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputDirectory + File.separator + fileName, StandardCharsets.UTF_8))) {
                // Function to write a sequence to the file with a given header.
                Consumer<String> dump = header -> {
                    try {
                        String sequence = strip ? content.toString().replaceAll(Constants.gapString, Constants.EMPTY) : content.toString();
                        writer.write("%s%s\n".formatted(header, String.join("\n", Splitter.fixedLength(80).split(sequence))));
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                };

                // Write the reference sequence if requested.
                if (reference) {
                    IntStream.rangeClosed(feature.start, feature.end).forEach(resolveReference::accept);
                    dump.accept(">reference allelic_frequency=%s\n".formatted(feature.getAttribute("reference_frequency")));
                }

                // Write sequences for each allele.
                Feature.Allele allele;
                String alt;
                boolean positionConserved;
                int deletedPositions;
                Tuple<String, Integer> context;
                Tuple<String, Integer> nullContext = new Tuple<>(Constants.EMPTY, 0);
                for (String alleleUid : alleleUids) {
                    allele = feature.getAllele(alleleUid);
                    content.setLength(0);
                    deletedPositions = 0;
                    for (int position = feature.start; position <= feature.end; position++) {
                        if (deletedPositions > 0) {
                            context = positionalContext.get(position);
                            content.append(SequenceOperations.padGaps(Constants.gapString, 1 + context.b));
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
                                if (positionConserved && !conserved) {
                                    content.append(SequenceOperations.padGaps(Constants.EMPTY, context.b));
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
                        dump.accept(allele.getFastaHeader(feature.name, alleleUid));
                    } else {
                        for (String sampleName : allele.getOccurrence()) {
                            allele.getFastaHeader(feature.name, sampleName);
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
         * @param feature     The genomic feature for which sequences are exported.
         * @param contig      The contig containing the feature and its variants.
         * @param sampleNames A set of sample names to filter proteoforms for sequence generation.
         * @param conserved   If true, generates sequences with conserved reference content.
         * @param merge       If true, merges sequences for all samples into a single output.
         * @param strip       If true, un-aligns sequences by removing gaps.
         * @param reference   If true, includes the reference sequence in the output.
         * @throws IOException     If an I/O error occurs during file writing.
         * @throws MusialException If a MUSIAL-specific error occurs.
         */
        private static void exportAaSequences(Feature feature, Contig contig, Set<String> sampleNames,
                                              boolean conserved, boolean merge, boolean strip, boolean reference) throws IOException, MusialException {
            // Check if the contig has reference sequence information; required for amino acid export.
            if (!contig.hasSequence()) {
                Logging.logWarning("Skip feature %s; contig %s has to have reference sequence information for amino acid export."
                        .formatted(feature.name, contig.name));
                return;
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

            // Translate the reference nucleotide sequence to amino acid sequence.
            final char[] referenceContent = SequenceOperations
                    .translateSequence(contig.getSubsequence(feature.start, feature.end), feature.isReverse()).toCharArray();

            // Collect all variants associated with the selected proteoforms.
            Set<Tuple<Integer, String>> variantsSet = new HashSet<>();
            for (String proteoformUid : proteoformUids) {
                Feature.Proteoform proteoform = feature.getProteoform(proteoformUid);
                proteoform.getVariants().forEach((key, value) -> variantsSet.add(new Tuple<>(key, value)));
            }
            ArrayList<Tuple<Integer, String>> variants = new ArrayList<>(variantsSet);
            variants.sort(Comparator.comparingInt(i -> i.a));

            // Map to store positional context for variants.
            HashMap<Integer, Tuple<String, Integer>> positionalContext = new HashMap<>();

            // Function to update the positional context with variant information.
            BiConsumer<Tuple<Integer, String>, Integer> updatePositionalContext = (context, insertionLength) ->
                    positionalContext.merge(context.a, new Tuple<>(context.b, insertionLength), (e1, e2) -> {
                        if (!Objects.equals(e2.a, Constants.EMPTY) && !Objects.equals(e1.a, Constants.EMPTY) && !Objects.equals(e1.a, e2.a)) {
                            Logging.logWarning("Reference content conflict at position %d (%s and %s).".formatted(context.a, e1.a, e2.a));
                        }
                        return new Tuple<>(e1.a.equals(Constants.EMPTY) ? e2.a : e1.a, Math.max(e1.b, e2.b));
                    });

            // Process each variant to populate the positional context.
            for (Tuple<Integer, String> variant : variants) {
                if (VariantInformation.isSubstitution(variant.b)) {
                    updatePositionalContext.accept(new Tuple<>(variant.a, String.valueOf(referenceContent[variant.a - 1])), 0);
                } else {
                    boolean mixed = variant.b.charAt(0) != referenceContent[variant.a - 1];
                    if (VariantInformation.isDeletion(variant.b)) {
                        if (mixed) {
                            updatePositionalContext.accept(new Tuple<>(variant.a, String.valueOf(referenceContent[variant.a - 1])), 0);
                        }
                        for (int i = 1; i < variant.b.length(); i++) {
                            updatePositionalContext.accept(new Tuple<>(variant.a + i, String.valueOf(referenceContent[variant.a - 1 + i])), 0);
                        }
                    } else if (VariantInformation.isInsertion(variant.b)) {
                        int length = variant.b.length() - 1;
                        if (mixed) {
                            updatePositionalContext.accept(new Tuple<>(variant.a, String.valueOf(referenceContent[variant.a - 1])), length);
                        } else {
                            updatePositionalContext.accept(new Tuple<>(variant.a, Constants.EMPTY), length);
                        }
                    }
                }
            }

            // StringBuilder to construct the sequence content.
            StringBuilder content = new StringBuilder(conserved ? Math.ceilDiv((feature.end - feature.start + 1), 3) : variants.size());

            // Function to resolve reference content for a given position.
            Consumer<Integer> resolveReference = position -> {
                if (positionalContext.containsKey(position)) {
                    Tuple<String, Integer> context = positionalContext.get(position);
                    String referenceBase = context.a.isEmpty() && conserved ? String.valueOf(referenceContent[position - 1]) : context.a;
                    content.append(SequenceOperations.padGaps(referenceBase, referenceBase.length() + context.b));
                } else if (conserved) {
                    content.append(referenceContent[position - 1]);
                }
            };

            // Write the sequences to a FASTA file.
            String fileName = String.format("%s_c%d_m%d_x%d_r%d.fasta",
                    feature.name,
                    conserved ? 1 : 0,
                    merge ? 1 : 0,
                    strip ? 1 : 0,
                    reference ? 1 : 0
            );
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputDirectory + File.separator + fileName, StandardCharsets.UTF_8))) {

                // Function to write a sequence to the file with a given header.
                Consumer<String> dump = header -> {
                    try {
                        String sequence = strip ? content.toString().replaceAll(Constants.gapString, Constants.EMPTY) : content.toString();
                        writer.write("%s%s\n".formatted(header, String.join("\n", Splitter.fixedLength(80).split(sequence))));
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                };

                // Write the reference sequence if requested.
                if (reference) {
                    IntStream.rangeClosed(feature.start, feature.end).forEach(resolveReference::accept);
                    dump.accept(">reference allelic_frequency=%s\n".formatted(feature.getAttribute("reference_frequency")));
                }

                // Write sequences for each proteoform.
                for (String proteoformUid : proteoformUids) {
                    Feature.Proteoform proteoform = feature.getProteoform(proteoformUid);
                    content.setLength(0);
                    int deletedPositions = 0;
                    for (int position = 1; position <= referenceContent.length; position++) {
                        if (deletedPositions > 0) {
                            Tuple<String, Integer> context = positionalContext.get(position);
                            content.append(SequenceOperations.padGaps(Constants.gapString, 1 + context.b));
                            deletedPositions--;
                            if (proteoform.hasVariant(position))
                                Logging.logWarning("Conflict with variant %s at deleted position %d for allele %s of feature %s."
                                        .formatted(proteoform.getVariant(position), position, proteoformUid, feature.name));
                        } else if (proteoform.hasVariant(position)) {
                            String alt = proteoform.getVariant(position);
                            Tuple<String, Integer> context = positionalContext.getOrDefault(position, new Tuple<>(Constants.EMPTY, 0));
                            boolean positionConserved = context.a.isEmpty();
                            if (VariantInformation.isSubstitution(alt)) {
                                content.append(SequenceOperations.padGaps(alt, 1 + context.b));
                            } else if (VariantInformation.isDeletion(alt)) {
                                if (positionConserved && !conserved) {
                                    content.append(SequenceOperations.padGaps(Constants.EMPTY, context.b));
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
                        dump.accept(proteoform.getFastaHeader(feature.name, proteoformUid));
                    } else {
                        for (String alleleUid : proteoform.getOccurrence()) {
                            for (String sampleName : feature.getAllele(alleleUid).getOccurrence()) {
                                proteoform.getFastaHeader(feature.name, sampleName);
                            }
                        }
                    }
                }
            }
        }

    }

}