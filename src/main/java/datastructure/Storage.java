package datastructure;

import com.google.gson.GsonBuilder;
import com.google.gson.internal.LinkedTreeMap;
import exceptions.MusialException;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Tuple;
import htsjdk.variant.utils.BinomialCoefficientUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import main.Musial;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.RandomStringUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import utility.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * The `Storage` class is a central component of the project, designed to manage genomic data, including contigs, features, and samples.
 * <p>
 * It provides methods for adding, retrieving, and processing genomic information, as well as handling variant calls and annotations.
 * The class is structured to support efficient storage and manipulation of data, leveraging Java collections and utility classes.
 * It integrates external tools like SnpEff for annotation and ensures data integrity through validation and compliance with Sequence Ontology rules.
 * <p>
 * <b>Core Data Structures:</b>
 * <ul>
 * <li>Contigs: Stored in a `Map (String, {@link Contig})`, contigs represent chromosomes or plasmids. Each contig can store its sequence and associated variants.</li>
 * <li>Features: Stored in a `Map (String, {@link Feature})`, features represent genomic elements like genes or mRNA. These are validated and processed using Sequence Ontology (SO) terms.</li>
 * <li>Samples: Stored in a `Map (String, {@link Sample}), samples represent variant calls from distinct biological samples. Metadata and variant calls are associated with each sample.</li>
 * </ul>
 * <b>Key Functionalities:</b>
 * <ul>
 * <li>Adding and Managing Contigs: Methods like `addContig` and `addContigIfAbsent` allow adding contigs with sequences or as placeholders.</li>
 * <li>Handling Features: Features are processed and validated using Sequence Ontology terms, with methods like `transferFeatureInformation` and `validateFeatures`.</li>
 * <li>Managing Samples and Variants: Samples are added using `addSample`, and variant calls are managed with methods like `addVariantCallToSample`.</li>
 * <li>Variant Processing and Annotation: The `updateVariants` method processes VCF files to extract variant information, while `runSnpEffAnalysis` integrates with SnpEff to annotate novel variants.</li>
 * </ul>
 */
public class Storage {

    /**
     * Represents the parameters used for configuring the storage of variant data.
     * <p>
     * This record encapsulates various parameters that control the behavior of the storage system,
     * including thresholds for coverage and frequency, as well as exclusions for specific positions
     * and variants. These parameters are immutable once set.
     *
     * @param minimalCoverage         The minimal coverage of a variant call to be accepted. Must be greater than or equal to 0.
     *                                This ensures that only variant calls with sufficient read depth are considered.
     * @param minimalFrequency        The minimal frequency of a variant call to be accepted. Must be between 0.0 and 1.0, inclusive.
     *                                This parameter filters out low-frequency variants that may be due to sequencing errors.
     * @param storeFiltered           Whether to retain filtered calls as ambiguous bases (N). If true, filtered calls are retained
     *                                in the storage, allowing downstream analysis to consider them as ambiguous data.
     * @param skipSnpEff              Whether to skip the SnpEff annotation process. If true, the annotation step is bypassed,
     *                                which can save time if annotation is not required.
     * @param skipProteoformInference Whether to skip proteoform inference. If true, the inference of proteoforms
     *                                (protein isoforms) is not performed, which can be useful for non-coding regions.
     * @param excludedPositions       A map associating contig names with sets of positions to exclude from storage. Cannot be null
     *                                but can be empty. This allows specific genomic positions to be ignored during analysis.
     * @param excludedVariants        A map associating contig names with sets of alternative variants to exclude from storage. Cannot
     *                                be null but can be empty. This allows specific variants to be excluded from consideration.
     */
    private record Parameters(
            double minimalCoverage, // Minimum read depth required for a variant call to be accepted.
            double minimalFrequency, // Minimum allele frequency required for a variant call to be accepted.
            boolean storeFiltered, // Flag to determine whether filtered calls should be retained as ambiguous bases.
            boolean skipSnpEff, // Flag to determine whether SnpEff annotation should be skipped.
            boolean skipProteoformInference, // Flag to determine whether proteoform inference should be skipped.
            Map<String, Set<Integer>> excludedPositions, // Map of contig names to sets of positions to exclude from analysis.
            Map<String, Set<String>> excludedVariants // Map of contig names to sets of alternative variants to exclude.
    ) {
        // Compact constructor with validation logic omitted for simplicity.
    }

    /**
     * Static parameters used by this storage.
     * <p>
     * This field holds an instance of {@link Parameters}, which contains the configuration
     * for the storage system. The parameters are immutable and define the behavior of the
     * storage, such as thresholds and exclusions.
     */
    private final Parameters parameters;

    /**
     * Contigs/chromosomes/plasmids of the reference sequence.
     */
    private final Map<String, Contig> contigs;

    /**
     * Maintained genomic features.
     */
    private final Map<String, Feature> features;

    /**
     * Individual samples, i.e., variant calls from one distinct biological sample.
     */
    private final Map<String, Sample> samples;

    /**
     * Transient list of novel variants. <i>This is automatically filled during variant call processing.</i>
     */
    private transient ArrayList<Triple<String, Integer, String>> novelVariants = new ArrayList<>(128);

    /**
     * Transient sample information. This is only relevant for the `build` task. <i>This should only be set by the
     * {@link Factory#samplesFromCLI} method during initialization of a storage.</i>
     */
    private transient Map<String, Map<String, String>> sampleInfo = new LinkedTreeMap<>();

    /**
     * Transient list of VCF files to collect data from to expand this storage. <i>This should only be set by the
     * {@link Factory#samplesFromCLI} method during initialization of a storage.</i>
     */
    private transient List<File> vcfFiles = new ArrayList<>();

    /**
     * Transient accessor to (indexed) reference sequences. <i>This should only be set by the
     * {@link Factory} class during initialization of a storage.</i>
     */
    private transient ReferenceSequenceFile reference = null;

    /**
     * Transient accessor to the VCF handler.
     */
    private transient VcfHandler vcfHandler = new VcfHandler();

    /**
     * Map of sequence ontology (SO) terms and their respective hierarchy levels as used by MUSIAL.
     * TODO: Optional extension to support UTRs, etc.?
     */
    public final static Map<String, Integer> SO = new HashMap<>() {{
        put("region", 0);
        put("gene", 1);
        put("pseudogene", 1);
        put("mRNA", 2);
        put("tRNA", 2);
        put("rRNA", 2);
        put("tmRNA", 2);
        put("ncRNA", 2);
        put("SRP_RNA", 2);
        put("RNase_P_RNA", 2);
        put("CDS", 3);
        put("exon", 3);
    }};

    /**
     * Factory class for creating and managing instances of {@link Storage}.
     * <p>
     * This class provides methods to load storage from command line interface (CLI) parameters,
     * files, and to save storage to files. It also handles the initialization of transient properties
     * and manages the loading of reference sequences, features, and sample information.
     */
    public static final class Factory {

        /**
         * Initializes a {@link Storage} from CLI parameters.
         * <p>
         * This method initializes a new {@link Storage} instance using parameters from the command line interface (CLI).
         * It also loads reference sequences, features, and sample information from the CLI parameters.
         *
         * @return A {@link Storage} object representing the loaded data.
         * @throws MusialException If an error occurs while loading or validating the data.
         * @throws IOException     If an error occurs while reading files or parsing data.
         */
        public static Storage fromCLI() throws MusialException, IOException {
            if (CLI.parameters.isEmpty()) {
                throw new MusialException("CLI parameters are empty; run CLI.parse() before initializing the storage.");
            }
            Parameters parameters = parametersFromCLI();
            Storage storage = new Storage(parameters);
            referenceFromCLI(storage);
            featuresFromCLI(storage);
            samplesFromCLI(storage);
            storage.validateFeatures();
            return storage;
        }

        /**
         * Initializes a {@link Storage} from a file.
         * <p>
         * This method loads the storage from a specified file in JSON format. It handles both compressed
         * and uncompressed files. The method also initializes transient properties after loading.
         *
         * @param file The file to load the storage from.
         * @return A {@link Storage} object representing the loaded data.
         * @throws IOException If an error occurs while reading the file or parsing the JSON data.
         */
        public static Storage deserialize(File file) throws IOException {
            Validation.checkFile(file);

            Function<BufferedReader, Storage> storageFromReader = reader ->
                    new GsonBuilder().setPrettyPrinting().create().fromJson(reader, Storage.class);

            try (BufferedReader bufferedReader = new BufferedReader(
                    new InputStreamReader(file.getAbsolutePath().endsWith(".gz")
                            ? new GZIPInputStream(Files.newInputStream(file.toPath()))
                            : Files.newInputStream(file.toPath())))) {
                Storage storage = storageFromReader.apply(bufferedReader);
                storage.setTransientProperties();
                return storage;
            } catch (IOException e) {
                throw new IOException("Failed to load MUSIAL storage from file %s; %s"
                        .formatted(file.getAbsolutePath(), e.getMessage()));
            }
        }

        /**
         * Writes the given `Storage` object to a file in JSON format.
         * <p>
         * This method ensures that the file has the correct extension (`.json` or `.json.gz` for GZIP-compressed files),
         * converts the `Storage` object to a JSON string, and writes it to the specified file. If the file path ends
         * with `.gz`, the JSON data is compressed using GZIP before writing.
         *
         * @param storage The `Storage` object to be serialized and written to the file.
         * @param file    The `File` object representing the target file.
         * @throws IOException If an error occurs during file operations, such as writing or compression.
         */
        public static void serialize(Storage storage, File file) throws IOException {
            // Ensure the file has the correct extension
            if (!(file.getAbsolutePath().endsWith(".json") || file.getAbsolutePath().endsWith(".json.gz"))) {
                file = new File(file.getAbsolutePath() + Musial.storageExtension);
            }

            // Convert the storage object to a JSON string
            String jsonData = new GsonBuilder().setPrettyPrinting().create().toJson(storage, Storage.class);

            // Write the JSON data to the file, using GZIP if necessary
            try (Writer writer = file.getAbsolutePath().endsWith(".gz")
                    ? new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file)))
                    : new FileWriter(file, StandardCharsets.UTF_8)) {
                writer.write(jsonData);
            } catch (IOException e) {
                // Throw a new IOException with a detailed error message if writing fails
                throw new IOException(String.format("Failed to write MUSIAL storage to file %s; %s.", file.getAbsolutePath(), e.getMessage()));
            }
        }

        /**
         * Sets sample information in the storage from a specified file.
         * <p>
         * This method reads a tabular file containing sample metadata and populates the `sampleInfo` map
         * in the provided `Storage` instance. The file is expected to be in a format that can be parsed
         * into a nested map structure, where the outer map keys represent sample names and the inner map
         * contains attribute-value pairs for each sample.
         *
         * @param storage The {@link Storage} instance where the sample information will be stored.
         * @param file    The {@link File} object representing the input file containing sample metadata.
         * @throws IOException If an error occurs while reading or parsing the file.
         */
        public static void setSampleInformation(Storage storage, File file) throws IOException {
            Validation.checkFile(file);
            storage.sampleInfo.putAll(IO.readTabularFileAsNestedMap(file));
        }

        /**
         * Sets the VCF files in the storage from a list of file paths.
         * <p>
         * This method processes a list of file paths, adding valid VCF files or directories containing VCF files
         * to the storage's `vcfFiles` list. It supports both uncompressed `.vcf` files and compressed `.vcf.gz` files.
         * If no valid VCF files are found, an exception is thrown.
         *
         * @param storage The {@link Storage} instance where the VCF files will be stored.
         * @param paths   A list of file paths to process. Each path can be a file or a directory.
         * @throws MusialException If no valid VCF files are found in the provided paths.
         * @throws IOException     If an error occurs while accessing the file system.
         */
        public static void setVcfFiles(Storage storage, List<String> paths) throws MusialException, IOException {
            for (String pathName : paths) {
                Path path = Paths.get(pathName);
                Function<Path, Boolean> isVcf = p -> (p.toString().endsWith(FileExtensions.VCF) ||
                        p.toString().endsWith(FileExtensions.COMPRESSED_VCF));

                if (Files.isDirectory(path)) {
                    try (Stream<Path> stream = Files.list(path)) {
                        List<File> files = stream.filter(isVcf::apply).map(Path::toFile).toList();
                        for (File file : files) {
                            Validation.checkFile(file);
                            storage.vcfFiles.add(file);
                        }
                    }
                } else if (isVcf.apply(path)) {
                    File file = path.toFile();
                    Validation.checkFile(file);
                    storage.vcfFiles.add(path.toFile());
                }
            }

            // Throw an exception if no valid VCF files were added to the storage.
            if (storage.vcfFiles.isEmpty()) {
                throw new MusialException("Expect at least one VCF file, but none was provided.");
            }
        }

        /**
         * Loads parameters from CLI.
         * <p>
         * This method loads parameters from the command line interface (CLI) and validates them. It checks for
         * the presence of required parameters and sets default values for optional ones. The method also handles
         * exceptions related to invalid parameter values.
         *
         * @return A {@link Parameters} object containing the loaded and validated parameters.
         * @throws IOException If an error occurs while reading or validating the parameters for excluded positions/variants.
         */
        private static Parameters parametersFromCLI() throws IOException {
            double minimalCoverage = 3.0; // Default value
            if (CLI.parameters.containsKey("minimalCoverage")) {
                String value = String.valueOf(CLI.parameters.get("minimalCoverage"));
                if (Validation.isPositiveDouble(value)) {
                    minimalCoverage = Double.parseDouble(value);
                } else {
                    Logging.logWarning("Invalid value for `minimalCoverage`; expected a positive integer. Defaulting to 3.");
                }
            } else {
                Logging.logWarning("No value for `minimalCoverage` specified; defaulting to 3.");
            }

            double minimalFrequency = 0.65; // Default value
            if (CLI.parameters.containsKey("minimalFrequency")) {
                String value = String.valueOf(CLI.parameters.get("minimalFrequency"));
                if (Validation.isPercentage(value)) {
                    minimalFrequency = Double.parseDouble(value);
                } else {
                    Logging.logWarning("Invalid value for `minimalFrequency`; expected a percentage between 0.0 and 1.0. Defaulting to 0.65.");
                }
            } else {
                Logging.logWarning("No value for `minimalFrequency` specified; defaulting to 0.65.");
            }

            boolean storeFiltered = true; // Default value
            if (CLI.parameters.containsKey("storeFiltered")) {
                String value = String.valueOf(CLI.parameters.get("storeFiltered"));
                if ("true".equalsIgnoreCase(value) || "false".equalsIgnoreCase(value)) {
                    storeFiltered = Boolean.parseBoolean(value);
                } else {
                    Logging.logWarning("Invalid value for `storeFiltered`; expected `true` or `false`. Defaulting to true.");
                }
            } else {
                Logging.logWarning("No value for `storeFiltered` specified; defaulting to true.");
            }

            boolean skipSnpEff = false; // Default value
            if (CLI.parameters.containsKey("skipSnpEff")) {
                String value = String.valueOf(CLI.parameters.get("skipSnpEff"));
                if ("true".equalsIgnoreCase(value) || "false".equalsIgnoreCase(value)) {
                    skipSnpEff = Boolean.parseBoolean(value);
                } else {
                    Logging.logWarning("Invalid value for `skipSnpEff`; expected `true` or `false`. Defaulting to false.");
                }
            } else {
                Logging.logWarning("No value for `skipSnpEff` specified; defaulting to false.");
            }

            boolean skipProteoformInference = true; // Default value
            if (CLI.parameters.containsKey("skipProteoformInference")) {
                String value = String.valueOf(CLI.parameters.get("skipProteoformInference"));
                if ("true".equalsIgnoreCase(value) || "false".equalsIgnoreCase(value)) {
                    skipProteoformInference = Boolean.parseBoolean(value);
                } else {
                    Logging.logWarning("Invalid value for `skipProteoformInference`; expected `true` or `false`. Defaulting to true.");
                }
            } else {
                Logging.logWarning("No value for `skipProteoformInference` specified; defaulting to true.");
            }

            Map<String, Set<Integer>> excludedPositions = excludedPositionsFromCLI();
            Map<String, Set<String>> excludedVariants = excludedVariantsFromCLI();

            return new Parameters(minimalCoverage, minimalFrequency, storeFiltered, skipSnpEff, skipProteoformInference, excludedPositions,
                    excludedVariants);
        }

        /**
         * Loads excluded positions from CLI parameters.
         * <p>
         * This method reads a file containing excluded positions and populates a map where the key is the contig name
         * and the value is a set of excluded positions. The file is expected to have rows with at least three columns:
         * the contig name, the start position, and the end position. The positions are inclusive.
         * <p>
         * If the file path is null, blank, or the file is invalid, the method returns an empty map.
         * If a line has fewer than three columns, an IOException is thrown.
         *
         * @return A map where the key is the contig name and the value is a set of excluded positions.
         * @throws IOException If the file is invalid or a line has fewer than three columns.
         */
        private static Map<String, Set<Integer>> excludedPositionsFromCLI() throws IOException {
            // Initialize a map to store excluded positions, with the contig name as the key.
            Map<String, Set<Integer>> excludedPositions = new LinkedHashMap<>();

            String path = (String) CLI.parameters.get("excludedPositions");
            // Check if the path is null or blank; if so, return an empty map.
            if (path != null && !path.isBlank()) {
                File file = new File(path);
                // Validate the file's existence and readability.
                Validation.checkFile(file);
                // Read the file content into a list of strings.
                List<String> content = IO.readFile(file);
                // Detect the separator used in the file (e.g., tab or comma).
                String separator = IO.detectSeparator(content);
                if (separator.isEmpty()) {
                    throw new IOException("Invalid field separator in file %s. Use \\t (tab) or , (comma)."
                            .formatted(file.getAbsolutePath()));
                }

                // Iterate over each line in the file.
                for (String line : content) {
                    // Skip lines that start with a comment sign.
                    if (!line.startsWith(Constants.SIGN)) {
                        String[] parts = line.split(separator);
                        if (parts.length < 3) {
                            throw new IOException("Invalid number of columns in row %d of file %s."
                                    .formatted(content.indexOf(line) + 1, file.getAbsolutePath()));
                        }
                        // Add the range of positions to the map for the corresponding contig.
                        excludedPositions
                                .computeIfAbsent(parts[0], k -> new HashSet<>()) // Create a new set if the contig is not already in the map.
                                .addAll(IntStream.rangeClosed(Integer.parseInt(parts[1]), Integer.parseInt(parts[2])) // Generate a range of positions.
                                        .boxed() // Box the primitive int values into Integer objects.
                                        .collect(Collectors.toSet())); // Collect the range into a set.
                    }
                }
            }

            if (!excludedPositions.isEmpty()) {
                Logging.logInfo("Set %d positions as excluded.".formatted(excludedPositions.values().stream().mapToInt(Set::size).sum()));
            }

            // Return the populated map of excluded positions.
            return excludedPositions;
        }

        /**
         * Loads excluded variants from CLI parameters.
         * <p>
         * This method reads a file containing excluded variants and populates a map where the key is the contig name
         * and the value is a set of excluded variants. Each excluded variant is represented as a string in the format:
         * "position:reference:alternate".
         * <p>
         * If the file path is null, blank, or the file is invalid, the method returns an empty map.
         *
         * @return A map where the key is the contig name and the value is a set of excluded variants.
         * @throws IOException If the file is invalid or cannot be read.
         */
        private static Map<String, Set<String>> excludedVariantsFromCLI() throws IOException {
            // Initialize a map to store excluded variants, with the contig name as the key.
            Map<String, Set<String>> excludedVariants = new LinkedHashMap<>();

            String path = (String) CLI.parameters.get("excludedVariants");
            // Check if the path is null or blank; if so, return an empty map.
            if (path != null && !path.isBlank()) {
                File file = new File(path);
                // Validate the file's existence and readability.
                Validation.checkFile(file);
                // Read the file content into a list of strings.
                List<String> content = IO.readFile(file);
                // Detect the separator used in the file (e.g., tab or comma).
                String separator = IO.detectSeparator(content);
                if (separator.isEmpty()) {
                    throw new IOException("Invalid field separator in file %s. Use \\t (tab) or , (comma)."
                            .formatted(file.getAbsolutePath()));
                }

                // Iterate over each line in the file.
                for (String line : content) {
                    // Skip lines that start with a comment sign.
                    if (!line.startsWith(Constants.SIGN)) {
                        String[] parts = line.split(separator);
                        if (parts.length < 4) {
                            throw new IOException("Invalid number of columns in row %d of file %s."
                                    .formatted(content.indexOf(line) + 1, file.getAbsolutePath()));
                        }
                        String chrom = parts[0]; // Contig name.
                        String pos = parts[1];   // Position.
                        String ref = parts[2];   // Reference base.
                        String alt = parts[3];   // Alternate base.
                        // Add the excluded variant to the map for the corresponding contig.
                        excludedVariants.computeIfAbsent(chrom, k -> new HashSet<>())
                                .add(pos + Constants.COLON + ref + Constants.COLON + alt);
                    }
                }
            }

            if (!excludedVariants.isEmpty()) {
                Logging.logInfo("Set %d variants as excluded.".formatted(excludedVariants.values().stream().mapToInt(Set::size).sum()));
            }

            // Return the populated map of excluded variants.
            return excludedVariants;
        }

        /**
         * Loads sample information, i.e. metadata and variant information, from CLI parameters.
         * <p>
         * This method loads sample information from the command line interface (CLI) parameters. See the
         * {@link #setSampleInformation} and {@link #setVcfFiles} methods for details on
         * how the sample information is loaded.
         *
         * @param storage The {@link Storage} instance to load sample information into.
         * @throws IOException     If an error occurs while reading the sample metadata file or VCF files.
         * @throws MusialException If an error occurs while validating the sample information or VCF files.
         */
        private static void samplesFromCLI(Storage storage) throws IOException, MusialException {
            // Parse sample metadata file if available.
            String path = (String) CLI.parameters.get("vcfMeta");
            if (path != null && !path.isBlank()) {
                File file = new File(path);
                // File validation is handled in the setSampleInformation method.
                setSampleInformation(storage, file);
            }

            // Collect VCF file paths from the specified list.
            //noinspection unchecked
            List<String> paths = (List<String>) CLI.parameters.get("vcfInput");
            // File validation is handled in the setVcfFiles method.
            setVcfFiles(storage, paths);
        }

        /**
         * Loads the reference sequence from CLI parameters.
         * <p>
         * This method loads the reference sequence from the command line interface (CLI) parameters. It checks
         * if the reference sequence file is specified and validates its accessibility. If the file is valid,
         * it creates an {@link IndexedFastaSequenceFile} instance and sets it as the reference for the storage.
         *
         * @param storage The {@link Storage} instance to load the reference into.
         * @throws MusialException If an error occurs while loading or parsing the reference sequence.
         * @throws IOException     If an error occurs while reading the reference file or creating the index.
         */
        private static void referenceFromCLI(Storage storage) throws MusialException, IOException {
            String path = (String) CLI.parameters.get("reference");
            if (path != null && !path.isBlank()) {
                File file = new File(path);
                Validation.checkFile(file);
                try {
                    storage.setReference(new IndexedFastaSequenceFile(file.toPath(),
                            FastaSequenceIndexCreator.buildFromFasta(file.toPath())));
                } catch (IOException e) {
                    throw new MusialException("Failed to create index for reference sequence %s; %s"
                            .formatted(file.getAbsolutePath(), e.getMessage()));
                }
            }
        }

        /**
         * Loads features from CLI parameters.
         * <p>
         * This method loads features from the command line interface (CLI) parameters. It supports loading
         * and validating reference features from a file.
         *
         * @param storage The {@link Storage} instance to load features into.
         * @throws MusialException If an error occurs while loading or parsing the features.
         * @throws IOException     If an error occurs while reading the feature file or parsing the data.
         */
        private static void featuresFromCLI(Storage storage) throws MusialException, IOException {
            // Reference annotation requires a reference sequence.
            if (!storage.hasReference() && CLI.parameters.containsKey("annotation")) {
                throw new MusialException("Annotation (GFF3) specified without a sequence (FASTA).");
            }

            // Features require a reference annotation.
            if (CLI.parameters.containsKey("features") && !CLI.parameters.containsKey("annotation")) {
                throw new MusialException("Features specified without an annotation (GFF3).");
            }

            // Read reference features file and parse available reference features, if specified.
            FeatureList featureList = new FeatureList();
            if (CLI.parameters.containsKey("annotation")) {
                String path = (String) CLI.parameters.get("annotation");
                if (path != null && !path.isBlank()) {
                    File file = new File(path);
                    Validation.checkFile(file);
                    featureList = GFF3Reader.read(file.getCanonicalPath());
                }
            }

            // Helper to process attributes
            Consumer<Map<String, String>> reprocessAttributes = attributes -> attributes.replaceAll((k, v) -> Arrays.stream(v.split(Constants.COMMA))
                    .map(s -> s.split("\\|")[0])
                    .filter(s -> !s.isEmpty())
                    .collect(Collectors.joining(Constants.COMMA)));

            // Process reference features.
            String path = (String) CLI.parameters.get("features");
            if (path != null && !path.isBlank()) {
                File file = new File(path);
                Validation.checkFile(file);

                List<String> content = IO.readFile(file);
                String separator = IO.detectSeparator(content);
                if (separator.isEmpty()) {
                    throw new IOException("Invalid field separator in file %s. Use \\t (tab) or , (comma)."
                            .formatted(file.getAbsolutePath()));
                }

                for (String line : content) {
                    if (line.startsWith(Constants.SIGN)) continue;

                    String[] parts = line.split(separator);
                    if (parts.length < 2) {
                        throw new IOException("Invalid row format in reference features file %s."
                                .formatted(file.getAbsolutePath()));
                    }

                    String key = parts[0];
                    String value = parts[1];
                    Map<String, String> customAttributes = new HashMap<>();
                    if (parts.length > 2) {
                        Arrays.stream(parts[2].split(Constants.SEMICOLON))
                                .map(entry -> entry.split(Constants.EQUAL))
                                .filter(info -> info.length == 2)
                                .forEach(info -> customAttributes.put(info[0], info[1]));
                    }

                    FeatureList matchedFeatures = featureList.selectByAttribute(key, value);
                    if (matchedFeatures.isEmpty()) {
                        Logging.logWarning("No match for value %s under key %s in reference features file."
                                .formatted(value, key));
                        continue;
                    }

                    for (FeatureI matchedFeature : matchedFeatures) {
                        Map<String, String> attributes = matchedFeature.getAttributes();
                        reprocessAttributes.accept(attributes);
                        storage.addFeature(matchedFeature, value, attributes);
                    }

                    reprocessAttributes.accept(customAttributes);
                    storage.getFeature(value).addAttributesIfAbsent(customAttributes);
                }
            }

            if (!CLI.parameters.containsKey("features") && !featureList.isEmpty()) {
                Logging.logInfo("No explicit features are specified; read all annotated features from file %s."
                        .formatted(CLI.parameters.get("annotation")));
                for (FeatureI featureI : featureList) {
                    if (!"region".equals(featureI.type())) {
                        Map<String, String> attributes = featureI.getAttributes();
                        reprocessAttributes.accept(attributes);
                        String name = attributes.getOrDefault("Name", "feature_" + (storage.features.size() + 1));
                        storage.addFeature(featureI, name, attributes);
                    }
                }
            }

            // Parse features from reference sequence if none are available
            if (storage.features.isEmpty() && storage.hasReference()) {
                Logging.logWarning("No features were specified; one feature of type region per reference contig is derived.");
                ReferenceSequence entry = storage.reference.nextSequence();
                while (entry != null) {
                    Map<String, String> attributes = new HashMap<>(1);
                    attributes.put("ID", "%s:%d..%d".formatted(entry.getName(), 1, entry.length()));
                    storage.addFeature(entry.getName(), entry.getName(), 1, entry.length(), '+', "region", attributes);
                    entry = storage.reference.nextSequence();
                }
                storage.reference.reset();
            }
        }
    }

    /**
     * Constructs a new {@link Storage} instance with the specified parameters.
     * <p>
     * This constructor initializes the {@link Storage} object with the provided configuration
     * parameters. It also initializes empty containers for contigs, features, and samples
     * using {@link LinkedTreeMap}, ensuring that the data is stored in a sorted and efficient manner.
     *
     * @param parameters The {@link Parameters} object containing the configuration for the storage.
     *                   This includes thresholds, exclusions, and other settings for managing
     *                   genomic data.
     */
    private Storage(Parameters parameters) {
        this.parameters = parameters;
        this.contigs = new LinkedTreeMap<>(); // Initialize empty contig container.
        this.features = new LinkedTreeMap<>(); // Initialize empty feature container.
        this.samples = new LinkedTreeMap<>(); // Initialize empty sample container.
    }

    /**
     * Retrieve a collection of sequence ontology terms for a given level.
     *
     * @param level The level to retrieve the sequence ontology terms for.
     * @return A collection of sequence ontology terms for the specified level.
     */
    public static Collection<String> getSOTerms(int level) {
        return SO.entrySet().stream().filter(e -> e.getValue().equals(level)).map(Map.Entry::getKey).collect(Collectors.toSet());
    }

    /**
     * Return whether {@link #reference} is set.
     *
     * @return Get the reference sequence.
     */
    public boolean hasReference() {
        return Objects.nonNull(this.reference);
    }

    /**
     * Retrieves the minimum coverage parameter.
     *
     * @return Get the minimum coverage parameter.
     */
    public double minimumCoverage() {
        return this.parameters.minimalCoverage;
    }

    /**
     * Retrieves the minimum frequency parameter.
     *
     * @return Get the minimum frequency parameter.
     */
    public double minimumFrequency() {
        return this.parameters.minimalFrequency;
    }

    /**
     * Retrieves whether filtered calls should be retained as ambiguous bases (N).
     *
     * @return {@code true} if filtered calls are retained as ambiguous bases, {@code false} otherwise.
     */
    public boolean storeFiltered() {
        return this.parameters.storeFiltered;
    }

    /**
     * Checks if the SnpEff annotation process should be skipped.
     *
     * @return {@code true} if the SnpEff annotation process is skipped, {@code false} otherwise.
     */
    public boolean skipSnpEff() {
        return this.parameters.skipSnpEff;
    }

    /**
     * Checks if proteoform inference should be run.
     *
     * @return {@code true} if proteoform inference should be run, {@code false} otherwise.
     */
    public boolean runProteoformInference() {
        return !this.parameters.skipProteoformInference;
    }

    /**
     * Whether {@code position} is excluded on {@code contig}.
     *
     * @param contig   Contig (name) to check for exclusion.
     * @param position Position to check for exclusion.
     * @return True if {@code position} on {@code contig} is excluded from analysis.
     */
    public boolean isPositionExcluded(String contig, int position) {
        return this.parameters.excludedPositions.containsKey(contig) && this.parameters.excludedPositions.get(contig).contains(position);
    }

    /**
     * Whether {@code variant} is excluded on {@code contig} at {@code position}.
     *
     * @param contig    Contig (name) to check for exclusion.
     * @param position  Position to check for exclusion.
     * @param reference Reference base at {@code position}.
     * @param variant   Variant base at {@code position}.
     * @return True if {@code variant} on {@code contig} at {@code position} is excluded from analysis.
     */
    public boolean isVariantExcluded(String contig, int position, String reference, String variant) {
        return this.parameters.excludedVariants.containsKey(contig) &&
                this.parameters.excludedVariants.get(contig).contains(reference + Constants.COLON + position + Constants.COLON + variant);
    }

    /**
     * Retrieves the number of processed genotype records.
     * <p>
     * This method returns the value of the `processedGenotypes` field from the `VcfHandler` class.
     * The field keeps track of the total number of genotype records that have been processed
     * during the analysis of VCF files.
     *
     * @return The total number of processed genotype records as a {@code long}.
     */
    public long getProcessedGenotypes() {
        return this.vcfHandler.processedGenotypes;
    }

    /**
     * Adds a contig to the storage with the specified name and sequence.
     * <p>
     * This method compresses the provided sequence using GZIP and calculates its length.
     * If the sequence is null or empty, it assigns an empty string as the compressed sequence
     * and sets the length to 0. The contig is then added to the storage with its attributes.
     *
     * @param name     The name of the contig to add.
     * @param sequence The sequence of the contig. Can be null or empty.
     * @throws IOException If an error occurs during sequence compression.
     */
    public void addContig(String name, String sequence) throws IOException {
        if (!hasContig(name)) {
            String compressedSequence;
            int length;
            if (Objects.nonNull(sequence) && !sequence.isEmpty()) {
                compressedSequence = IO.gzipCompress(sequence); // Compress the sequence using GZIP.
                length = sequence.length(); // Calculate the length of the sequence.
            } else {
                compressedSequence = Constants.EMPTY; // Assign an empty string if the sequence is null or empty.
                length = 0; // Set the length to 0 for an empty sequence.
            }
            this.contigs.put(name, new Contig(name, compressedSequence)); // Add the contig to the storage.
            this.contigs.get(name).addAttributeIfAbsent(Constants.$Contig_length, String.valueOf(length)); // Add the length attribute.
        }
    }

    /**
     * Retrieves a contig by its name.
     * <p>
     * This method searches for a contig in the storage by its name. If the contig exists,
     * it returns the corresponding {@link Contig} object. If the contig does not exist,
     * it returns {@code null}.
     *
     * @param name The name of the contig to retrieve.
     * @return The {@link Contig} object associated with the specified name, or {@code null} if no such contig exists.
     */
    public Contig getContig(String name) {
        return this.contigs.getOrDefault(name, null);
    }

    /**
     * Query whether a contig is stored in this instance by its name.
     *
     * @param name The name of the contig.
     * @return True if a contig is stored for {@code name}.
     */
    public boolean hasContig(String name) {
        return this.contigs.containsKey(name);
    }

    /**
     * Checks if all contigs in the storage have a non-empty sequence.
     * <p>
     * This method iterates through all contigs in the storage and checks if all of them have
     * a sequence that is not empty.
     *
     * @return {@code true} if any contig has an empty sequence, {@code false} otherwise.
     */
    public boolean hasMissingContigSequences() {
        return !this.contigs.values().stream().allMatch(Contig::hasSequence);
    }

    /**
     * Returns a collection view of the contigs stored in the storage.
     *
     * @return Collection of stored {@link #contigs}.
     */
    public Collection<Contig> getContigs() {
        return this.contigs.values();
    }

    /**
     * Returns the total number of variants stored.
     *
     * @return Number of all stored {@link Contig#variants} across all {@link #contigs}.
     */
    public long getVariantsCount() {
        return contigs.values().stream().mapToLong(Contig::getVariantsCount).sum();
    }

    /**
     * Checks if there are any novel variants stored in the storage.
     * <p>
     * This method verifies whether the `novelVariants` list contains any entries.
     * Novel variants are those that have been identified during variant call processing
     * but are not yet annotated or processed further.
     *
     * @return {@code true} if there are novel variants in the storage, {@code false} otherwise.
     */
    public boolean hasNovelVariants() {
        return !this.novelVariants.isEmpty();
    }

    /**
     * Retrieve a {@link Feature} by its name.
     *
     * @param name The name of the feature.
     * @return The queried feature or null, if no feature is stored with {@code name}.
     */
    public Feature getFeature(String name) {
        return this.features.getOrDefault(name, null);
    }

    /**
     * Returns a collection view of the features stored in the storage.
     *
     * @return Collection of {@link #features}.
     */
    public Collection<Feature> getFeatures() {
        return this.features.values();
    }

    /**
     * Query whether a feature is stored in this instance by its name.
     *
     * @param name The name of the feature.
     * @return True if a feature is stored for {@code name}.
     */
    public boolean hasFeature(String name) {
        return this.features.containsKey(name);
    }

    /**
     * Retrieve a {@link Sample} by its name.
     *
     * @param name The name of the sample.
     * @return The queried sample.
     */
    public Sample getSample(String name) {
        return this.samples.get(name);
    }

    /**
     * Returns a collection view of the samples stored in the storage.
     *
     * @return Collection of {@link #samples}.
     */
    public Collection<Sample> getSamples() {
        return this.samples.values();
    }

    /**
     * Retrieves a collection of samples that need to be updated based on the presence of variant records.
     * <p>
     * This method filters the samples stored in the `samples` map and returns only those samples
     * whose names are present as keys in the `vcfAnalysis.records` map. These samples are considered
     * to have associated variant records and require updates.
     *
     * @return A collection of {@link Sample} objects that need to be updated.
     */
    public Collection<Sample> getSamplesToUpdate() {
        return samples.values().stream()
                .filter(sample -> this.vcfHandler.records.containsKey(sample.name))
                .collect(Collectors.toList());
    }

    /**
     * Query whether a sample is stored in this instance by its name.
     *
     * @param name The name of the sample.
     * @return True if a sample is stored for {@code name}.
     */
    public boolean hasSample(String name) {
        return this.samples.containsKey(name);
    }

    /**
     * Updates the variants in the storage by processing VCF files and transferring the relevant information.
     * <p>
     * This method performs the following steps:
     * <ul>
     *   <li>Clears the existing variant records in the {@link VcfHandler}.</li>
     *   <li>Analyzes each VCF file in the {@code vcfFiles} list to extract variant information.</li>
     *   <li>Clears the {@code vcfFiles} list after processing.</li>
     *   <li>Transfers sample information from the processed VCF records to the storage.</li>
     *   <li>Transfers sample attributes from the {@code sampleInfo} map to the corresponding samples in the storage.</li>
     *   <li>Transfers variant information from the sample variant calls to the storage.</li>
     * </ul>
     *
     * @throws IOException If an error occurs while analyzing VCF files or transferring data.
     */
    public void updateVariants() throws IOException {
        vcfHandler.clearRecords(); // Clear existing variant records.

        if (this.features.isEmpty()) {
            Logging.logWarning("No features available; infer features from variant calls.");
            for (Map.Entry<String, Integer> entry : vcfHandler.inferContigs().entrySet()) {
                addContig(entry.getKey(), Constants.EMPTY);
                String id = "%s:%d..%d".formatted(entry.getKey(), 1, entry.getValue());
                Map<String, String> attributes = new HashMap<>(1);
                attributes.put("ID", id);
                addFeature(entry.getKey(), entry.getKey(), 1, entry.getValue(), '+', "region", attributes);
            }
        }

        for (File file : vcfFiles) {
            vcfHandler.analyzeVcf(file); // Analyze each VCF file to extract variant information.
        }
        vcfFiles.clear(); // Clear the list of VCF files after processing.
        transferSampleInformation(); // Transfer sample information to the storage.
        transferSampleAttributes(); // Transfer sample attributes to the storage.
        transferVariantsInformation(); // Transfer variant information to the storage.
    }

    /**
     * Annotates novel variants in the storage using SnpEff.
     * <p>
     * This method invokes the SnpEff analysis process to annotate novel variants stored in the `Storage` instance.
     * It delegates the annotation task to the `runSnpEffAnalysis` method of the `VcfHandler` class.
     * The results of the annotation are integrated back into the storage.
     *
     * @throws IOException     If an error occurs during file operations required for the SnpEff analysis.
     * @throws MusialException If the SnpEff annotation process encounters an error.
     */
    public void annotateVariants() throws IOException, MusialException {
        this.vcfHandler.runSnpEffAnalysis();
    }

    /**
     * Updates sequence types for all samples and features in the storage.
     * <p>
     * This method iterates through all samples that need to be updated and all features in the storage.
     * For each feature, it retrieves the associated contig and filters the variants for the sample
     * within the feature's start and end positions. The filtered variants are reduced to a map
     * containing the variant positions and their corresponding alternative allele base strings.
     * <p>
     * If the filtered variants are not empty, the method updates the allele for the feature using
     * the contig, variants, and sample. If the feature is coding and the contig has a sequence,
     * the proteoform for the feature is also updated.
     * <p>
     * Finally, the method performs HDBSCAN clustering for alleles and proteoforms per feature,
     * updating their attributes with the clustering results.
     *
     * @throws IOException     If an error occurs during sequence processing.
     * @throws MusialException If an error occurs during allele or proteoform updates.
     */
    public void updateSequenceTypes() throws IOException, MusialException {
        // Iterate through all samples that need to be updated.
        for (Sample sample : getSamplesToUpdate()) {
            // Iterate through all features in the storage.
            for (Feature feature : getFeatures()) {
                // Retrieve the contig associated with the feature.
                Contig contig = getContig(feature.contig);
                // Filter variants for the sample within the feature's start and end positions.
                ArrayList<Tuple<Integer, String>> variants = contig.getVariantsBySampleAndLocation(sample.name, feature.start, feature.end);
                if (variants.isEmpty()) continue; // Skip if no variants are found.
                // Update allele information for the feature with respect to the sample.
                String alleleUid = feature.updateAllele(contig, variants, sample);
                // If the feature is coding and the contig has a sequence, update the proteoform.
                if (feature.isCoding() && contig.hasSequence() && runProteoformInference()) {
                    feature.updateProteoform(contig, alleleUid);
                }
            }
        }

        // Perform clustering for alleles and proteoforms per feature.
        for (Feature feature : getFeatures()) {
            // Reset clustering and add alleles to the dataset.
            Clustering.reset();
            feature.getAlleles().forEach(allele ->
                    Clustering.addToDataset(allele._uid, allele.getVariants())
            );
            // Train the clustering model.
            Clustering.train();

            // Process clustering results for alleles.
            Clustering.getClusteringResult().forEach(entry -> {
                // Retrieve the allele by its name and update its attributes with clustering results.
                Feature.Allele allele = feature.getAllele(entry.name());
                allele.setAttribute(Constants.$SequenceType_clusterLabel, String.valueOf(entry.label()));
                allele.setAttribute(Constants.$SequenceType_clusterIndex, String.valueOf(entry.idx()));
                allele.setAttribute(Constants.$SequenceType_clusterOutlierScore, String.format(Locale.US, "%.3f", entry.outlierScore()));
                // Update the allele's name to include clustering information.
                allele.setName("%s.a%s.%s".formatted(feature.name, entry.label(), entry.idx()));
            });

            // If the feature is coding, process proteoforms.
            if (feature.isCoding() && runProteoformInference()) {
                // Reset clustering and add proteoforms to the dataset.
                Clustering.reset();
                feature.getProteoforms().forEach(proteoform ->
                        Clustering.addToDataset(proteoform._uid, proteoform.getVariants())
                );
                // Train the clustering model.
                Clustering.train();

                // Process clustering results for proteoforms.
                Clustering.getClusteringResult().forEach(entry -> {
                    // Retrieve the proteoform by its name and update its attributes with clustering results.
                    Feature.Proteoform proteoform = feature.getProteoform(entry.name());
                    proteoform.setAttribute(Constants.$SequenceType_clusterLabel, String.valueOf(entry.label()));
                    proteoform.setAttribute(Constants.$SequenceType_clusterIndex, String.valueOf(entry.idx()));
                    proteoform.setAttribute(Constants.$SequenceType_clusterOutlierScore, String.format(Locale.US, "%.3f", entry.outlierScore()));
                    // Update the proteoform's name to include clustering information.
                    proteoform.setName("%s.p%s.%s".formatted(feature.name, entry.label(), entry.idx()));
                });
            }
        }
    }

    /**
     * Updates various statistics for samples, contigs, and features in the storage.
     * <p>
     * This method performs the following tasks:
     * <ul>
     *   <li>Calculates and updates statistics for each sample, including total calls, filtered calls, mean coverage, and mean quality.</li>
     *   <li>Determines the frequency of disrupted and modified proteoforms for coding features in each sample.</li>
     *   <li>Updates variant frequency for each contig and aggregates substitution and indel counts per sample.</li>
     *   <li>Calculates and updates allelic frequencies, reference frequencies, and proteoform statistics for each feature.</li>
     * </ul>
     */
    public void updateStatistics() {
        // Initialize decimal formats for frequency and general statistics.
        DecimalFormat frequencyFormat = new DecimalFormat("0.##E0", DecimalFormatSymbols.getInstance(Locale.US));
        DecimalFormat decimalFormat = new DecimalFormat("##.#", DecimalFormatSymbols.getInstance(Locale.US));

        // Count the number of coding features in the storage.
        float noCodingFeatures = features.values().stream().filter(Feature::isCoding).count();

        // Initialize maps to store substitution and indel counts per sample.
        Map<String, Integer> perSampleSubstitutions = new HashMap<>();
        Map<String, Integer> perSampleInDels = new HashMap<>();

        // Iterate over each sample to calculate and update sample-specific statistics.
        for (Sample sample : samples.values()) {
            int totalCalls = 0, filteredCalls = 0;
            List<Integer> coverages = new ArrayList<>(), qualities = new ArrayList<>();
            perSampleSubstitutions.put(sample.name, 0);
            perSampleInDels.put(sample.name, 0);

            // Process variant calls for the sample to calculate coverage and quality statistics.
            for (Map<Integer, String> variantCalls : sample.variantCalls.values()) {
                totalCalls += variantCalls.size();
                for (String variantCall : variantCalls.values()) {
                    String[] callParts = variantCall.split(Constants.SEMICOLON);
                    coverages.add(Integer.parseInt(callParts[1]));
                    if (variantCall.startsWith(Constants.lowCoverageCallPrefix) || variantCall.startsWith(Constants.lowFrequencyCallPrefix)) {
                        filteredCalls++;
                    } else {
                        qualities.add(Integer.parseInt(callParts[2]));
                    }
                }
            }

            // Update sample attributes with calculated statistics.
            sample.setAttribute("#calls", String.valueOf(totalCalls));
            sample.setAttribute("#filtered", String.valueOf(filteredCalls));
            sample.setAttribute("mean_coverage", decimalFormat.format(coverages.stream().mapToInt(Integer::intValue).average().orElse(0)));
            sample.setAttribute("mean_quality", decimalFormat.format(qualities.stream().mapToInt(Integer::intValue).average().orElse(0)));
            sample.setAttribute("reference_frequency", frequencyFormat.format(1 - (sample.alleles.size() / (float) features.size())));

            // Calculate proteoform statistics for coding features if proteoform inference is not skipped.
            if (!parameters.skipProteoformInference()) {
                int disrupted = 0, modified = 0;
                for (var alleleEntry : sample.alleles.entrySet()) {
                    Feature feature = getFeature(alleleEntry.getKey());
                    if (feature.isCoding()) {
                        String proteoformUid = feature.getAllele(alleleEntry.getValue()).getAttribute(Constants.$Allele_proteoform);
                        if (!Objects.equals(proteoformUid, Constants.synonymous)) {
                            Collection<String> effects = feature.getProteoform(proteoformUid).getAttributeAsCollection(Constants.$SequenceType_effects);
                            if (effects.contains("start_lost") || effects.contains("stop_gained")) {
                                disrupted++;
                            } else {
                                modified++;
                            }
                        }
                    }
                }
                sample.setAttribute("proteoform_disrupted_frequency", frequencyFormat.format(disrupted / noCodingFeatures));
                sample.setAttribute("proteoform_modified_frequency", frequencyFormat.format(modified / noCodingFeatures));
            }
        }

        // Update variant frequency and aggregate substitution/indel counts for each contig.
        for (Contig contig : contigs.values()) {
            contig.variants.forEach((position, innerMap) -> innerMap.forEach((altBases, variantInfo) -> {
                int sampleCount = variantInfo.getSampleOccurrence().size();
                variantInfo.setAttribute("variant_frequency", frequencyFormat.format((float) sampleCount / samples.size()));

                variantInfo.getSampleOccurrence().forEach(sampleName -> {
                    Map<String, Integer> targetMap = switch (variantInfo.type) {
                        case SNV -> perSampleSubstitutions;
                        case INSERTION, DELETION -> perSampleInDels;
                    };
                    targetMap.put(sampleName, targetMap.get(sampleName) + 1);
                });
            }));
        }

        // Update sample attributes with aggregated substitution and indel counts.
        perSampleSubstitutions.forEach((sample, count) -> samples.get(sample).setAttribute("#substitutions", String.valueOf(count)));
        perSampleInDels.forEach((sample, count) -> samples.get(sample).setAttribute("#indels", String.valueOf(count)));

        // Initialize a map to store proteoform occurrences for each feature.
        Map<String, Integer> perProteoformOccurrence = new HashMap<>();

        // Iterate over each feature to calculate and update feature-specific statistics.
        for (Feature feature : features.values()) {
            float nonReferenceOccurrence = 0;
            int proteoformsModified = 0, proteoformsDisrupted = 0;
            perProteoformOccurrence.clear();

            // Process alleles for the feature to calculate allelic frequencies and proteoform statistics.
            for (SequenceType allele : feature.getAlleles()) {
                int alleleOccurrence = allele.getOccurrence().size();
                allele.setAttribute("allelic_frequency", frequencyFormat.format(alleleOccurrence / (float) samples.size()));
                nonReferenceOccurrence += alleleOccurrence;

                if (!parameters.skipProteoformInference() && feature.isCoding()) {
                    String proteoformUid = allele.getAttribute(Constants.$Allele_proteoform);
                    if (!Objects.equals(proteoformUid, Constants.synonymous)) {
                        Collection<String> effects = feature.getProteoform(proteoformUid).getAttributeAsCollection(Constants.$SequenceType_effects);
                        if (effects.contains("start_lost") || effects.contains("stop_gained")) {
                            proteoformsDisrupted++;
                        } else {
                            proteoformsModified++;
                        }
                        perProteoformOccurrence.merge(proteoformUid, alleleOccurrence, Integer::sum);
                    }
                }
            }

            // Update feature attributes with calculated statistics.
            feature.setAttribute("reference_frequency", frequencyFormat.format(1 - (nonReferenceOccurrence / samples.size())));
            feature.setAttribute("#alleles", decimalFormat.format(feature.alleles.size()));

            if (!parameters.skipProteoformInference() && feature.isCoding()) {
                int proteoformCount = feature.getProteoforms().size();
                float productDisruptedFrequency = proteoformCount == 0 ? 0 : proteoformsDisrupted / (float) proteoformCount;
                float productModifiedFrequency = proteoformCount == 0 ? 0 : proteoformsModified / (float) proteoformCount;
                feature.setAttribute("product_disrupted_frequency", frequencyFormat.format(productDisruptedFrequency));
                feature.setAttribute("product_modified_frequency", frequencyFormat.format(productModifiedFrequency));
                feature.setAttribute("#proteoforms", decimalFormat.format(feature.proteoforms.size()));
                perProteoformOccurrence.forEach((proteoformUid, count) ->
                        feature.getProteoform(proteoformUid).setAttribute("allelic_frequency", frequencyFormat.format(count / (float) samples.size()))
                );
            }
        }
    }

    /**
     * Initializes transient properties of the {@link Storage} instance.
     * <p>
     * This method initializes transient properties such as {@link #novelVariants}, {@link #sampleInfo},
     * and {@link #vcfFiles}. It also ensures that the contigs have their sequence caches initialized.
     */
    private void setTransientProperties() {
        if (this.vcfHandler == null)
            this.vcfHandler = new VcfHandler();
        if (this.novelVariants == null)
            this.novelVariants = new ArrayList<>(64);
        if (this.sampleInfo == null)
            this.sampleInfo = new LinkedTreeMap<>();
        if (this.vcfFiles == null)
            this.vcfFiles = new ArrayList<>();
        contigs.values().forEach(contig -> {
            if (contig.sequenceCache == null) contig.sequenceCache = new HashMap<>();
        });
    }

    /**
     * Set the reference to use to the passed {@link IndexedFastaSequenceFile}.
     *
     * @param indexedFastaSequenceFile Instance of {@link IndexedFastaSequenceFile}. Can also be null.
     * @throws IOException If reference sequence compression fails.
     */
    private void setReference(IndexedFastaSequenceFile indexedFastaSequenceFile) throws IOException {
        this.reference = indexedFastaSequenceFile;
        ReferenceSequence referenceSequence = this.reference.nextSequence();
        while (Objects.nonNull(referenceSequence)) {
            if (!this.hasContig(referenceSequence.getName())) {
                this.addContig(referenceSequence.getName(), referenceSequence.getBaseString());
            }
            referenceSequence = this.reference.nextSequence();
        }
        this.reference.reset();
    }

    /**
     * Transfers feature information from a {@link FeatureI} object to the storage.
     * <p>
     * This method extracts the necessary details from the provided {@link FeatureI} object, such as
     * the sequence name, start and end positions, strand, and type, and delegates the processing
     * to the overloaded {@code transferFeatureInformation} method.
     *
     * @param featureI   The {@link FeatureI} object containing the feature information to transfer.
     * @param name       The name of the feature.
     * @param attributes A map of attributes associated with the feature.
     */
    private void addFeature(FeatureI featureI, String name, Map<String, String> attributes) {
        addFeature(name, featureI.seqname(), featureI.location().bioStart(), featureI.location().bioEnd(), featureI.location().bioStrand(),
                featureI.type(), attributes);
    }

    /**
     * Transfers feature information to the storage.
     * <p>
     * This method processes and validates the provided feature information, including its type, location, and attributes.
     * It checks if the feature type is supported by the Sequence Ontology (SO) map and determines a unique identifier (UID)
     * for the feature. If a feature with the same UID already exists, it validates compatibility with the parent feature
     * and updates the "children" attribute if applicable. Otherwise, it creates a new feature and adds it to the storage.
     * Processed attributes are removed from the attributes map, and the remaining attributes are extended for the feature.
     *
     * @param name       The name of the feature.
     * @param chrom      The chromosome or contig name where the feature is located.
     * @param start      The start position of the feature.
     * @param end        The end position of the feature.
     * @param strand     The strand of the feature ('+' or '-').
     * @param type       The type of the feature (e.g., "gene", "mRNA").
     * @param attributes A map of attributes associated with the feature.
     */
    private void addFeature(String name, String chrom, Number start, Number end, char strand, String type, Map<String, String> attributes) {
        // Validate the feature type against the Sequence Ontology (SO) map.
        if (!SO.containsKey(type)) {
            if (Logging.logDump.add("INVALID_FEATURE_TYPE_%s".formatted(type))) {
                Logging.logError("Features of type %s are currently not supported and will be ignored (%s)."
                        .formatted(type, name));
            }
            return;
        }

        // Validate the feature location.
        if ((int) start >= (int) end || (int) start < 1) {
            Logging.logError("Feature %s has an invalid location (%s:%d..%d) and will be ignored."
                    .formatted(name, chrom, (int) start, (int) end));
            return;
        }

        // Validate the compatibility with stored contigs.
        if (!this.hasContig(chrom)) {
            Logging.logError("Feature %s specifies an unknown parent locus (%s) and will be ignored."
                    .formatted(name, chrom));
            return;
        } else {
            int chromLength = Integer.parseInt(this.getContig(chrom).getAttribute(Constants.$Contig_length));
            if (chromLength != 0 && chromLength < (int) end) {
                Logging.logError("Feature %s specifies a location (%s:%d..%d) that exceeds the length of its parent locus (%s)."
                        .formatted(name, chrom, (int) start, (int) end, chromLength));
                return;
            }
        }

        // Determine the unique identifier (UID) for the feature.
        final String uid;
        if (attributes.containsKey("Parent")) {
            uid = attributes.get("Parent").matches("^.*-.*$") ? attributes.get("Parent").split("-")[1] : attributes.get("Parent");
        } else if (attributes.containsKey("ID")) {
            uid = attributes.get("ID").matches("^.*-.*$") ? attributes.get("ID").split("-")[1] : attributes.get("ID");
        } else {
            uid = attributes.containsKey("locus_tag") ? attributes.get("locus_tag") : "%s:%d..%d".formatted(chrom, (int) start, (int) end);
        }

        // Check if a feature with the same UID already exists.
        Optional<Feature> optionalFeature = this.features.values().stream().filter(feature -> feature._uid.equals(uid)).findFirst();
        Feature feature;

        if (optionalFeature.isPresent()) {
            feature = optionalFeature.get();

            // Validate compatibility with the parent feature if the "Parent" attribute is present.
            if (attributes.containsKey("Parent") && attributes.get("Parent").matches("^.*-%s$".formatted(uid))) {
                if ((int) start < feature.start || (int) end > feature.end || !feature.contig.equals(chrom) || feature.strand != strand) {
                    Logging.logError("Feature update failed; feature %s with UID %s has an incompatible location with its parent feature %s."
                            .formatted(feature.name, feature._uid, name));
                    return;
                }

                // Add or extend the "children" attribute for the feature.
                if (!feature.hasAttribute("children")) {
                    feature.setAttribute("children", "%s:%d:%d".formatted(type, (int) start, (int) end));
                } else {
                    feature.extendAttribute("children", "%s:%d:%d".formatted(type, (int) start, (int) end));
                }
            } else {
                Logging.logError("Feature update failed; feature %s with UID %s already exists, but is not the parent of feature %s."
                        .formatted(feature.name, feature._uid, name));
                return;
            }
        } else {
            // Create a new feature if it does not already exist.
            feature = new Feature(name, chrom, start, end, strand, type, uid);
            feature.setAttribute("children", ""); // Ensure the "children" attribute is initialized.
            this.features.put(name, feature);
        }

        // Remove processed attributes from the attributes map.
        attributes.remove("ID");
        attributes.remove("Parent");
        attributes.remove("Name");
        attributes.remove("old_locus_tag");
        attributes.remove("gbkey");

        // Extend the feature's attributes with the remaining attributes.
        feature.extendAttributes(attributes);
    }

    /**
     * Validates the features stored in the `features` map to ensure compliance with Sequence Ontology (SO) rules.
     * <p>
     * This method performs the following validations and adjustments:
     * <ul>
     *     <li>Removes children for features with level 0 SO terms.</li>
     *     <li>Ensures only one level 1 SO term exists for a feature or its children.</li>
     *     <li>Ensures only one level 2 SO term exists for a feature or its children.</li>
     *     <li>Adjusts feature type and location for features with lower-level SO terms.</li>
     *     <li>Imputes missing children for specific SO term relationships (e.g., CDS to mRNA).</li>
     * </ul>
     * Features that fail validation are either adjusted or removed from the `features` map.
     */
    private void validateFeatures() {
        for (String featureName : new ArrayList<>(features.keySet())) {
            Feature feature = features.get(featureName);
            SortedMap<String, List<Tuple<Integer, Integer>>> children = feature.getChildren();

            // Remove children for level 0 SO term types.
            if (SO.get(feature.type) == 0 && !children.isEmpty()) {
                feature.setAttribute("children", "");
                Logging.logWarningOnce("REMOVE_SO0_CHILDREN",
                        "Features of type(s) %s are not supported to have children and associated children of %s will be removed."
                                .formatted(String.join(", ", getSOTerms(0)), feature.name));
                return;
            }

            // Ensure only one level 1 SO term exists.
            if (countSOLevel(children, feature, 1) > 1) {
                Logging.logWarning("Only one of %s is allowed as the type of the feature or its children; %s is removed."
                        .formatted(String.join(", ", getSOTerms(1)), feature.name));
                features.remove(featureName);
            }

            // Ensure only one level 2 SO term exists.
            if (countSOLevel(children, feature, 2) > 1) {
                Logging.logWarning("Only one of %s is allowed as the type of the feature or its children; %s is removed."
                        .formatted(String.join(", ", getSOTerms(2)), feature.name));
                features.remove(featureName);
            }

            // Adjust feature type and location for lower-level SO terms.
            if (SO.get(feature.type) > 1)
                adjustFeatureToGene(feature, children);

            // Impute missing children.
            inferFeatureChild(children, feature, "CDS", "mRNA");
            inferFeatureChild(children, feature, "mRNA", "CDS");
            for (String soTerm : getSOTerms(2)) {
                if (!"mRNA".equals(soTerm)) {
                    inferFeatureChild(children, feature, soTerm, "exon");
                }
            }
        }
    }

    /**
     * Counts the number of Sequence Ontology (SO) terms at a specific hierarchy level
     * for a given feature and its children.
     * <p>
     * This method calculates the total count of SO terms at the specified level by:
     * <ul>
     *     <li>Counting the children whose SO term matches the specified level.</li>
     *     <li>Adding 1 if the feature's own SO term matches the specified level.</li>
     * </ul>
     *
     * @param children A sorted map of child features, where the key is the SO term and the value is a list of position ranges.
     * @param feature  The {@link Feature} object whose SO term is being evaluated.
     * @param level    The hierarchy level of the SO terms to count.
     * @return The total count of SO terms at the specified level for the feature and its children.
     */
    private int countSOLevel(SortedMap<String, List<Tuple<Integer, Integer>>> children, Feature feature, int level) {
        return (int) children.keySet().stream().filter(k -> SO.get(k) == level).count() + (SO.get(feature.type) == level ? 1 : 0);
    }

    /**
     * Adjusts a feature to represent a "gene" type by updating its type, location, and children.
     * <p>
     * This method recalculates the start and end positions of the feature based on its children's ranges,
     * removes children with level 1 Sequence Ontology (SO) terms, and updates the feature's type to "gene".
     * The updated feature is then stored back into the `features` map.
     *
     * @param feature  The {@link Feature} object to adjust.
     * @param children A sorted map of child features, where the key is the SO term and the value is a list of position ranges.
     */
    private void adjustFeatureToGene(Feature feature, SortedMap<String, List<Tuple<Integer, Integer>>> children) {
        // Calculate the new start position as the minimum of the children's start positions and the feature's current start.
        int start = Math.min(children.values().stream()
                .flatMapToInt(ls -> ls.stream().mapToInt(t -> t.a))
                .min().orElse(Integer.MAX_VALUE), feature.start);

        // Calculate the new end position as the maximum of the children's end positions and the feature's current end.
        int end = Math.max(children.values().stream()
                .flatMapToInt(ls -> ls.stream().mapToInt(t -> t.b))
                .max().orElse(Integer.MIN_VALUE), feature.end);

        // Remove children with level 1 SO terms.
        children.keySet().removeIf(k -> SO.get(k) == 1);

        // Add the current feature's type and location as a child.
        children.put(feature.type, Collections.singletonList(new Tuple<>(feature.start, feature.end)));

        // Create a new feature with the updated type, location, and attributes.
        Feature updatedFeature = new Feature(feature.name, feature.contig, start, end, feature.strand, "gene", feature._uid);
        updatedFeature.setAttributes(feature.getAttributes());
        updatedFeature.setChildren(children);

        // Store the updated feature back into the features map.
        features.put(feature.name, updatedFeature);
    }

    /**
     * Imputes a feature's children (e.g. mRNA of a gene) based on the location ranges of an existing source children type.
     * I.e., if the source type exists in the children of the given feature and the target type does not,
     * this method calculates the minimum start and maximum end positions from the source type's ranges
     * and creates a new source children of the specified type.
     * <p>
     * This procedure is used to correct, for example, missing mRNA children of a gene, if a CDS child is present.
     *
     * @param children   A sorted map of child features, where the key is the source type and the value is a list of position ranges.
     * @param feature    The {@link Feature} object whose children are being processed.
     * @param sourceType The type of the source (e.g., "CDS") to derive ranges from.
     * @param targetType The type of the target (e.g., "mRNA") to impute.
     */
    private void inferFeatureChild(SortedMap<String, List<Tuple<Integer, Integer>>> children, Feature feature,
                                   String sourceType, String targetType) {
        // Check if the source type exists and the target type does not.
        if (children.containsKey(sourceType) && !children.containsKey(targetType)) {
            // Get the ranges of the source type.
            List<Tuple<Integer, Integer>> sourceRanges = children.get(sourceType);

            // Calculate the minimum start and maximum end positions from the source ranges.
            int start = sourceRanges.stream().mapToInt(t -> t.a).min().orElse(0);
            int end = sourceRanges.stream().mapToInt(t -> t.b).max().orElse(0);

            // If valid start and end positions are found, create a new range for the target type.
            if (start > 0 && end > 0) {
                children.put(targetType, Collections.singletonList(new Tuple<>(start, end)));

                // Store the updated children back into the feature.
                feature.setChildren(children);
            }
        }
    }

    /**
     * Transfers sample attributes from the `sampleInfo` map to the corresponding samples in the `samples` map.
     * <p>
     * This method iterates through the entries in the `sampleInfo` map. For each entry, it checks if a sample
     * with the corresponding name exists in the `samples` map. If the sample exists, it adds any attributes
     * from the `sampleInfo` entry that are not already present in the sample.
     */
    private void transferSampleAttributes() {
        for (Map.Entry<String, Map<String, String>> entry : this.sampleInfo.entrySet()) {
            if (hasSample(entry.getKey())) {
                samples.get(entry.getKey()).addAttributesIfAbsent(entry.getValue());
            }
        }
    }

    /**
     * Adds a sample to the storage with the specified name, if not already present.
     * <p>
     * This method creates a new {@link Sample} object with the given name and the current number of features.
     * It then adds the sample to the {@link #samples} map. If there is any sample-specific metadata available
     * in {@link #sampleInfo}, it is added to the sample as attributes.
     *
     * @param name The name of the sample to add.
     */
    private void addSample(String name) {
        if (!hasSample(name)) {
            this.samples.put(name, new Sample(name, features.size())); // Create and add a new sample to the map.
            this.samples.get(name).addAttributesIfAbsent(this.sampleInfo.getOrDefault(name, Collections.emptyMap())); // Add metadata if available.
        }
    }

    /**
     * Removes the occurrences of the specified samples from the storage.
     * <p>
     * This method iterates through all contigs and samples provided in the input collection.
     * For each sample, it removes the sample's occurrences from the contig's variants and
     * the feature's alleles. If a variant or allele no longer has any occurrences after
     * removal, it is deleted from the storage.
     *
     * @param samples A collection of {@link Sample} objects whose occurrences are to be removed.
     */
    public void removeSampleOccurrence(Collection<Sample> samples) {
        for (Contig contig : contigs.values()) {
            for (Sample sample : samples) {
                // Remove the sample's occurrences from the contig's variants.
                contig.getVariantsBySample(sample.name).forEach(variant -> {
                    VariantInformation variantInformation = contig.getVariantInformation(variant.a, variant.b);
                    variantInformation.occurrence.get(Constants.$Attributable_samplesOccurrence).remove(sample.name);
                    if (variantInformation.occurrence.get(Constants.$Attributable_samplesOccurrence).isEmpty())
                        contig.variants.get(variant.a).remove(variant.b);
                    if (contig.variants.get(variant.a).isEmpty())
                        contig.variants.remove(variant.a);
                });

                // Remove the sample's occurrences from the feature's alleles.
                for (Map.Entry<String, String> e : sample.getAlleles()) {
                    features.get(e.getKey()).alleles.get(e.getValue()).occurrence.remove(sample.name);
                    if (features.get(e.getKey()).alleles.get(e.getValue()).occurrence.isEmpty())
                        features.get(e.getKey()).alleles.remove(e.getValue());
                }
            }
        }
    }

    /**
     * Transfers sample information from variant records to the storage.
     * <p>
     * This method processes variant records stored in {@link VcfHandler#records} and updates the storage
     * with variant calls for each sample, contig, and position. It calculates the total depth of coverage (DP),
     * determines the best allele based on phred-scaled likelihoods (PL) or allele depth (AD), and builds
     * a variant call string. The method also handles exclusions for low frequency, low coverage, and specific
     * variants, and skips passing reference calls.
     */
    public void transferSampleInformation() {
        // Iterate over each sample in the variant records.
        vcfHandler.records.forEach((sampleName, contigs) -> {
            // Iterate over each contig for the current sample.
            contigs.forEach((contigName, variants) -> {
                // Iterate over each variant context for the current contig.
                variants.forEach((POS, alleles) -> {
                    // Calculate the total depth of coverage (DP) for the alleles.
                    int contextDP = alleles.stream().mapToInt(VcfHandler.Allele::AD).sum();
                    String callPrefix = Constants.EMPTY;
                    StringBuilder callBuilder = new StringBuilder();
                    int callIdx, callGQ = -1;
                    float callFrequency;

                    // Determine the best allele based on PL (phred-scaled likelihoods) or AD (allele depth).
                    if (alleles.stream().noneMatch(allele -> allele.PL() == -1)) {
                        // Use PL values to determine the best allele.
                        int minPL = alleles.stream().mapToInt(VcfHandler.Allele::PL).min().orElse(Integer.MAX_VALUE);
                        callIdx = alleles.stream().mapToInt(VcfHandler.Allele::PL).boxed().toList().indexOf(minPL);
                        callGQ = Math.min(
                                alleles.stream().mapToInt(VcfHandler.Allele::PL).sorted().skip(1).findFirst().orElse(0) - minPL,
                                99);
                        callFrequency = alleles.get(callIdx).AD() / (float) contextDP;
                    } else {
                        // Use AD values to determine the best allele.
                        int maxAD = alleles.stream().mapToInt(VcfHandler.Allele::AD).max().orElse(-1);
                        callIdx = alleles.stream().mapToInt(VcfHandler.Allele::AD).boxed().toList().indexOf(maxAD);
                        callFrequency = maxAD / (float) contextDP;
                    }

                    // Skip excluded variants.
                    if (callIdx > 0 && isVariantExcluded(
                            contigName, POS,
                            SequenceOperations.stripGaps(alleles.get(callIdx).REF()),
                            SequenceOperations.stripGaps(alleles.get(callIdx).ALT())
                    )) return;

                    // Set call prefix for low frequency or coverage.
                    if (callFrequency < minimumFrequency()) callPrefix = Constants.lowFrequencyCallPrefix;
                    if (contextDP < minimumCoverage()) callPrefix = Constants.lowCoverageCallPrefix;

                    // Skip passing reference calls.
                    if (callIdx == 0 && callPrefix.equals(Constants.EMPTY)) return;

                    // Build the call string with allele information.
                    callBuilder.append(callPrefix).append(callIdx).append(Constants.SEMICOLON)
                            .append(contextDP).append(Constants.SEMICOLON).append(callGQ).append(Constants.SEMICOLON);
                    alleles.forEach(allele -> callBuilder.append(allele.REF()).append(Constants.COLON)
                            .append(allele.ALT()).append(Constants.COLON).append(allele.AD())
                            .append(Constants.COLON).append(allele.PL()).append(Constants.COMMA));
                    callBuilder.deleteCharAt(callBuilder.length() - 1);

                    // Add the variant call to the sample.
                    addVariantCallToSample(sampleName, contigName, POS, callBuilder.toString());
                });
            });
        });
    }

    /**
     * Transfers variant information from sample variant calls to the storage.
     * <p>
     * This method processes variant calls for each sample and contig, resolving conflicts and handling
     * deletions, insertions, and mixed InDels. It ensures that variants are stored in a canonical format
     * and accounts for the effects of upstream deletions on downstream variants. Variants are added to
     * the contig's variant map, and warnings are logged for conflicts or unhandled cases.
     */
    public void transferVariantsInformation() {
        for (Sample sample : getSamplesToUpdate()) {
            for (String contigName : sample.variantCalls.keySet()) {
                // Establish a sorted list of canonical variants for the sample and contig.
                TreeMap<Integer, Tuple<String, String>> variants = new TreeMap<>();

                // Add variants to the list, resolving conflicts by keeping the more specific reference content.
                BiConsumer<Integer, Tuple<String, String>> addVariant = (position, content) -> {
                    Tuple<String, String> previousContent = variants.get(position);
                    if (previousContent == null) {
                        variants.put(position, content);
                    } else if (!previousContent.equals(content) && !previousContent.a.contains(content.a)) {
                        if (content.a.contains(previousContent.a)) {
                            variants.put(position, content);
                        } else {
                            Logging.logWarning("Conflict of variant records at position %d on contig %s for sample %s. %s (stored) and %s."
                                    .formatted(position, contigName, sample.name, previousContent, content));
                        }
                    }
                };

                // Process each variant call for the contig.
                for (Map.Entry<Integer, String> variantCall : sample.variantCalls.get(contigName).entrySet()) {
                    int position = variantCall.getKey();
                    String[] context = variantCall.getValue().split(Constants.SEMICOLON);
                    if (context[0].startsWith("0")) continue; // Skip reference calls.

                    boolean isAmbiguous = (context[0].startsWith(Constants.lowFrequencyCallPrefix)
                            || context[0].startsWith(Constants.lowCoverageCallPrefix));
                    if (!storeFiltered() && isAmbiguous) continue; // Skip ambiguous calls if not allowed.

                    int callIndex = Integer.parseInt(isAmbiguous ? context[0].substring(1) : context[0]);
                    String[] genotype = context[3].split(Constants.COMMA)[callIndex].split(Constants.COLON);
                    String referenceContent = genotype[0];
                    String alternativeContent = isAmbiguous ? (Constants.anyNucleotide.repeat(referenceContent.length())) : genotype[1];

                    if (VariantInformation.isPaddedCanonicalVariant(referenceContent, alternativeContent)) {
                        addVariant.accept(position, new Tuple<>(referenceContent, alternativeContent));
                    } else {
                        for (Triple<Integer, String, String> canonicalVariant : SequenceOperations.getCanonicalVariants(referenceContent, alternativeContent)) {
                            addVariant.accept(position + canonicalVariant.getLeft(), new Tuple<>(canonicalVariant.getMiddle(), canonicalVariant.getRight()));
                        }
                    }
                }

                // Process variants to account for deletions and mixed InDels.
                StringBuilder referenceBuilder = new StringBuilder();
                StringBuilder alternativeBuilder = new StringBuilder();
                int variantStartPosition = 0;
                int deletionExtension = 0;

                // Helper function to resolve and add a variant to the contig.
                Consumer<Integer> resolveVariant = (position) -> {
                    if (!VariantInformation.isPaddedCanonicalVariant(referenceBuilder.toString(), alternativeBuilder.toString())) {
                        Tuple<String, String> realignedMixedIndel =
                                SequenceOperations.globalNucleotideSequenceAlignment(
                                        SequenceOperations.stripGaps(referenceBuilder.toString()),
                                        SequenceOperations.stripGaps(alternativeBuilder.toString()),
                                        Math.max(referenceBuilder.length(), alternativeBuilder.length()),
                                        2,
                                        SequenceOperations.MarginalGaps.FORBID,
                                        SequenceOperations.MarginalGaps.PENALIZE,
                                        null
                                );
                        ArrayList<Triple<Integer, String, String>> resolvedVariants = SequenceOperations.getCanonicalVariants(realignedMixedIndel.a, realignedMixedIndel.b);
                        for (Triple<Integer, String, String> resolvedVariant : resolvedVariants) {
                            addVariantToContig(contigName, sample.name, position + resolvedVariant.getLeft(), resolvedVariant.getMiddle(), resolvedVariant.getRight());
                        }
                    } else {
                        addVariantToContig(contigName, sample.name, position, referenceBuilder.toString(), alternativeBuilder.toString());
                    }
                };

                // Iterate through the sorted variants and handle deletions and insertions.
                for (Map.Entry<Integer, Tuple<String, String>> variant : variants.entrySet()) {
                    int position = variant.getKey();
                    Tuple<String, String> variantContent = variant.getValue();
                    String referenceContent = variantContent.a;
                    String alternativeContent = variantContent.b;

                    if (position > deletionExtension && referenceBuilder.length() > 0 && alternativeBuilder.length() > 0) {
                        resolveVariant.accept(variantStartPosition);
                        referenceBuilder.setLength(0);
                        alternativeBuilder.setLength(0);
                        deletionExtension = 0;
                    }

                    if (deletionExtension == 0 && referenceBuilder.length() == 0 && alternativeBuilder.length() == 0) {
                        if (alternativeContent.startsWith("*")) {
                            alternativeContent = alternativeContent.replaceFirst("\\*", String.valueOf(referenceContent.charAt(0)));
                            if (Objects.equals(referenceContent, alternativeContent)) continue;
                        }
                        if (VariantInformation.isDeletion(referenceContent, alternativeContent, true)) {
                            referenceBuilder.append(referenceContent);
                            alternativeBuilder.append(alternativeContent);
                            variantStartPosition = position;
                            deletionExtension = position + alternativeContent.length() - 1;
                        } else {
                            referenceBuilder.append(referenceContent);
                            alternativeBuilder.append(alternativeContent);
                            resolveVariant.accept(position);
                            referenceBuilder.setLength(0);
                            alternativeBuilder.setLength(0);
                        }
                        continue;
                    }

                    if (position <= deletionExtension) {
                        if (VariantInformation.isSubstitution(referenceContent, alternativeContent)) {
                            continue;
                        }
                        if (VariantInformation.isDeletion(referenceContent, alternativeContent, true)) {
                            int updatedDeletionExtension = position + alternativeContent.length() - 1;
                            if (updatedDeletionExtension > deletionExtension) {
                                referenceBuilder.append(StringUtils.right(referenceContent, updatedDeletionExtension - deletionExtension));
                                alternativeBuilder.append(StringUtils.right(alternativeContent, updatedDeletionExtension - deletionExtension));
                                deletionExtension = updatedDeletionExtension;
                            }
                            continue;
                        }
                        if (VariantInformation.isInsertion(referenceContent, alternativeContent, true)) {
                            int offset = position - variantStartPosition;
                            alternativeBuilder.replace(offset, offset + 1, alternativeBuilder.charAt(offset) + alternativeContent.substring(1));
                            referenceBuilder.replace(offset, offset + 1, referenceBuilder.charAt(offset) + referenceContent.substring(1));
                            continue;
                        }
                    }

                    Logging.logWarning("Failed to handle variant %s at position %d on contig %s for sample %s."
                            .formatted(referenceContent + ">" + alternativeContent, position, contigName, sample.name));
                }

                if (referenceBuilder.length() > 0 && alternativeBuilder.length() > 0) {
                    resolveVariant.accept(variantStartPosition);
                }
            }
        }
    }

    /**
     * Adds a variant call to the sample stored in {@link #samples} with the key {@code sampleName}. The call is expected to be in the
     * format specified by {@link Sample#variantCallPattern}, where {@code CALL_INDEX} is one of {@code f} (low frequency), {@code x} (low coverage),
     * or the index of the alternative call (starting at 0 for the reference allele). {@code DP, AD} as defined in VCF specification
     * (<a href="https://samtools.github.io/hts-specs/VCFv4.2.pdf">samtools.github.io/hts-specs/VCFv4.2.pdf</a>).
     * <p>
     * If a call is already present at the specified position, it will be overwritten!
     *
     * @param sampleName The name of the sample.
     * @param contigName The name of the contig.
     * @param position   The position of the variant call on the contig.
     * @param call       The variant call to add.
     */
    private void addVariantCallToSample(String sampleName, String contigName, int position, String call) {
        if (!Sample.variantCallPattern.matcher(call).matches())
            throw new IllegalArgumentException("Invalid call format %s for sample %s at position %d on contig %s. Expected: %s."
                    .formatted(call, sampleName, position, contigName, Sample.variantCallPattern.pattern()));
        if (this.hasSample(sampleName)) {
            Sample sample = this.getSample(sampleName);
            if (!sample.variantCalls.containsKey(contigName)) sample.variantCalls.put(contigName, new TreeMap<>());
            sample.variantCalls.get(contigName).put(position, call);
        } else {
            throw new IllegalArgumentException(String.format("Failed to add variant call; Sample %s not found.", sampleName));
        }
    }

    /**
     * Adds a variant to the specified contig in the storage.
     * <p>
     * This method validates the variant to ensure it is in a canonical padded format. It then retrieves all features
     * annotated for the given position and checks if the contig, sample, and features exist in the storage. If valid,
     * the variant is added to the contig's variant map, and its occurrences in the sample and features are updated.
     * If the variant is novel, it is added to the `novelVariants` list.
     *
     * @param contigName         The name of the contig to which the variant belongs.
     * @param sampleName         The name of the sample associated with the variant.
     * @param position           The position of the variant on the contig.
     * @param referenceContent   The reference base(s) for the variant.
     * @param alternativeContent The alternative base(s) for the variant.
     * @throws IllegalArgumentException If the variant is not in canonical format, or if the contig, sample, or features are not found.
     */
    private void addVariantToContig(String contigName, String sampleName, int position, String referenceContent, String alternativeContent) {
        if (!VariantInformation.isPaddedCanonicalVariant(referenceContent, alternativeContent)) {
            throw new IllegalArgumentException("Failed to add non-canonical variant %s > %s at position %d on contig %s."
                    .formatted(referenceContent, alternativeContent, position, contigName));
        }
        // Get all features that are annotated for the position.
        Set<String> featureNames = features.values().stream()
                .filter(feature -> feature.start <= position && feature.end >= position)
                .map(feature -> feature.name).collect(Collectors.toSet());
        if (this.hasContig(contigName) && this.hasSample(sampleName) && featureNames.stream().allMatch(this::hasFeature)) {
            VariantInformation variantInformation;
            Contig contig = this.getContig(contigName);
            contig.variants.putIfAbsent(position, new HashMap<>());
            variantInformation = contig.variants.get(position)
                    .putIfAbsent(alternativeContent, new VariantInformation(referenceContent, alternativeContent));
            contig.variants.get(position).get(alternativeContent).addSampleOccurrence(sampleName);
            featureNames.forEach(featureName -> contig.variants.get(position).get(alternativeContent).addFeatureOccurrence(featureName));
            if (Objects.nonNull(variantInformation) && !variantInformation.reference.equals(referenceContent)) {
                Logging.logWarning("Variant call %s at position %d on contig %s occurred with different reference bases %s and %s."
                        .formatted(alternativeContent, position, contigName, variantInformation.reference, referenceContent));
            }
            if (Objects.isNull(variantInformation))
                novelVariants.add(new ImmutableTriple<>(contigName, position, alternativeContent));
        } else {
            throw new IllegalArgumentException(String.format("Failed to add variant to contig %s as either the contig, sample (%s), or feature (%s) were not found.",
                    contigName, sampleName, String.join(", ", featureNames)));
        }
    }

    /**
     * Encapsulates utility methods for analyzing and processing variant call format (VCF) files.
     * <p>
     * This class implements methods and data structures for handling VCF file analysis,
     * including allele information storage, variant context processing, and integration with
     * external tools like SnpEff for annotation.
     */
    private class VcfHandler {

        /**
         * Used to store information about a specific allele,
         * including its reference sequence (REF), alternate sequence (ALT),
         * depth of coverage (DP), and phred-scaled likelihood (PL).
         * <p>
         * If the reference and alternate sequences represent a canonical variant,
         * the constructor ensures that the sequences are padded with gaps to align
         * their lengths.
         *
         * @param REF The reference allele sequence.
         * @param ALT The alternate allele sequence.
         * @param AD  The depth of coverage for the allele.
         * @param PL  The phred-scaled likelihood for the allele.
         */
        public record Allele(String REF, String ALT, int AD, int PL) {
            /**
             * Constructs an Allele instance and ensures that the reference and
             * alternate sequences are padded with gaps if they represent a
             * canonical variant.
             */
            public Allele {
                if (VariantInformation.isCanonicalVariant(REF, ALT)) {
                    REF = SequenceOperations.padGaps(REF, ALT.length());
                    ALT = SequenceOperations.padGaps(ALT, REF.length());
                } else if (!Objects.equals(ALT, Constants.DOT)) {
                    // If no hom. reference call or canonical variant is represented, a re-alignment is performed.
                    Tuple<String, String> alignment = SequenceOperations.globalNucleotideSequenceAlignment(REF, ALT, 3, 2,
                            SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, null);
                    REF = alignment.a;
                    ALT = alignment.b;
                }
            }
        }

        /**
         * A nested HashMap structure used to store allele information for VCF records.
         * <p>
         * The structure is organized as follows:
         * <ul>
         *   <li>The outermost key is a {@link String} representing the sample name.</li>
         *   <li>The second-level key is a {@link String} representing the contig name.</li>
         *   <li>The third-level key is an {@link Integer} representing the start position of the variant.</li>
         *   <li>The innermost value is an {@link ArrayList} of {@link Allele} objects, representing the alleles at the given position.</li>
         * </ul>
         * This structure is used to store and organize allele data for efficient access and processing.
         */
        private final HashMap<String, HashMap<String, HashMap<Integer, ArrayList<Allele>>>> records = new HashMap<>();

        /**
         * Counter for processed genotype records.
         * <p>
         * This field is used to keep track of the number of genotype records that have been processed
         * during the analysis of VCF files. It is initialized to 0 and can be incremented as records
         * are processed.
         */
        private long processedGenotypes = 0; // Counter for processed records.

        /**
         * Clears the nested HashMap structure storing allele information for variant records.
         * <p>
         * This method removes all entries from the {@link #records} map, effectively resetting it.
         * It is typically used to free up memory or prepare for a new analysis.
         */
        private void clearRecords() {
            records.clear();
        }

        /**
         * Analyzes a VCF file and processes its records for the storage.
         * <p>
         * This method reads a VCF file, queries features from the storage, and processes
         * variant records within the specified regions. It filters out excluded positions
         * and processes each genotype in the record. Genotypes that are no-calls or lack
         * required attributes (AD and DP for alternatives, DP for references) are skipped.
         * Valid genotypes are further processed using {@link #transferGenotype}.
         *
         * @param vcfFile The VCF file to be analyzed.
         * @throws IOException If an error occurs while reading the VCF file.
         */
        private void analyzeVcf(File vcfFile) throws IOException {
            try (VCFFileReader vcfFileReader = IO.initializeVCFFileReader(vcfFile)) {
                // Iterate over features in the storage.
                for (Feature feature : Storage.this.getFeatures()) {
                    vcfFileReader
                            .query(feature.contig, feature.start, feature.end) // Query VCF for the feature region.
                            .stream()
                            .filter(record -> !Storage.this.isPositionExcluded(record.getContig(), record.getStart())) // Exclude specific positions.
                            .forEach(record -> {
                                // Process each genotype in the variant context.
                                for (Genotype genotype : record.getGenotypes()) {
                                    // Skip genotypes marked as no-call.
                                    if (genotype.isNoCall()) {
                                        continue;
                                    }
                                    // Log an error and skip if alternative calls lack AD or DP attributes.
                                    if (genotype.hasAltAllele() && !(genotype.hasAD() && genotype.hasDP())) {
                                        Logging.logError(String.format("VCF record analysis failed; AD and DP genotype attributes are required to process alternative calls."
                                                        + " At least one is missing in file %s at position %s for genotype %s and the variant will be skipped.",
                                                vcfFile.getAbsolutePath(), record.getStart(), genotype.getSampleName()));
                                        continue;
                                    }
                                    // Log an error and skip if reference calls lack the DP attribute.
                                    if (genotype.isHomRef() && !genotype.hasDP()) {
                                        Logging.logError(String.format("VCF record analysis failed; DP genotype attribute is required to process reference calls."
                                                        + " The attribute is missing in file %s at position %s for genotype %s and the variant will be skipped.",
                                                vcfFile.getAbsolutePath(), record.getStart(), genotype.getSampleName()));
                                        continue;
                                    }
                                    // Transfer valid genotype information to the storage.
                                    transferGenotype(record, genotype);
                                }
                            });
                }
            }
        }

        /**
         * Computes the maximum end position for each contig across multiple VCF files.
         * <p>
         * This method iterates through a list of VCF files, reads their intervals, and updates a map
         * where the key is the contig name and the value is the maximum end position observed for that
         * contig. If a contig already exists in the map, the method updates its value to the maximum
         * of the current value and the new interval's end position.
         *
         * @return A {@link HashMap} where the keys are contig names and the values are the maximum end positions.
         * @throws IOException If an error occurs while reading the VCF files.
         */
        private HashMap<String, Integer> inferContigs() throws IOException {
            // Initialize a map to store the maximum end position for each contig.
            HashMap<String, Integer> contigs = new HashMap<>();

            // Iterate over the list of VCF files.
            for (File vcfFile : vcfFiles) {
                // Open a VCF file reader for the current file.
                try (VCFFileReader reader = IO.initializeVCFFileReader(vcfFile)) {
                    // Iterate over the intervals in the VCF file and update the map.
                    reader.toIntervalList().getIntervals().forEach(interval ->
                            contigs.merge(interval.getContig(), interval.getEnd(), Math::max) // Merge the current end position with the maximum.
                    );
                }
            }

            // Return the map containing the maximum end positions for each contig.
            return contigs;
        }

        /**
         * Transfers genotype information from a VCF record to the storage.
         * <p>
         * This method processes a genotype from a VCF record, ensuring that the associated contig and sample
         * exist in the storage. It extracts reference and alternative alleles, computes phred-scaled likelihoods (PL),
         * and updates the storage with allele information for the given sample, contig, and position.
         *
         * @param context  The {@link VariantContext} object containing the variant information.
         * @param genotype The {@link Genotype} object representing the sample's genotype for the variant.
         */
        private void transferGenotype(VariantContext context, Genotype genotype) {
            // Extract sample name and ensure contig and sample exist in storage.
            // TODO: Sample names in the VCF can have a "$" suffix to be merged within one sample in musial.
            String sampleName = genotype.getSampleName().split("\\$")[0];
            Storage.this.addSample(sampleName);

            // Process reference and alternative alleles.
            String REF = context.getReference().getBaseString();
            int[] ADs = genotype.isHomRef() ? new int[]{genotype.getDP()} : genotype.getAD(); // Allele depths.
            int[] PLs = computePLs(genotype, ADs.length); // Phred-scaled likelihoods.

            // Iterate through alleles and add or update their context in the storage.
            for (int i = 0; i < ADs.length; i++) {
                String ALT = (i == 0) ? Constants.DOT : context.getAlleles().get(i).getBaseString(); // Alternate allele.
                updateRecords(sampleName, context, REF, ALT, ADs[i], PLs[i]); // Update storage with allele information.
            }

            // Count processed genotype records.
            processedGenotypes++;
        }

        /**
         * Adds or updates allele information in the nested {@link #records} structure for a given sample, contig, and start position.
         * <p>
         * This method ensures that the nested structure for storing allele information exists for the specified
         * sample, contig, and start position. It then stores the provided reference and alternate sequences,
         * depth of coverage, and phred-scaled likelihood as an {@link Allele} - if an allele with
         * the same reference and alternate sequences already exists, its depth of coverage and likelihood are
         * updated. Otherwise, the new allele is added to the list.
         *
         * @param sampleName The name of the sample.
         * @param context    The {@link VariantContext} containing contig and start position information.
         * @param REF        The reference allele sequence.
         * @param ALT        The alternate allele sequence.
         * @param AD         The depth of coverage for the allele.
         * @param PL         The phred-scaled likelihood for the allele.
         */
        private void updateRecords(String sampleName, VariantContext context, String REF, String ALT, int AD, int PL) {
            // Ensure the nested structure exists for the given sample, contig, and start position.
            ArrayList<Allele> alleles = records
                    .computeIfAbsent(sampleName, k -> new HashMap<>())
                    .computeIfAbsent(context.getContig(), k -> new HashMap<>())
                    .computeIfAbsent(context.getStart(), k -> new ArrayList<>());

            // Create a new Allele object.
            Allele allele = new Allele(REF, ALT, AD, PL);

            // Check if the allele already exists and update it, otherwise add it.
            int idx = alleles.indexOf(allele);
            if (idx > -1) {
                Allele existingAllele = alleles.get(idx);
                alleles.set(idx, new Allele(
                        REF, ALT,
                        AD + existingAllele.AD,
                        Math.min(PL, existingAllele.PL)
                ));
            } else {
                alleles.add(allele);
            }
        }

        /**
         * Computes the phred-scaled likelihoods (PLs) for a given genotype.
         * <p>
         * This method calculates the PL values for each allele in the genotype based on its ploidy and
         * either the PL or GL (genotype likelihood) attributes. If neither PL nor GL is present, the
         * method returns an array filled with -1.
         *
         * @param genotype The {@link Genotype} object containing allele and likelihood information.
         * @param n        The number of alleles in the VCF record the genotype is associated with.
         * @return An array of integers representing the phred-scaled likelihoods for each allele.
         */
        private int[] computePLs(Genotype genotype, int n) {
            // Get the number of alleles in the genotype.
            int[] PL = new int[n];

            // If neither PL nor GL attributes are present, return an array filled with -1.
            if (!genotype.hasPL() && !genotype.hasAnyAttribute("GL")) {
                Arrays.fill(PL, -1);
                return PL;
            }

            // Define a function to compute the PL index based on the genotype's ploidy.
            Function<Integer, Integer> getPLIndex = switch (genotype.getPloidy()) {
                case 1 -> i -> i; // Haploid case: index is the allele index.
                case 2 -> i -> (i * (i + 1)) / 2 + i; // Diploid case: triangular index calculation.
                default -> i -> { // Higher ploidy case: calculate using binomial coefficients.
                    int index = 0;
                    for (int m = 1; m <= genotype.getPloidy(); m++) {
                        index += BinomialCoefficientUtil.binomialCoefficient(i + m + 1, m);
                    }
                    return index;
                };
            };

            // If the genotype has PL values, retrieve them using the computed index.
            if (genotype.hasPL()) {
                for (int i = 0; i < n; i++) {
                    PL[i] = genotype.getPL()[getPLIndex.apply(i)];
                }
            } else {
                // If the genotype has GL values, convert them to PL values.
                List<Float> GL = ((List<?>) genotype.getAnyAttribute("GL")).stream()
                        .map(e -> (Float) e)
                        .toList();
                for (int i = 0; i < n; i++) {
                    PL[i] = -10 * Math.round(GL.get(getPLIndex.apply(i)));
                }
            }

            return PL;
        }

        /**
         * Runs a SnpEff analysis on the outer {@link Storage} and stores the results in the storage.
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
         *   <li>Handles errors during the SnpEff build and annotation processes, logging them and saving error logs to the output directory.</li>
         *   <li>Cleans up the temporary directory after the analysis is complete.</li>
         * </ul>
         *
         * @throws MusialException If the SnpEff annotation process fails.
         * @throws IOException     If an error occurs while reading or writing files.
         */
        private void runSnpEffAnalysis() throws MusialException, IOException {
            // Generate a temporary directory for snpEff.
            String prefix = "%s-%s".formatted("snpEff", RandomStringUtils.randomAlphanumeric(6));
            Path temp = Files.createTempDirectory(prefix);
            try {
                // Create map to store variant pointers and write storage variants to temporary VCF file.
                ArrayList<Tuple<Triple<String, Integer, String>, VariantInformation>> variants = new ArrayList<>(Storage.this.novelVariants.size());
                for (Triple<String, Integer, String> variant : Storage.this.novelVariants) {
                    variants.add(new Tuple<>(variant, Storage.this.getContig(variant.getLeft()).getVariantInformation(variant.getMiddle(), variant.getRight())));
                }
                IO.writeFile(Path.of(temp + "/variants" + FileExtensions.VCF), IO.generateVcfContent(variants));

                // Write reference .gff and .fasta to temp. target directory.
                IO.writeFile(Path.of(temp + "/data/reference/genes.gff"), IO.generateGffContent(Storage.this));
                IO.writeFile(Path.of(temp + "/data/genomes/reference.fa"), IO.generateReferenceFastaContent(Storage.this));

                // Copy snpEff config and JAR to temp. target directory.
                Path snpEffConfigPath = Path.of(temp + "/snpEff.config");
                IO.copyResourceToFile("/snpEff/snpEff.config", snpEffConfigPath);
                IO.copyResourceToFile("/snpEff/snpEff.jar", Path.of(temp + "/snpEff.jar"));

                // Add reference .fasta and .gff information to snpEff.config.
                String codonTableConfig = Storage.this.getContigs()
                        .stream()
                        .map(contig -> "reference.genome.%s : Bacterial_and_Plant_Plastid".formatted(contig.name))
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
                try {
                    List<String> lines = IO.readFile(new File(temp + "/annotation" + FileExtensions.VCF))
                            .stream()
                            .filter(s -> !s.startsWith(Constants.SIGN))
                            .toList();
                    int index = 0;
                    for (String line : lines) {
                        String[] annotationFields = line.split("\t");
                        if (!annotationFields[7].equals(".")) {
                            annotationFields = annotationFields[7].replace("ANN=", "").split(Constants.COMMA)[0].split("\\|");
                            for (int i = 0; i < annotationFields.length; i++) {
                                if (i == 1 || i == 2 || i == 5 || i == 7 || i == 12 || i == 13) {
                                    variants.get(index).b.addAttributeIfAbsent(
                                            Constants.snpEffAttributeKeyPrefix + Constants.snpEffKeys.get(i),
                                            i == 1 ? annotationFields[i].replaceAll("&", Constants.COMMA) : annotationFields[i]
                                    );
                                } else if (i == 6) {
                                    variants.get(index).b.addAttributeIfAbsent(
                                            Constants.snpEffAttributeKeyPrefix + Constants.snpEffKeys.get(i),
                                            annotationFields[i].split("-")[1]
                                    );
                                }
                            }
                        }
                        index++;
                    }
                } catch (FileNotFoundException e) {
                    throw new MusialException(String.format("Failed to read SnpEff annotation. %s", e.getMessage()));
                }
            } finally {
                File buildErrorFile = new File(temp + "/snpEff.build.err");
                if (buildErrorFile.exists() && buildErrorFile.length() != 0) {
                    Logging.logError("SnpEff `build` has raised an error or warning; a copy of the log file is in the output directory - the annotations may be incorrect.");
                    FileUtils.copyFile(buildErrorFile, new File(Musial.outputDirectory.getAbsolutePath()
                            + "/musial_snpeff_build_%s.error".formatted(Musial.runId)));
                }
                File annErrorFile = new File(temp + "/snpEff.ann.err");
                if (annErrorFile.exists() && annErrorFile.length() != 0) {
                    Logging.logError("SnpEff `ann` has raised an error or warning; a copy of the log file is in the output directory - the annotations may be incorrect.");
                    FileUtils.copyFile(annErrorFile, new File(Musial.outputDirectory.getAbsolutePath()
                            + "/musial_snpeff_ann_%s.error".formatted(Musial.runId)));
                }
                // Clean up temporary directory.
                FileUtils.deleteDirectory(temp.toFile());
            }
        }

    }
}