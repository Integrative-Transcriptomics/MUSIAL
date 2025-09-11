package cli;

import com.fasterxml.jackson.databind.JsonNode;
import com.github.fge.jackson.JsonLoader;
import com.github.fge.jsonschema.core.exceptions.ProcessingException;
import com.github.fge.jsonschema.core.report.ProcessingReport;
import com.github.fge.jsonschema.main.JsonSchema;
import com.github.fge.jsonschema.main.JsonSchemaFactory;
import com.google.gson.Gson;
import exceptions.MusialException;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import main.Musial;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.io.file.PathUtils;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import util.IO;
import util.Logging;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Handles the {@code build} task CLI parameters.
 * <p>
 * This class defines the command-line options and validation logic for the {@code build} task. It allows users to specify a JSON file
 * containing the task parameters for MUSIAL. Specifically, it is a direct reflection of the build configuration.
 */
public class CLIBuild implements CLI {

    /**
     * Creates and configures the command-line options for the Build task.
     * <p>
     * This method defines the available command-line options for the Build task, including their descriptions, argument requirements, and
     * whether they are mandatory.
     *
     * @return An {@link Options} object containing the defined command-line options.
     */
    public static Options options() {
        Options options = new Options();
        options.addOption(Option.builder("C")
                .longOpt("configuration")
                .desc("Path to a JSON file specifying the build task parameter configuration for MUSIAL. Visit the documentation for " +
                        "details.")
                .hasArg()
                .required()
                .build());
        return options;
    }

    /**
     * The minimal depth of coverage wrt. reads of a variant to be accepted.
     */
    public final int minimalCoverage;

    /**
     * The minimal frequency wrt. reads supporting a nucleotide variant call for a variant to be accepted.
     */
    public final double minimalFrequency;

    /**
     * If filtered variants are stored as ambiguous nucleotides in the storage.
     */
    public final boolean storeFiltered;

    /**
     * If annotation of variants is skipped during the build process.
     */
    public final boolean skipAnnotation;

    /**
     * If typing of samples is skipped during the build process.
     */
    public final boolean skipTyping;

    /**
     * A map containing the masked positions, with the contig ID as the key and a set of excluded positions as the value.
     */
    public final Map<String, Set<Integer>> maskedPositions;

    /**
     * The reference sequence file used for variant annotation and typing.
     */
    public final IndexedFastaSequenceFile reference;

    /**
     * The list of features parsed from the annotation file.
     */
    public final FeatureList featureList;

    /**
     * The output path where the MUSIAL storage will be saved.
     */
    public final Path output;

    /**
     * The list of VCF files to be processed.
     */
    public final List<Path> vcfFiles;

    /**
     * A map containing the VCF metadata, with the sample ID as the key and another map of metadata attributes and their values as the
     */
    public final Map<String, Map<String, String>> vcfMeta;

    /**
     * A map containing the features to be matched, with the feature ID as the key and another map of at least a `key` and `value` attribute
     * to match the feature with.
     */
    public final Map<String, Map<String, String>> features;

    /**
     * Constructs a CLIBuild instance by parsing and initializing configuration parameters.
     * <p>
     * This constructor processes the command-line arguments to extract and validate the build task configuration.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @throws IOException         If an I/O error occurs during file reading.
     * @throws MusialException     If the configuration validation fails or required parameters are missing.
     * @throws ProcessingException If the JSON schema validation fails.
     */
    public CLIBuild(CommandLine arguments) throws IOException, MusialException, ProcessingException {
        Map<String, Object> configuration = parseConfiguration(arguments);
        this.minimalCoverage = parseMinimalCoverage(configuration);
        this.minimalFrequency = parseMinimalFrequency(configuration);
        this.storeFiltered = parseStoreFiltered(configuration);
        this.skipAnnotation = parseSkipAnnotation(configuration);
        this.skipTyping = parseSkipTyping(configuration);
        this.maskedPositions = parseMaskedPositions(configuration);
        this.reference = parseReference(configuration);
        this.featureList = parseAnnotation(configuration);
        this.output = Common.parseOutput(configuration);
        this.vcfFiles = Common.parseVcfFiles(configuration);
        this.vcfMeta = Common.parseVcfMeta(configuration);
        this.features = parseFeatures(configuration);
    }

    /**
     * Parses the build configuration from the specified command-line arguments.
     * <p>
     * This method validates the build configuration file against a predefined JSON schema to ensure correctness. If the validation is
     * successful, the configuration is parsed into a JSON object and returned as a map.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Loads the JSON schema from the classpath resource.</li>
     *     <li>Validates the configuration file against the schema.</li>
     *     <li>Parses the configuration file into a map using Gson.</li>
     * </ul>
     * <p>
     * If the validation fails or the file cannot be parsed, appropriate exceptions are thrown.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments. The configuration file path is expected
     *                  to be provided as the value of the "C" option.
     * @return A {@link Map} containing the parsed configuration parameters.
     * @throws IOException         If an I/O error occurs while reading the schema or configuration file.
     * @throws MusialException     If the configuration file is invalid or deserialization fails.
     * @throws ProcessingException If the JSON schema validation fails.
     */
    private Map<String, Object> parseConfiguration(CommandLine arguments) throws IOException, MusialException, ProcessingException {
        // Validate the build configuration against schema.
        try (InputStream schemaStream = Objects.requireNonNull(Musial.class.getResourceAsStream("/buildConfigurationSchema.json"))) {
            JsonSchema schema = JsonSchemaFactory.byDefault()
                    .getJsonSchema(JsonLoader.fromReader(new InputStreamReader(schemaStream, StandardCharsets.UTF_8)));
            JsonNode config = JsonLoader.fromPath(arguments.getOptionValue("C"));
            ProcessingReport report = schema.validate(config);
            if (!report.isSuccess()) {
                throw new MusialException("Invalid build configuration file:\n%s".formatted(report.toString()));
            }
        }

        // Parse the build configuration into a JSON object.
        try (Reader reader = Files.newBufferedReader(Path.of(arguments.getOptionValue("C")), StandardCharsets.UTF_8)) {
            //noinspection unchecked
            return new Gson().fromJson(reader, HashMap.class);
        } catch (Exception e) {
            throw new MusialException("Deserialization of build configuration file failed; %s".formatted(e.getMessage()));
        }
    }

    /**
     * Parses the minimal coverage value from the configuration map.
     * <p>
     * This method retrieves the "minimalCoverage" value from the provided configuration map. If the value is present, it is parsed as a
     * double and converted to an integer. The value must be a positive number; otherwise, a {@link NumberFormatException} is thrown. If the
     * value is not specified, the default value of 3 is used.
     * <p>
     * The method logs the parsed or default value for "minimalCoverage".
     *
     * @param configuration A {@link Map} containing the configuration parameters.
     * @return The parsed minimal coverage value as an integer.
     * @throws NumberFormatException If the "minimalCoverage" value is not a positive number.
     */
    private int parseMinimalCoverage(Map<String, Object> configuration) {
        int minimalCoverage = 3; // Default value for minimal coverage
        if (configuration.containsKey("minimalCoverage")) {
            double value = (double) configuration.get("minimalCoverage");
            if (value > 0) {
                minimalCoverage = (int) value;
                Logging.logConfig("`minimalCoverage` set to %d.".formatted(minimalCoverage));
            } else {
                throw new NumberFormatException("Value for `minimalCoverage` must be a positive integer, not %s.".formatted(value));
            }
        } else {
            Logging.logConfig("No value for `minimalCoverage` specified; defaulting to 3.");
        }
        return minimalCoverage;
    }

    /**
     * Parses the minimal frequency value from the configuration map.
     * <p>
     * This method retrieves the "minimalFrequency" value from the provided configuration map. If the value is present, it is parsed as a
     * double. The value must be within the range (0, 1]; otherwise, a {@link NumberFormatException} is thrown. If the value is not
     * specified, the default value of 0.65 is used.
     * <p>
     * The method logs the parsed or default value for "minimalFrequency".
     *
     * @param configuration A {@link Map} containing the configuration parameters.
     * @return The parsed minimal frequency value as a double.
     * @throws NumberFormatException If the "minimalFrequency" value is not within the range (0, 1].
     */
    private double parseMinimalFrequency(Map<String, Object> configuration) {
        double minimalFrequency = 0.65; // Default value for minimal frequency
        if (configuration.containsKey("minimalFrequency")) {
            double value = (double) configuration.get("minimalFrequency");
            if (value > 0 && value <= 1) {
                minimalFrequency = value;
                Logging.logConfig("`minimalFrequency` set to %.2f.".formatted(minimalFrequency));
            } else {
                throw new NumberFormatException("Value for `minimalFrequency` must be in (0, 1], not %s.".formatted(value));
            }
        } else {
            Logging.logConfig("No value for `minimalFrequency` specified; defaulting to 0.65.");
        }
        return minimalFrequency;
    }

    /**
     * Parses the "storeFiltered" value from the configuration map.
     * <p>
     * This method retrieves the "storeFiltered" value from the provided configuration map. If the value is present, it is validated to
     * ensure it is either "true" or "false" (case-insensitive). If valid, the value is parsed as a boolean. If the value is not specified,
     * the default value of false is used.
     * <p>
     * The method logs the parsed or default value for "storeFiltered". If the value is invalid, a warning is logged.
     *
     * @param configuration A {@link Map} containing the configuration parameters.
     * @return The parsed "storeFiltered" value as a boolean.
     */
    private boolean parseStoreFiltered(Map<String, Object> configuration) {
        boolean storeFiltered = false; // Default value for storeFiltered
        if (configuration.containsKey("storeFiltered")) {
            storeFiltered = (boolean) configuration.get("storeFiltered");
            Logging.logConfig("`storeFiltered` set to %s.".formatted(storeFiltered));
        } else {
            Logging.logConfig("No value for `storeFiltered` specified; defaulting to false.");
        }
        return storeFiltered;
    }

    /**
     * Parses the "skipAnnotation" value from the configuration map.
     * <p>
     * This method retrieves the "skipAnnotation" value from the provided configuration map. If the value is present, it is validated to
     * ensure it is either "true" or "false" (case-insensitive). If valid, the value is parsed as a boolean. If the value is not specified,
     * the default value of false is used.
     * <p>
     * The method logs the parsed or default value for "skipAnnotation". If the value is invalid, a warning is logged.
     *
     * @param configuration A {@link Map} containing the configuration parameters.
     * @return The parsed "skipAnnotation" value as a boolean.
     */
    private boolean parseSkipAnnotation(Map<String, Object> configuration) {
        boolean skipAnnotation = false; // Default value for skipAnnotation
        if (configuration.containsKey("skipAnnotation")) {
            skipAnnotation = (boolean) configuration.get("skipAnnotation");
            Logging.logConfig("`skipAnnotation` set to %s.".formatted(skipAnnotation));
        } else {
            Logging.logConfig("No value for `skipAnnotation` specified; defaulting to false.");
        }
        return skipAnnotation;
    }

    /**
     * Parses the "skipTyping" value from the configuration map.
     * <p>
     * This method retrieves the "skipTyping" value from the provided configuration map. If the value is present, it is validated to ensure
     * it is either "true" or "false" (case-insensitive). If valid, the value is parsed as a boolean. If the value is not specified, the
     * default value of false is used.
     * <p>
     * The method logs the parsed or default value for "skipTyping". If the value is invalid, a warning is logged.
     *
     * @param configuration A {@link Map} containing the configuration parameters.
     * @return The parsed "skipTyping" value as a boolean.
     */
    private boolean parseSkipTyping(Map<String, Object> configuration) {
        boolean skipTyping = false; // Default value for skipTyping
        if (configuration.containsKey("skipTyping")) {
            skipTyping = (boolean) configuration.get("skipTyping");
            Logging.logConfig("`skipTyping` set to %s.".formatted(skipTyping));
        } else {
            Logging.logConfig("No value for `skipTyping` specified; defaulting to false.");
        }
        return skipTyping;
    }

    /**
     * Parses the masked positions from the configuration map.
     * <p>
     * This method reads a BED file specified in the "mask" parameter of the configuration map and extracts the masked positions. The
     * positions are stored in a map where the key is the contig ID, and the value is a set of excluded positions. If the "mask" parameter
     * is not specified or the file is invalid, an empty map is returned.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Checks if the "mask" parameter exists in the configuration map.</li>
     *     <li>Validates the file path to ensure it is a regular, non-empty file.</li>
     *     <li>Reads the BED file line by line, decodes each line into a {@link BEDFeature}, and extracts
     *         the start and end positions.</li>
     *     <li>Populates the map with the extracted positions, grouped by contig ID.</li>
     * </ul>
     * <p>
     * If the file is invalid or empty, a {@link MusialException} is thrown.
     *
     * @param configuration A {@link Map} containing the configuration parameters. The "mask" parameter specifies the path to the BED file.
     * @return A {@link Map} where the key is the contig ID, and the value is a set of excluded positions.
     * @throws IOException     If an I/O error occurs while reading the file.
     * @throws MusialException If the "mask" file is invalid or not a regular file.
     */
    private Map<String, Set<Integer>> parseMaskedPositions(Map<String, Object> configuration) throws IOException, MusialException {
        // Initialize a map to store excluded positions, with the contig id as the key.
        Map<String, Set<Integer>> maskedPositions = new HashMap<>();
        BEDCodec codec = new BEDCodec();
        int m = 0;

        // Check if a mask file is specified in the configuration; if not, return an empty map.
        if (!configuration.containsKey("mask"))
            return maskedPositions;

        Path path = Path.of((String) configuration.get("mask"));
        // Check if the path is null or blank; if so, return an empty map.
        if (PathUtils.isRegularFile(path) && !PathUtils.isDirectory(path) && !PathUtils.isEmptyFile(path)) {
            File file = path.toFile();
            try (BufferedReader reader = new BufferedReader(new FileReader(file, StandardCharsets.UTF_8))) {
                String line = reader.readLine();
                while (line != null) {
                    BEDFeature bedFeature = codec.decode(line);
                    if (bedFeature != null) {
                        maskedPositions
                                .computeIfAbsent(bedFeature.getContig(), k -> new HashSet<>())
                                .addAll(IntStream.rangeClosed(bedFeature.getStart(), bedFeature.getEnd())
                                        .boxed()
                                        .collect(Collectors.toSet()));
                        m += bedFeature.getEnd() - bedFeature.getStart() + 1;
                    }
                    line = reader.readLine();
                }
            }
            Logging.logConfig("`mask` set to %d positions.".formatted(m));
        } else {
            throw new MusialException("File specified for `mask` %s is empty or no regular file.".formatted(path));
        }
        // Return the populated map of excluded positions.
        return maskedPositions;
    }

    /**
     * Parses the reference sequence file from the configuration map.
     * <p>
     * This method retrieves the "reference" parameter from the configuration map and validates the file path. If the file exists, is not a
     * directory, and is not empty, it creates an {@link IndexedFastaSequenceFile} using the specified file. If the "reference" parameter is
     * not provided, the method returns null.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Checks if the "reference" parameter exists in the configuration map.</li>
     *     <li>Validates the file path to ensure it is a regular, non-empty file.</li>
     *     <li>Creates and returns an {@link IndexedFastaSequenceFile} for the specified file.</li>
     * </ul>
     * <p>
     * If the file is invalid, a {@link MusialException} is thrown.
     *
     * @param configuration A {@link Map} containing the configuration parameters. The "reference" parameter specifies the path to the
     *                      reference file.
     * @return A {@link IndexedFastaSequenceFile} object representing the reference sequence file, or null if not specified.
     * @throws IOException     If an I/O error occurs while accessing the file.
     * @throws MusialException If the "reference" file is invalid or not a regular file.
     */
    private IndexedFastaSequenceFile parseReference(Map<String, Object> configuration) throws IOException, MusialException {
        // Check if a reference file is specified in the configuration; if not, return null.
        if (!configuration.containsKey("reference")) {
            return null;
        }

        Path path = Path.of((String) configuration.get("reference"));
        if (PathUtils.isRegularFile(path) && !PathUtils.isDirectory(path) && !PathUtils.isEmptyFile(path)) {
            File file = path.toFile();
            return new IndexedFastaSequenceFile(file.toPath(),
                    FastaSequenceIndexCreator.buildFromFasta(file.toPath()));
        } else {
            throw new MusialException("File specified for `reference` %s is empty or no regular file.".formatted(path));
        }
    }

    /**
     * Parses the annotation file from the configuration map.
     * <p>
     * This method retrieves the "annotation" parameter from the configuration map and validates the file path. If the file exists, is not a
     * directory, and is not empty, it reads the file using the {@link GFF3Reader} and parses it into a {@link FeatureList}. If the
     * "annotation" parameter is not provided, the method returns an empty list.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Checks if the "annotation" parameter exists in the configuration map.</li>
     *     <li>Validates the file path to ensure it is a regular, non-empty file.</li>
     *     <li>Reads and parses the file into a {@link FeatureList}.</li>
     * </ul>
     * <p>
     * If the file is invalid, a {@link MusialException} is thrown.
     *
     * @param configuration A {@link Map} containing the configuration parameters. The "annotation" parameter specifies the path to the
     *                      annotation file.
     * @return A {@link FeatureList} object representing the parsed annotation features; the list may be empty.
     * @throws IOException     If an I/O error occurs while accessing the file.
     * @throws MusialException If the "annotation" file is invalid or not a regular file.
     */
    private FeatureList parseAnnotation(Map<String, Object> configuration) throws IOException, MusialException {
        // Check if an annotation file is specified in the configuration; if not, return null.
        if (!configuration.containsKey("annotation")) {
            return new FeatureList();
        }

        Path path = Path.of((String) configuration.get("annotation"));
        if (PathUtils.isRegularFile(path) && !PathUtils.isDirectory(path) && !PathUtils.isEmptyFile(path)) {
            File file = path.toFile();
            FeatureList featureList = GFF3Reader.read(file.getCanonicalPath());
            return featureList;
        } else {
            throw new MusialException("File specified for `annotation` %s is empty or no regular file.".formatted(path));
        }
    }

    /**
     * Parses the features from the configuration map.
     * <p>
     * This method retrieves the "features" parameter from the configuration map and validates the file path. If the file exists, is not a
     * directory, and is not empty, it reads the file as a nested map using the {@link IO#readTabularFileAsNestedMap(File)} utility method.
     * If the "features" parameter is not provided or the file is invalid, an empty map is returned.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Checks if the "features" parameter exists in the configuration map.</li>
     *     <li>Validates the file path to ensure it is a regular, non-empty file.</li>
     *     <li>Reads and parses the file into a nested map.</li>
     * </ul>
     * <p>
     * If the file is invalid, an empty map is returned.
     *
     * @param configuration A {@link Map} containing the configuration parameters. The "features" parameter specifies the path to the
     *                      features file.
     * @return A {@link Map} where the key is the feature category, and the value is another map containing feature attributes and their
     * values.
     * @throws IOException     If an I/O error occurs while accessing the file.
     * @throws MusialException If the features do not contain at least a `key` and `value` attribute.
     */
    private Map<String, Map<String, String>> parseFeatures(Map<String, Object> configuration) throws IOException, MusialException {
        Map<String, Map<String, String>> features = Collections.emptyMap();
        // Check if an annotation file is specified in the configuration; if not, return an empty map.
        if (configuration.containsKey("features")) {
            Path path = Path.of((String) configuration.get("features"));
            if (PathUtils.isRegularFile(path) && !PathUtils.isDirectory(path) && !PathUtils.isEmptyFile(path)) {
                features = IO.readTabularFileAsNestedMap(path.toFile());
            }
            if (features.values().stream().anyMatch(m -> !m.containsKey("key") || !m.containsKey("value")))
                throw new MusialException("Each feature in `features` must contain at least a `key` and `value` attribute.");
        }
        return features;
    }

}
