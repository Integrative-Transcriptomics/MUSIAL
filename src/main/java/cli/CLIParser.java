package cli;

import com.fasterxml.jackson.databind.JsonNode;
import com.github.fge.jackson.JsonLoader;
import com.github.fge.jsonschema.core.exceptions.ProcessingException;
import com.github.fge.jsonschema.core.report.ProcessingReport;
import com.github.fge.jsonschema.main.JsonSchema;
import com.github.fge.jsonschema.main.JsonSchemaFactory;
import com.google.gson.Gson;
import exceptions.MusialException;
import main.Musial;
import main.Tasks;
import org.apache.commons.cli.*;
import org.apache.commons.io.IOUtils;
import utility.Validation;

import java.io.*;
import java.lang.reflect.Array;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Objects;

/**
 * Parses command line interface arguments.
 * <p>
 * Used to parse and validate passed command line options or print help information to the user.
 *
 * @author Simon Hackl
 * @version 2.3
 * @since 2.1
 */
@SuppressWarnings("DuplicatedCode")
public class CLIParser {

    /**
     * Description of single MUSIAL tasks used to construct help message. If new tasks are implemented, a description should be supplemented here.
     */
    @SuppressWarnings("FieldCanBeLocal")
    private final String musialTaskDescription = """
            \n
             build            Build a local database/MUSIAL storage (brotli compressed binary JSON format) from variant calls.
             view_features    View annotated features from a built MUSIAL storage in a tabular format.
             view_samples     View annotated samples from a built MUSIAL storage in a tabular format.
             view_variants    View annotated variants from a built MUSIAL storage in a tabular format.
             export_table     Export variants from a built MUSIAL storage into a matrix-like .tsv file.
             export_sequence  Generate sequences in .fasta format from a built MUSIAL storage.
            \n
            """;

    /**
     * Internal representation of the build task JSON configuration file.
     */
    public HashMap<String, Object> buildConfiguration;

    /**
     * Arguments parsed from command line.
     */
    public CommandLine arguments;

    /**
     * Constructor of the {@link CLIParser} class.
     *
     * @param args {@link String} {@link Array} of passed command line arguments.
     */
    public CLIParser(String[] args) throws MusialException {
        // Init. `Option` object.
        Options options = new Options();

        // Add explicit options dependent on task to execute.
        if (Musial.TASK.equalsIgnoreCase(Tasks.BUILD.name())) {
            addBuildOptions(options);
        } else if (Musial.TASK.equalsIgnoreCase(Tasks.VIEW_FEATURES.name())) {
            addViewFeaturesOptions(options);
        } else if (Musial.TASK.equalsIgnoreCase(Tasks.VIEW_SAMPLES.name())) {
            addViewSamplesOptions(options);
        } else if (Musial.TASK.equalsIgnoreCase(Tasks.VIEW_VARIANTS.name())) {
            addViewVariantsOptions(options);
        } else if (Musial.TASK.equalsIgnoreCase(Tasks.EXPORT_TABLE.name())) {
            addExportTableOptions(options);
        } else if (Musial.TASK.equalsIgnoreCase(Tasks.EXPORT_SEQUENCE.name())) {
            addExportSequenceOptions(options);
        } else {
            exitNotRecognized(args[0]);
        }
        // Init. formatter for help message and default command line parser.
        HelpFormatter helpformatter = new HelpFormatter();
        CommandLineParser parser = new DefaultParser();
        // Check if help argument was specified.
        boolean isFirstArgument = true;
        for (String arg : args) {
            if (Objects.equals(arg, "-h") || Objects.equals(arg, "--help")) {
                String headerText;
                if (isFirstArgument) {
                    headerText =
                            "MUSIAL provides distinct tasks to computes prokaryotic genome, gene and protein sequence alignments from variantInformation call datasets from multiple samples of one species."
                                    + " Available tasks are:"
                                    + musialTaskDescription
                                    + "Call `java -jar " + Musial.NAME + "-" + Musial.VERSION + ".jar <task> [-h|--help]` for more information.";
                } else {
                    headerText = "command line arguments of task " + Musial.TASK;
                }
                helpformatter.printHelp(
                        150,
                        "java -jar " + Musial.NAME + "-" + Musial.VERSION + ".jar",
                        headerText,
                        options,
                        "",
                        true
                );
                System.exit(0);
            }
            isFirstArgument = false;
        }
        // Else, command line interface arguments are parsed wrt. specified task.
        try {
            this.arguments = parser.parse(options, args);
            if (Musial.TASK.equalsIgnoreCase(Tasks.BUILD.name())) {
                this.buildConfiguration = validateBuildOptions();
            } else if (Musial.TASK.equalsIgnoreCase(Tasks.VIEW_FEATURES.name())) {
                validateViewFeaturesOptions();
            } else if (Musial.TASK.equalsIgnoreCase(Tasks.VIEW_SAMPLES.name())) {
                validateViewSamplesOptions();
            } else if (Musial.TASK.equalsIgnoreCase(Tasks.VIEW_VARIANTS.name())) {
                validateViewVariantsOptions();
            } else if (Musial.TASK.equalsIgnoreCase(Tasks.EXPORT_TABLE.name())) {
                validateExportTableOptions();
            } else if (Musial.TASK.equalsIgnoreCase(Tasks.EXPORT_SEQUENCE.name())) {
                validateExportSequenceOptions();
            }
        } catch (ParseException | IOException | ProcessingException e) {
            throw new MusialException("(CLI Argument Parsing Failed) " + e.getMessage());
        }
    }

    /**
     * Fills in cli options of the build task.
     *
     * @param options {@link Option} instance to add options to.
     */
    private void addBuildOptions(Options options) {
        options.addOption(
                Option.builder("c")
                        .longOpt("configuration")
                        .desc("Path to a .json file specifying the build configuration for MUSIAL. Please visit https://github.com/Integrative-Transcriptomics/MUSIAL for a detailed explanation on how to specify the MUSIAL build configuration file.")
                        .hasArg()
                        .required()
                        .build()
        );
    }

    /**
     * Fills in cli options of the view features task.
     *
     * @param options {@link Option} instance to add options to.
     */
    private void addViewFeaturesOptions(Options options) {
        options.addOption(
                Option.builder("I")
                        .longOpt("storage")
                        .desc("Path to a .json file generated with the build task of MUSIAL to view.")
                        .hasArg()
                        .required()
                        .build()
        );
        Option optionFeatures = Option.builder("f")
                .longOpt("features")
                .desc("Explicit space separated list of features to view (Default: all).")
                .hasArgs()
                .build();
        optionFeatures.setArgs(Option.UNLIMITED_VALUES);
        options.addOption(optionFeatures);
        options.addOption(
                Option.builder("o")
                        .longOpt("output")
                        .desc("Path to a file to write the output to (Default: stdout).")
                        .hasArg()
                        .build()
        );
    }

    /**
     * Fills in cli options of the view samples task.
     *
     * @param options {@link Option} instance to add options to.
     */
    private void addViewSamplesOptions(Options options) {
        options.addOption(
                Option.builder("I")
                        .longOpt("storage")
                        .desc("Path to a .json file generated with the BUILD task of MUSIAL to view.")
                        .hasArg()
                        .required()
                        .build()
        );
        Option optionSamples = Option.builder("s")
                .longOpt("samples")
                .desc("Explicit space separated list of samples to view (Default: all).")
                .hasArgs()
                .build();
        optionSamples.setArgs(Option.UNLIMITED_VALUES);
        options.addOption(optionSamples);
        options.addOption(
                Option.builder("o")
                        .longOpt("output")
                        .desc("Path to a file to write the output to (Default: stdout).")
                        .hasArg()
                        .build()
        );
    }

    /**
     * Fills in cli options of the view variants task.
     *
     * @param options {@link Option} instance to add options to.
     */
    private void addViewVariantsOptions(Options options) {
        options.addOption(
                Option.builder("I")
                        .longOpt("storage")
                        .desc("Path to a .json file generated with the BUILD task of MUSIAL to view.")
                        .hasArg()
                        .required()
                        .build()
        );
        Option optionFeatures = Option.builder("f")
                .longOpt("features")
                .desc("Explicit space separated list of features to restrict variants to (Default: all).")
                .hasArgs()
                .build();
        optionFeatures.setArgs(Option.UNLIMITED_VALUES);
        options.addOption(optionFeatures);
        Option optionSamples = Option.builder("s")
                .longOpt("samples")
                .desc("Explicit space separated list of samples to restrict variants to (Default: all).")
                .hasArgs()
                .build();
        optionSamples.setArgs(Option.UNLIMITED_VALUES);
        options.addOption(optionSamples);
        options.addOption(
                Option.builder("c")
                        .longOpt("content")
                        .desc("Sets the content type of viewed variants; One of `nucleotide` or `aminoacid` (Default: nucleotide).")
                        .hasArg()
                        .build()
        );
        options.addOption(
                Option.builder("o")
                        .longOpt("output")
                        .desc("Path to a file to write the output to (Default: stdout).")
                        .hasArg()
                        .build()
        );
    }

    /**
     * Fills in cli options of the export table task.
     *
     * @param options {@link Option} instance to add options to.
     */
    private void addExportTableOptions(Options options) {
        options.addOption(
                Option.builder("I")
                        .longOpt("storage")
                        .desc("Path to a .json file generated with the BUILD task of MUSIAL to view.")
                        .hasArg()
                        .required()
                        .build()
        );
        options.addOption(
                Option.builder("O")
                        .longOpt("output")
                        .desc("Path to a file to write the output to.")
                        .hasArg()
                        .required()
                        .build()
        );
        Option optionFeatures = Option.builder("F")
                .longOpt("feature")
                .desc("The feature for which variants should be exported.")
                .hasArgs()
                .required()
                .build();
        optionFeatures.setArgs(Option.UNLIMITED_VALUES);
        options.addOption(optionFeatures);
        Option optionSamples = Option.builder("s")
                .longOpt("samples")
                .desc("Explicit space separated list of samples to restrict variants to (Default: all).")
                .hasArgs()
                .build();
        optionSamples.setArgs(Option.UNLIMITED_VALUES);
        options.addOption(optionSamples);
        options.addOption(
                Option.builder("c")
                        .longOpt("content")
                        .desc("Sets the content type of viewed variants; One of `nucleotide` or `aminoacid` (Default: nucleotide).")
                        .hasArg()
                        .build()
        );
    }

    /**
     * Fills in cli options of the export sequence task.
     *
     * @param options {@link Option} instance to add options to.
     */
    private void addExportSequenceOptions(Options options) {
        options.addOption(
                Option.builder("I")
                        .longOpt("storage")
                        .desc("Path to a .json file generated with the BUILD task of MUSIAL to view.")
                        .hasArg()
                        .required()
                        .build()
        );
        options.addOption(
                Option.builder("O")
                        .longOpt("output")
                        .desc("Path to a file to write the output to.")
                        .hasArg()
                        .required()
                        .build()
        );
        Option optionFeatures = Option.builder("F")
                .longOpt("feature")
                .desc("The feature for which variants should be exported.")
                .hasArgs()
                .required()
                .build();
        optionFeatures.setArgs(Option.UNLIMITED_VALUES);
        options.addOption(optionFeatures);
        options.addOption(
                Option.builder("k")
                        .longOpt("conserved")
                        .desc("Export conserved sites (Default: Only variantInformation sites).")
                        .build()
        );
        options.addOption(
                Option.builder("r")
                        .longOpt("reference")
                        .desc("Path to a .fasta file yielding the reference sequences with which the specified MUSIAL storage file was built. If the file is not indexed, this wil be done automatically. This option is only required for `content=nucleotide` and `conserved`.")
                        .hasArg()
                        .build()
        );
        Option optionSamples = Option.builder("s")
                .longOpt("samples")
                .desc("Explicit space separated list of samples to restrict variants to (Default: all).")
                .hasArgs()
                .build();
        optionSamples.setArgs(Option.UNLIMITED_VALUES);
        options.addOption(optionSamples);
        options.addOption(
                Option.builder("c")
                        .longOpt("content")
                        .desc("Sets the content type of viewed variants; One of `nucleotide` or `aminoacid` (Default: nucleotide).")
                        .hasArg()
                        .build()
        );
        options.addOption(
                Option.builder("a")
                        .longOpt("aligned")
                        .desc("Whether to align sequences (Default: No alignment).")
                        .build()
        );
        options.addOption(
                Option.builder("m")
                        .longOpt("merge")
                        .desc("Whether to merge samples by alleles and proteoforms (Default: No merging).")
                        .build()
        );
    }

    /**
     * Validates the configuration JSON file passed for the build task against the build configuration JSON schema.
     *
     * @return A {@link HashMap} storing the parsed build configuration.
     * @throws IOException         If any I/O operation during the validation process fails.
     * @throws MusialException     If the JSON schema validation fails.
     * @throws ProcessingException If the JSON schema setup fails.
     */
    private HashMap<String, Object> validateBuildOptions() throws IOException, MusialException, ProcessingException {
        // Validate passed build configuration file against resp. schema.
        try (FileOutputStream buildConfigurationSchemaOutputStream = new FileOutputStream("./musial_build_configuration.schema.json")) {
            // Copy build configuration file schema from class resources to current directory.
            final JsonSchemaFactory factory = JsonSchemaFactory.byDefault();
            IOUtils.copy(Objects.requireNonNull(Musial.class.getResourceAsStream("/musial_build_configuration.schema.json")), buildConfigurationSchemaOutputStream);
            final JsonNode buildConfigurationJsonSchema = JsonLoader.fromFile(new File("./musial_build_configuration.schema.json"));
            final JsonSchema buildConfigurationSchema = factory.getJsonSchema(buildConfigurationJsonSchema);
            // Validate passed configuration file with schema.
            final JsonNode buildConfiguration = JsonLoader.fromPath(arguments.getOptionValue("c"));
            ProcessingReport schemaValidationReport = buildConfigurationSchema.validate(buildConfiguration);
            if (!schemaValidationReport.isSuccess()) {
                throw new MusialException("(CLI Argument Parsing Failed) Validation of build configuration file failed; " + schemaValidationReport);
            }
        } finally {
            //noinspection ResultOfMethodCallIgnored
            new File("./musial_build_configuration.schema.json").delete();
        }
        // Parse specified build configuration into JSON object.
        try (
                BufferedReader bufferedReader =
                        new BufferedReader(
                                new InputStreamReader(
                                        Files.newInputStream(Path.of(arguments.getOptionValue("c"))), StandardCharsets.UTF_8
                                )
                        )
        ) {
            //noinspection unchecked
            return new Gson().fromJson(bufferedReader, HashMap.class);
        } catch (Exception e) {
            throw new MusialException(
                    "(CLI Argument Parsing Failed) " + e.getMessage());
        }
    }

    /**
     * Validates all user options specified for the view features task.
     *
     * @throws MusialException If the validation fails.
     */
    private void validateViewFeaturesOptions() throws MusialException {
        String EXCEPTION_PREFIX = "(Task `view_features` Argument Validation)";
        if (!Validation.isFile(new File(arguments.getOptionValue("I"))))
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("I") + "' for parameter `storage`; failed to access file.");
        if (arguments.hasOption("o")) {
            if (Validation.isFile(new File(arguments.getOptionValue("o"))))
                throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("o") + "' for parameter `output`; specified file already exists.");
        }
    }

    /**
     * Validates all user options specified for the view samples task.
     *
     * @throws MusialException If the validation fails.
     */
    private void validateViewSamplesOptions() throws MusialException {
        String EXCEPTION_PREFIX = "(Task `view_samples` Argument Validation)";
        if (!Validation.isFile(new File(arguments.getOptionValue("I"))))
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("I") + "' for parameter `storage`; failed to access file.");
        if (arguments.hasOption("o")) {
            if (Validation.isFile(new File(arguments.getOptionValue("o"))))
                throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("o") + "' for parameter `output`; specified file already exists.");
        }
    }

    /**
     * Validates all user options specified for the view variants task.
     *
     * @throws MusialException If the validation fails.
     */
    private void validateViewVariantsOptions() throws MusialException {
        String EXCEPTION_PREFIX = "(Task `view_variants` Argument Validation)";
        if (!Validation.isFile(new File(arguments.getOptionValue("I"))))
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("I") + "' for parameter `storage`; failed to access file.");
        if (arguments.hasOption("o")) {
            if (Validation.isFile(new File(arguments.getOptionValue("o"))))
                throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("o") + "' for parameter `output`; specified file already exists.");
        }
        if (arguments.hasOption("c")) {
            if (!Objects.equals(arguments.getOptionValue("c"), "nucleotide") && !Objects.equals(arguments.getOptionValue("c"), "aminoacid"))
                throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("c") + "' for parameter `content`; should be either `nucleotide` or `aminoacid`.");
        }

    }

    /**
     * Validates all user options specified for the export table task.
     *
     * @throws MusialException If the validation fails.
     */
    private void validateExportTableOptions() throws MusialException {
        String EXCEPTION_PREFIX = "(Task `export_table` Argument Validation)";
        if (!Validation.isFile(new File(arguments.getOptionValue("I"))))
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("I") + "' for parameter `storage`; failed to access file.");
        if (Validation.isFile(new File(arguments.getOptionValue("O"))))
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("O") + "' for parameter `output`; specified file already exists.");
        if (arguments.hasOption("c")) {
            if (!Objects.equals(arguments.getOptionValue("c"), "nucleotide") && !Objects.equals(arguments.getOptionValue("c"), "aminoacid"))
                throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("c") + "' for parameter `content`; should be either `nucleotide` or `aminoacid`.");
        }
        if (arguments.hasOption("k") && (!arguments.hasOption("c") || arguments.getOptionValue("c").equals("nucleotide")) && !arguments.hasOption("r"))
            throw new MusialException(EXCEPTION_PREFIX + " In order to export non conserved sites in content mode 'nucleotide' a reference sequence has to be provided.");
        if (arguments.hasOption("r") && !Validation.isFile(new File(arguments.getOptionValue("r"))))
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("r") + "' for parameter `reference`; failed to access file.");
    }

    /**
     * Validates all user options specified for the export sequence task.
     *
     * @throws MusialException If the validation fails.
     */
    private void validateExportSequenceOptions() throws MusialException {
        String EXCEPTION_PREFIX = "(Task `export_sequence` Argument Validation)";
        if (!Validation.isFile(new File(arguments.getOptionValue("I"))))
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("I") + "' for parameter `storage`; failed to access file.");
        if (Validation.isFile(new File(arguments.getOptionValue("O"))))
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("O") + "' for parameter `output`; specified file already exists.");
        if (arguments.hasOption("c")) {
            if (!Objects.equals(arguments.getOptionValue("c"), "nucleotide") && !Objects.equals(arguments.getOptionValue("c"), "aminoacid"))
                throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("c") + "' for parameter `content`; should be either `nucleotide` or `aminoacid`.");
        }
        if (arguments.hasOption("k") && (!arguments.hasOption("c") || arguments.getOptionValue("c").equals("nucleotide")) && !arguments.hasOption("r"))
            throw new MusialException(EXCEPTION_PREFIX + " In order to export non conserved sites in content mode 'nucleotide' a reference sequence has to be provided.");
        if (arguments.hasOption("r") && !Validation.isFile(new File(arguments.getOptionValue("r"))))
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + arguments.getOptionValue("r") + "' for parameter `reference`; failed to access file.");
    }

    /**
     * Exits the program, if an unknown task was specified.
     *
     * @param taskName The task name supplied by the user, i.e., the first cli argument.
     */
    private void exitNotRecognized(String taskName) {
        if (!Musial.TASK.equals("-h") && !Musial.TASK.equals("--help")) {
            System.out.println("Task `" + taskName + "` not recognized. Call `java -jar " + Musial.NAME + "-" + Musial.VERSION + ".jar [-h|--help]` for more information.");
            System.exit(0);
        }
    }

}
