package utility;

import com.fasterxml.jackson.databind.JsonNode;
import com.github.fge.jackson.JsonLoader;
import com.github.fge.jsonschema.core.exceptions.ProcessingException;
import com.github.fge.jsonschema.core.report.ProcessingReport;
import com.github.fge.jsonschema.main.JsonSchema;
import com.github.fge.jsonschema.main.JsonSchemaFactory;
import com.google.gson.Gson;
import exceptions.MusialException;
import main.Musial;
import org.apache.commons.cli.*;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Objects;

/**
 * Command Line Interface (CLI) utility class for parsing and validating command-line arguments.
 * <p>
 * This final class provides functionality to handle command-line arguments for various tasks
 * (e.g., build, update, view, sequence) in the MUSIAL application. It defines options, parses
 * arguments, validates them, and stores the parsed parameters for further processing.
 * </p>
 *
 * <p>
 * The CLI class supports displaying help messages, validating input files, and managing task-specific
 * logic through nested task classes. It ensures proper error handling and user guidance for invalid
 * or incomplete arguments.
 * </p>
 */
public final class CLI {

    /**
     * A map to store parsed command-line parameters.
     * <p>
     * This {@link HashMap} holds the validated and parsed parameters from the command-line arguments.
     * The keys represent parameter names, and the values represent their corresponding values.
     * The content of this map depends on the task being executed (e.g., build, update, view, sequence).
     * </p>
     */
    public static HashMap<String, Object> parameters = new HashMap<>();

    /**
     * Stores the command-line options for the CLI.
     * <p>
     * This {@link Options} object is used to define and manage the command-line options
     * available for the application. It includes options specific to the tasks
     * (e.g., build, update, view, sequence) and is populated during the initialization
     * of the {@link CLI} instance.
     * </p>
     */
    private static final Options options = new Options();

    /**
     * Stores the parsed command-line arguments.
     * <p>
     * This {@link CommandLine} object holds the arguments passed to the application
     * after being parsed and validated. It is used to access the values of the
     * specified options and arguments for further processing.
     * </p>
     */
    private static CommandLine arguments;

    public static void parse(String[] args) throws MusialException {
        // Fill options dependent on task to execute.
        switch (Musial.task) {
            case BUILD -> CLI.Build.options();
            case EXPAND -> CLI.Expand.options();
            case VIEW -> CLI.View.options();
            case SEQUENCE -> CLI.Sequence.options();
        }

        // Initialize help message formatter and default command line parser.
        HelpFormatter helpformatter = new HelpFormatter();
        CommandLineParser parser = new DefaultParser();

        // Check if a help argument (-h/--help) was specified.
        if (Arrays.stream(args).anyMatch(arg -> arg.equals("-h") || arg.equals("--help"))) {
            String helpText = """

                    MUSIAL aggregates variant calls from multiple samples of one prokaryotic species and provides an interface to compute a comprehensive overview as well as genome, gene and protein sequence alignments.

                    Available tasks are:
                    \033[47m\033[1;30m build    \033[0m : Build a local database file (storage) in JSON format from variant calls. This is mandatory to execute other tasks.
                    \033[47m\033[1;30m expand   \033[0m : Expand (or preview) an existing storage with sample data from variant call files.
                    \033[47m\033[1;30m view     \033[0m : View the content, i.e., entries (features, samples or variants) and their attributes, of a MUSIAL storage file.
                    \033[47m\033[1;30m sequence \033[0m : Export sequences of features from a MUSIAL storage file.

                    Call `java -jar %s-%s.jar <task> [-h|--help]` for more information.
                    """.formatted(Musial.softwareName, Musial.softwareVersion);

            // If a task was specified, but is not recognized, adjust the help text.
            if (Musial.task.equals(Musial.Task.UNDEFINED)) {
                helpText += """

                        Task \033[1;31m%s\033[0m not recognized.
                        """.formatted(args[0]);
            } else {
                // If a task was specified, add the command line arguments for that task.
                helpText += """

                        \033[47m\033[1;30m Command line arguments of task %s \033[0m
                        """.formatted(Musial.task.toString().toLowerCase());
            }

            // Print the help message with the specified options and arguments.
            helpformatter.printHelp(
                    150,
                    "java -jar %s-%s.jar %s".formatted(Musial.softwareName, Musial.softwareVersion,
                            Musial.task.equals(Musial.Task.UNDEFINED) ? "<task>" : Musial.task.toString().toLowerCase()),
                    helpText,
                    options,
                    "",
                    true
            );

            // Exit the program after displaying the help message.
            System.exit(0);
        }

        // Else, command line interface arguments are parsed wrt. specified task.
        try {
            arguments = parser.parse(options, args);
            switch (Musial.task) {
                case BUILD -> CLI.Build.transfer();
                case EXPAND -> CLI.Expand.transfer();
                case VIEW -> CLI.View.transfer();
                case SEQUENCE -> CLI.Sequence.transfer();
                default -> exitNotRecognized(args);
            }
        } catch (ParseException | IOException | ProcessingException e) {
            throw new MusialException("Command line argument parsing failed; %s".formatted(e.getMessage()));
        }
    }

    /**
     * Exits the program if the specified task is not recognized and no help option is provided.
     * <p>
     * This method checks if the command-line arguments do not include the help option (`-h` or `--help`).
     * If the task is unrecognized, it prints an error message with the unrecognized task and instructions
     * on how to access the help message. The program then terminates with a call to {@link System#exit(int)}.
     * </p>
     *
     * @param args An array of {@link String} representing the command-line arguments.
     */
    private static void exitNotRecognized(String[] args) {
        if (Arrays.stream(args).noneMatch(arg -> arg.equals("-h") || arg.equals("--help"))) {
            System.out.printf("Task \033[1;31m%s\033[0m not recognized. Call `java -jar %s-%s.jar [-h|--help]` for more information.%n",
                    args[0], Musial.softwareName, Musial.softwareVersion);
            System.exit(0);
        }
    }

    /**
     * Represents a task interface for defining CLI options and validation logic.
     * <p>
     * This interface provides two static methods:
     * <ul>
     *   <li>{@code options()}: Defines the command-line options specific to a task.</li>
     *   <li>{@code transfer()}: Validates the parsed command-line arguments for the task.</li>
     * </ul>
     * </p>
     *
     * <p>
     * Implementing classes should provide concrete implementations for these methods
     * to handle task-specific logic for CLI parsing. This is intended to only comprise
     * syntax parsing and no validation logic.
     * </p>
     *
     * @noinspection unused
     */
    private interface Task {

        /**
         * Defines the command-line options for the task.
         */
        private static void options() {
        }

        /**
         * Transfer the parsed command-line arguments for the task to the {@link #parameters} property.
         */
        private static void transfer() {

        }

    }

    /**
     * Handles the Build task for the CLI.
     * <p>
     * This class defines the command-line options and validation logic for the Build task.
     * It allows users to specify a JSON file containing the task parameters for MUSIAL.
     * </p>
     */
    private static class Build implements Task {

        /**
         * Defines the command-line options for the Build task.
         * <p>
         * This method adds the following option:
         * <ul>
         *   <li>`-P` or `--parameters`: Path to a `.json` file specifying the build task parameters for MUSIAL.</li>
         * </ul>
         * The `-P` option is required.
         * </p>
         */
        private static void options() {
            options.addOption(Option.builder("C")
                    .longOpt("configuration")
                    .desc("Path to a JSON file specifying the build task parameter configuration for MUSIAL. Visit the documentation for details.")
                    .hasArg()
                    .required()
                    .build());
        }

        /**
         * Transfers the command-line arguments for the build task.
         * <p>
         * This method performs the following steps:
         * <ul>
         *   <li>Validates the build parameters file against a predefined JSON schema located at `/buildConfigurationSchema.json`.</li>
         *   <li>Parses the validated file into a {@link HashMap} for further processing.</li>
         * </ul>
         * If the validation fails, a {@link MusialException} is thrown with the validation report.
         * If parsing fails, a {@link MusialException} is thrown with the error details.
         * </p>
         *
         * @throws IOException         If an I/O error occurs while reading the schema or the parameters file.
         * @throws MusialException     If the validation fails or the file cannot be parsed.
         * @throws ProcessingException If an error occurs during schema validation.
         */
        private static void transfer() throws IOException, MusialException, ProcessingException {
            // Validate the build configuration against schema.
            try (InputStream schemaStream = Objects.requireNonNull(Musial.class.getResourceAsStream("/buildConfigurationSchema.json"))) {
                JsonSchema schema = JsonSchemaFactory.byDefault()
                        .getJsonSchema(JsonLoader.fromReader(new InputStreamReader(schemaStream, StandardCharsets.UTF_8)));
                JsonNode config = JsonLoader.fromPath(arguments.getOptionValue("C"));
                ProcessingReport report = schema.validate(config);
                if (!report.isSuccess()) {
                    throw new MusialException("The specified build configuration file is invalid:\n%s".formatted(report.toString()));
                }
            }

            // Parse the build configuration into a JSON object.
            try (Reader reader = Files.newBufferedReader(Path.of(arguments.getOptionValue("C")), StandardCharsets.UTF_8)) {
                //noinspection unchecked
                parameters = new Gson().fromJson(reader, HashMap.class);
            } catch (Exception e) {
                throw new MusialException("Validation of task `build` CLI parameters failed; %s".formatted(e.getMessage()));
            }
        }

    }

    /**
     * Handles the Expand task for the CLI.
     * <p>
     * This class defines the command-line options and validation logic for the Expand task.
     * It allows users to specify input files, variant call files, sample annotations, and other options
     * for updating a MUSIAL storage file.
     * </p>
     */
    private static class Expand implements Task {

        /**
         * Defines the command-line options for the Expand task.
         * <p>
         * This method adds the following options:
         * <ul>
         *   <li>`-I` or `--input`: Path to the input `.json(.gz)` file generated with the build task.</li>
         *   <li>`-V` or `--files`: List of file or directory paths. Files must be in VCF format.</li>
         *   <li>`-m` or `--info`: Path to a `.tsv` or `.csv` file specifying sample annotations.</li>
         *   <li>`-o` or `--output`: Path to write the output file (default is to overwrite the input file).</li>
         *   <li>`-p` or `--preview`: Reports novel entries without writing the expanded storage to a file.</li>
         * </ul>
         * </p>
         */
        private static void options() {
            options.addOption(Option.builder("I")
                    .longOpt("storage")
                    .desc("Path to a .json(.gz) file generated with the build task of MUSIAL.")
                    .hasArg()
                    .required()
                    .build());
            options.addOption(Option.builder("V")
                    .longOpt("vcfInput")
                    .desc("List of file or directory paths. All files must be in VCF format.")
                    .hasArgs()
                    .required()
                    .build());
            options.addOption(Option.builder("m")
                    .longOpt("vcfMeta")
                    .desc("Path to a .tsv or .csv file specifying sample annotations.")
                    .hasArg()
                    .build());
            options.addOption(Option.builder("o")
                    .longOpt("output")
                    .desc("Path to write the output file (default: overwrite input file).")
                    .hasArg()
                    .build());
            options.addOption(Option.builder("p")
                    .longOpt("preview")
                    .desc("Only report on novel entries without writing the updated storage.")
                    .build());
        }

        /**
         * Transfers the command-line arguments for the update task.
         */
        private static void transfer() throws IOException {
            parameters = new HashMap<>();
            parameters.put("input", arguments.getOptionValue("I"));
            parameters.put("vcfInput", Arrays.asList(arguments.getOptionValues("V")));
            parameters.put("vcfMeta", arguments.getOptionValue("m"));
            parameters.put("output", arguments.getOptionValue("o", "overwrite"));
            parameters.put("write", !arguments.hasOption("p"));
        }

    }

    /**
     * Handles the view task for the CLI.
     * <p>
     * This class defines the command-line options and validation logic for the view task.
     * It allows users to specify input files, content types, filters, and output paths
     * for viewing the content of a MUSIAL storage file.
     * </p>
     */
    private static class View implements Task {

        /**
         * Defines the command-line options for the view task.
         * <p>
         * This method adds the following options:
         * <ul>
         *   <li>`-I` or `--input`: Path to the input `.json(.gz)` file generated with the build task.</li>
         *   <li>`-C` or `--content`: Specifies the content type (`features`, `samples`, `variants`, `alleles`, `sampleSequenceTypes`, or `variantCalls`).</li>
         *   <li>`-f` or `--filter`: List of feature/sample names or positions to filter the output (default is no filters).</li>
         *   <li>`-o` or `--output`: Path to write the output file (default is stdout).</li>
         * </ul>
         * </p>
         */
        private static void options() {
            options.addOption(Option.builder("I")
                    .longOpt("storage")
                    .desc("Path to a .json(.gz) file generated with the build task of MUSIAL.")
                    .hasArg()
                    .required()
                    .build());
            options.addOption(Option.builder("C")
                    .longOpt("content")
                    .desc("One of %s.".formatted(String.join(", ", Musial.View.content)))
                    .hasArg()
                    .required()
                    .build());
            options.addOption(Option.builder("f")
                    .longOpt("filter")
                    .desc("List of feature-, sample names, and/or positions for which the output is to be filtered (default: no filters). Entries may be ignored depending on the content.")
                    .hasArgs()
                    .build());
            options.addOption(Option.builder("o")
                    .longOpt("output")
                    .desc("Path to directory or file to write the output to (default: stdout).")
                    .hasArg()
                    .build());
        }

        /**
         * Transfers the command-line arguments for the view task.
         */
        private static void transfer() {
            parameters = new HashMap<>();
            parameters.put("input", arguments.getOptionValue("I"));
            parameters.put("content", arguments.getOptionValue("C").toLowerCase());
            parameters.put("filter", arguments.hasOption("f")
                    ? new HashSet<>(Arrays.asList(arguments.getOptionValues("f")))
                    : new HashSet<>());
            parameters.put("output", arguments.hasOption("o")
                    ? arguments.getOptionValues("o")
                    : "stdout");
        }

    }

    /**
     * Handles the sequence export task for the CLI.
     * <p>
     * This class defines the command-line options and validation logic for the sequence export task.
     * It allows users to specify input files, features, content types, samples, and other options
     * for exporting sequences from a MUSIAL storage file.
     * </p>
     */
    private static class Sequence implements Task {

        /**
         * Defines the command-line options for the sequence export task.
         * <p>
         * This method adds the following options:
         * <ul>
         *   <li>`-I` or `--input`: Path to the input `.json(.gz)` file generated with the build task.</li>
         *   <li>`-F` or `--features`: List of feature names to export data for. Non-coding features are skipped if `c` is `aa`.</li>
         *   <li>`-c` or `--content`: Specifies the content type (`nt` or `aa`, default is `nt`).</li>
         *   <li>`-s` or `--samples`: List of sample names to restrict the sequence export.</li>
         *   <li>`-m` or `--merge`: Exports sequences per sequence type instead of per sample.</li>
         *   <li>`-u` or `--unaligned`: Strips all gap characters from the sequences.</li>
         *   <li>`-k` or `--conserved`: Exports conserved sites.</li>
         *   <li>`-o` or `--output`: Path to a directory to write the output (default is the input storage directory).</li>
         * </ul>
         * </p>
         */
        private static void options() {
            options.addOption(Option.builder("I")
                    .longOpt("input")
                    .desc("Path to a .json(.gz) file generated with the build task of MUSIAL.")
                    .hasArg()
                    .required()
                    .build());
            options.addOption(Option.builder("F")
                    .longOpt("features")
                    .desc("List of feature names to export data for. Non-coding features are skipped if `content` is `aa`.")
                    .hasArgs()
                    .required()
                    .build());
            options.addOption(Option.builder("c")
                    .longOpt("content")
                    .desc("One of `nt` or `aa` (default: `nt`).")
                    .hasArg()
                    .build());
            options.addOption(Option.builder("s")
                    .longOpt("samples")
                    .desc("List of sample names to restrict the sequence export to.")
                    .hasArgs()
                    .build());
            options.addOption(Option.builder("m")
                    .longOpt("merge")
                    .desc("Export sequences per allele or proteoform instead of per sample.")
                    .build());
            options.addOption(Option.builder("x")
                    .longOpt("strip")
                    .desc("Strip all gap characters from the exported sequences.")
                    .build());
            options.addOption(Option.builder("r")
                    .longOpt("reference")
                    .desc("Include the reference sequence within the export.")
                    .build());
            options.addOption(Option.builder("k")
                    .longOpt("conserved")
                    .desc("Export conserved sites.")
                    .build());
            options.addOption(Option.builder("o")
                    .longOpt("output")
                    .desc("Path to a directory to write the output files to (default: parent of input).")
                    .hasArg()
                    .build());
        }

        /**
         * Transfers the command-line arguments for the sequence export task.
         */
        private static void transfer() {
            parameters = new HashMap<>();
            parameters.put("input", arguments.getOptionValue("I"));
            parameters.put("features", new HashSet<>(Arrays.asList(arguments.getOptionValues("F"))));
            parameters.put("content", arguments.getOptionValue("c", "nt").toLowerCase());
            parameters.put("samples", arguments.hasOption("s")
                    ? new HashSet<>(Arrays.stream(arguments.getOptionValues("s")).toList())
                    : new HashSet<>());
            parameters.put("merge", arguments.hasOption("m"));
            parameters.put("strip", arguments.hasOption("x"));
            parameters.put("conserved", arguments.hasOption("k"));
            parameters.put("reference", arguments.hasOption("r"));
            parameters.put("output", arguments.hasOption("o")
                    ? arguments.getOptionValue("o")
                    : "parent");
        }

    }

}
