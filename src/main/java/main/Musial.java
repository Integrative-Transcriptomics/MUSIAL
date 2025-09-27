package main;

import cli.*;
import exceptions.MusialException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.io.FileUtils;
import task.*;
import util.Logging;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Properties;
import java.util.logging.Level;

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
     * Private constructor to prevent instantiation of this utility class.
     */
    private Musial() {
    }

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
    public static MusialTask task;

    /**
     * File extension used for generated MUSIAL storage files.
     * <p>
     * This constant specifies the default file extension for storage files created or used by the MUSIAL application. It is marked as
     * `transient` to indicate that it should not be serialized as part of the class state.
     * <p>
     * The extension should either be `.json` or `.json.gz` depending on the compression method used. For production, use `.json.gz` for
     * compressed files.
     */
    public static final String OUTPUT_EXTENSION = ".json.gz";

    /**
     * Start time of the program.
     */
    public static long startTime;

    /**
     * Temporary directory for intermediate files.
     */
    public static File tempDir;

    /**
     * <b>Only relevant for development! For deployment this should be set to {@link Level#CONFIG}.</b>
     * <p>
     * The verbosity level used for logging. See {@link Level} for details.
     */
    private static final Level LOG_VERBOSITY = Level.CONFIG;

    /**
     * The main entry point of the MUSIAL application.
     * <p>
     * This method initializes the program, determines the task to execute based on the provided arguments, and executes the corresponding
     * functionality. It handles errors gracefully and logs relevant information.
     *
     * @param args Command-line arguments specifying the task and its parameters.
     *             <ul>
     *                 <li><b>args[0]</b>: The task to execute (e.g., BUILD, UPDATE, VIEW, SEQUENCE).</li>
     *                 <li>Additional arguments are parsed by implementations of {@link CLI}.</li>
     *             </ul>
     */
    public static void main(String[] args) {
        // Execution state of the program; 0 = OK, 1 = Error (internal), 2 = Error (unexpected).
        int status = 0;
        try {
            // Initialize the logging system.
            Logging.init(LOG_VERBOSITY);

            // Record the start time of the program.
            if (LOG_VERBOSITY.intValue() <= Level.FINE.intValue()) startTime = System.currentTimeMillis();

            // Load metadata such as software id, version, and contact information.
            loadMetadata();

            // Create temporary directory for intermediate files, if not present.
            tempDir = Files.createTempDirectory(Path.of(System.getProperty("user.dir")), Musial.name.toLowerCase()).toFile();

            // Check if any arguments were provided; if not, display usage information and exit.
            if (args.length == 0) {
                System.out.printf("No arguments were specified. Call `java -jar %s-%s.jar [-h|--help]` for more information.%n",
                        Musial.name, Musial.version);
                System.exit(0);
            }

            // Parse the first argument to determine the task to execute.
            try {
                task = MusialTask.valueOf(args[0].toUpperCase());
            } catch (IllegalArgumentException e) {
                // If the task is invalid, set it to UNDEFINED.
                task = MusialTask.UNDEFINED;
            }

            // Fill options dependent on task to execute.
            Options options = switch (Musial.task) {
                // NOTE: New tasks need to be added here.
                case BUILD -> CLIBuild.options();
                case EXPAND -> CLIExpand.options();
                case VIEW -> CLIView.options();
                case PROFILE -> CLIProfile.options();
                case SEQUENCE -> CLISequence.options();
                default -> new Options();
            };

            // Print help and exit if -h or --help is provided.
            if (Arrays.stream(args).sequential().anyMatch(a -> a.equals("-h") || a.equals("--help"))) {
                exitHelp(options, args);
            }

            // Parse arguments from command line.
            CommandLine arguments = new DefaultParser().parse(options, args);

            // Execute the task based on the parsed value.
            switch (task) {
                case BUILD -> {
                    Logging.logInfo("Execute task \033[1mbuild\033[0m");
                    CLIBuild cli = new CLIBuild(arguments);
                    ExecutorBuild executor = new ExecutorBuild(cli);
                    executor.run();
                }
                case EXPAND -> {
                    Logging.logInfo("Execute task \033[1mexpand\033[0m");
                    CLIExpand cli = new CLIExpand(arguments);
                    ExecutorExpand executor = new ExecutorExpand(cli);
                    executor.run();
                }
                case VIEW -> {
                    Logging.logInfo("Execute task \033[1mview\033[0m");
                    CLIView cli = new CLIView(arguments);
                    ExecutorView executor = new ExecutorView(cli);
                    executor.run();
                }
                case PROFILE -> {
                    Logging.logInfo("Execute task \033[1mprofile\033[0m");
                    CLIProfile cli = new CLIProfile(arguments);
                    ExecutorProfile executor = new ExecutorProfile(cli);
                    executor.run();
                }
                case SEQUENCE -> {
                    Logging.logInfo("Execute task \033[1msequence\033[0m");
                    CLISequence cli = new CLISequence(arguments);
                    ExecutorSequence executor = new ExecutorSequence(cli);
                    executor.run();
                }
                // Exit the program, if the task is undefined.
                default -> exitNotRecognized();
            }
        } catch (Exception e) {
            // Log the error message and stack trace, then exit with an error code.
            if (e.getClass().equals(MusialException.class)) {
                Logging.logExit("An internal error has occurred: %s".formatted(e.getMessage()));
                status = 1;
            } else {
                Logging.logExit("An unexpected error has occurred: %s".formatted(e.getMessage()));
                status = 2;
            }
            if (LOG_VERBOSITY.intValue() <= Level.FINE.intValue()) e.printStackTrace();
        } finally {
            // Log the total execution time if the verbosity level is set to FINE or lower.
            if (LOG_VERBOSITY.intValue() <= Level.FINE.intValue()) {
                long endTime = System.currentTimeMillis();
                long duration = endTime - startTime;
                Logging.logDebug("Total execution time: %.1f s (%.1f min)".formatted(duration / 1000.0, duration / 60000.0));
            }

            // Clean up temporary directory.
            try {
                FileUtils.deleteDirectory(Musial.tempDir);
            } catch (IOException e) {
                Logging.logSevere(e.getMessage());
            }

            System.exit(status);
        }
    }

    /**
     * Loads metadata, such as the software title and version from `/src/main/resources/info.properties` and prints the information to
     * stdout.
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
     * Exits the program when an unrecognized task is provided.
     * <p>
     * This method prints an error message to the console indicating that the specified task is not recognized. It provides instructions on
     * how to access help information for valid tasks. After displaying the message, the program terminates with an exit code of 0.
     */
    private static void exitNotRecognized() {
        System.out.printf("Task \033[1;31m%s\033[0m not recognized. Call `java -jar %s-%s.jar [-h|--help]` for more information.%n",
                Musial.task, Musial.name, Musial.version);
        System.exit(0);
    }

    /**
     * Displays the help message for the MUSIAL application and exits the program.
     * <p>
     * This method generates and prints a detailed help message based on the provided options and arguments. The help message includes a
     * description of the application, available tasks, and specific command-line arguments for the selected task. If the task is
     * unrecognized, the help message indicates this and provides general usage instructions. After displaying the help message, the program
     * terminates with an exit code of 0.
     *
     * @param options The {@link Options} object containing the command-line options for the specified task.
     * @param args    The command-line arguments provided by the user.
     */
    private static void exitHelp(Options options, String[] args) {
        // NOTE: New tasks need to be added here.
        String helpText = """

                MUSIAL aggregates and analyzes variant calls from multiple samples of a prokaryotic species and provides an interface to generate comprehensive statistics and alignments at the genome, gene and protein level.

                Available tasks are:
                \033[47m\033[1;30m build    \033[0m : Build a local database file (storage) in JSON format from variant calls; the mandatory input for other tasks.
                \033[47m\033[1;30m expand   \033[0m : Expand an existing storage file from variant call files and/or meta data.
                \033[47m\033[1;30m view     \033[0m : View the content (features, samples or variants; and their attributes) of a MUSIAL storage file.
                \033[47m\033[1;30m profile  \033[0m : Profile samples with respect to variants, alleles, or proteoforms.
                \033[47m\033[1;30m sequence \033[0m : Generate and write sequence data.

                Call `java -jar %s-%s.jar <task> [-h|--help]` for more information.
                """.formatted(Musial.name, Musial.version);

        // If a task was specified, but is not recognized, adjust the help text.
        if (Musial.task.equals(MusialTask.UNDEFINED)) {
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
        new HelpFormatter().printHelp(
                150,
                "java -jar %s-%s.jar %s".formatted(Musial.name, Musial.version,
                        Musial.task.equals(MusialTask.UNDEFINED) ? "<task>" : Musial.task.toString().toLowerCase()),
                helpText,
                options,
                "",
                true
        );

        // Exit the program after displaying the help message.
        System.exit(0);
    }

}