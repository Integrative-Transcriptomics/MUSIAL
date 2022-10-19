package cli;

import com.google.gson.Gson;
import components.IO;
import components.Logging;
import exceptions.MusialException;
import main.Musial;
import org.apache.commons.cli.*;
import org.json.simple.JSONObject;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.lang.reflect.Array;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Objects;
import java.util.function.Consumer;

/**
 * Parses command line interface arguments.
 * <p>
 * Used to parse and validate passed command line options or print help information to the user. Parsed arguments are
 * stored to be accessible for other components of the tool.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class CLIParser {

    public JSONObject configuration;

    public CommandLine arguments;

    /**
     * Constructor of the {@link CLIParser} class.
     *
     * @param args {@link String} {@link Array} of passed command line arguments.
     */
    public CLIParser(String[] args) throws MusialException {
        // Initialize `Option` object with all parameters.
        Options options = new Options();
        // Add configuration option that is used to specify configuration JSON file.
        options.addOption(
                Option.builder("c")
                        .longOpt("configuration")
                        .desc("Path to a .JSON file yielding the parameter sets to run MUSIAL modules. Please visit " + Logging.colorParameter("https://github.com/Integrative-Transcriptomics/MUSIAL") + " for more information about how to specify MUSIAL configuration files.")
                        .hasArg()
                        .required()
                        .build());
        // Add option to silence the software.
        options.addOption(
                Option.builder("s")
                        .longOpt("silent")
                        .desc("If set, no console output will be generated.")
                        .build()
        );
        // Add option to compress output.
        options.addOption(
                Option.builder("k")
                        .longOpt("compress")
                        .desc("If set, any output files will be compressed.")
                        .build()
        );
        // Instantiate a formatter for the help message and a default command line parser.
        HelpFormatter helpformatter = new HelpFormatter();
        CommandLineParser parser = new DefaultParser();
        Consumer<Options> checkHelp = (Options checkOptions) -> {
            for (String arg : args) {
                if (Objects.equals(arg, "-h")
                        || Objects.equals(arg, "--help")
                        || Objects.equals(arg, "h")
                        || Objects.equals(arg, "help")) {
                    Musial.printInfo();
                    helpformatter.printHelp(
                            160,
                            "java -jar " + Musial.NAME + "-" + Musial.VERSION + ".jar -configuration <FILE_PATH> [-silent|-compress]",
                            IO.LINE_SEPARATOR + "TODO",
                            options,
                            "",
                            true
                    );
                    System.exit(0);
                }
            }
        };
        // If `-h`/`--help` was specified with other options apply the same behaviour.
        checkHelp.accept(options);
        // Else, command line interface arguments are parsed.
        try {
            arguments = parser.parse(options, args);
            // Parse specified configuration `JSON` object.
            try (
                    BufferedReader bufferedReader =
                            new BufferedReader(
                                    new InputStreamReader(
                                            Files.newInputStream(Path.of(arguments.getOptionValue("c"))), StandardCharsets.UTF_8
                                    )
                            )
            ) {
                Gson gson = new Gson();
                this.configuration = gson.fromJson(bufferedReader, JSONObject.class);
            } catch (Exception e) {
                throw new MusialException(
                        "(CLI Argument Parsing Failed) " + e.getMessage());
            }
            // Set specified value for silent.
            Musial.SILENT = options.hasOption("s");
            // Set specified value for compress.
            Musial.COMPRESS = options.hasOption("k");
        } catch (ParseException e) {
            throw new MusialException("(CLI Argument Parsing Failed) " + e.getMessage());
        }
    }

}
