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
import org.apache.commons.cli.*;
import org.apache.commons.io.IOUtils;
import org.json.simple.JSONArray;

import java.io.*;
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

    public JSONArray configuration;

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
                        .desc(
                                """
                                        Path to a .json file specifying the run configurations for one or more MUSIAL modules. Please visit https://github.com/Integrative-Transcriptomics/MUSIAL for a detailed explanation on how to specify MUSIAL configuration files.
                                        Available modules are:
                                        BUILD: Generate a new variants dictionary .json file from multiple .vcf files wrt. one reference (.fasta + .gff).
                                        EXTRACT: Extract stored nucleotide or aminoacid variants of specified samples and features as a .tsv (tabular overview) or .fasta (sequences) file."""
                        )
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
                if (Objects.equals(arg, "-h") || Objects.equals(arg, "--help")) {
                    helpformatter.printHelp(
                            160,
                            "java -jar " + Musial.NAME + "-" + Musial.VERSION + ".jar",
                            "command line arguments:",
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
            // Validate input configuration.
            try (FileOutputStream schemaOut = new FileOutputStream("./MUSIAL_CONFIGURATION.schema.json")) {
                final JsonSchemaFactory factory = JsonSchemaFactory.byDefault();
                IOUtils.copy(Objects.requireNonNull(Musial.class.getResourceAsStream("/MUSIAL_CONFIGURATION.schema.json")), schemaOut);
                final JsonNode CONFIGURATION_SCHEMA = JsonLoader.fromFile(new File("./MUSIAL_CONFIGURATION.schema.json"));
                final JsonSchema SCHEMA = factory.getJsonSchema(CONFIGURATION_SCHEMA);
                final JsonNode good = JsonLoader.fromPath(arguments.getOptionValue("c"));
                ProcessingReport report = SCHEMA.validate(good);
                if (!report.isSuccess()) {
                    throw new MusialException("(Configuration Validation Failed) " + report);
                }
            } finally {
                //noinspection ResultOfMethodCallIgnored
                new File("./MUSIAL_CONFIGURATION.schema.json").delete();
            }
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
                this.configuration = gson.fromJson(bufferedReader, JSONArray.class);
            } catch (Exception e) {
                throw new MusialException(
                        "(CLI Argument Parsing Failed) " + e.getMessage());
            }
            // Set specified value for silent.
            Musial.SILENT = arguments.hasOption("s");
            // Set specified value for compress.
            Musial.COMPRESS = arguments.hasOption("k");
        } catch (ParseException e) {
            throw new MusialException("(CLI Argument Parsing Failed) " + e.getMessage());
        } catch (IOException | ProcessingException e) {
            throw new MusialException("(Configuration Validation Failed) " + e.getMessage());
        } finally {

        }
    }

}
