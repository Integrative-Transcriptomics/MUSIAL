package cli;

import components.IO;
import components.Logging;
import components.Validation;
import datastructure.VariantsDictionary;
import exceptions.MusialIOException;
import main.Musial;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Objects;
import java.util.function.Consumer;

/**
 * Parses command line interface arguments for the `MUSIAL inferSequences` module.
 * <p>
 * Used to parse and validate passed command line options or print help information to the user. Parsed arguments are
 * stored to be accessible for other components of the tool.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public final class CLIParametersInferSequences implements CLIParameters {

    /**
     * The unprocessed command line content passed by the user.
     */
    public String[] ARGUMENTS;
    /**
     * Whether to, in addition to sequence output, a multiple sequence alignment shall be built.
     */
    public boolean outputMSA = false;
    /**
     * Whether sequence information for optionally allocated proteins of specified features and samples shall be
     * extracted.
     * <p>
     * If `true`, proteoforms with any sample specified in {@link CLIParametersInferSequences#samples} will be
     * considered.
     */
    public boolean includePS = false;
    /**
     * {@link ArrayList} of {@link String}s yielding information about for which samples sequences shall be extracted.
     * <p>
     * If empty, sequences for all samples will be extracted.
     */
    public final ArrayList<String> samples = new ArrayList<>();
    /**
     * {@link ArrayList} of {@link String}s yielding information about for which features sequences shall be extracted.
     * <p>
     * If empty, sequences for all features will be extracted.
     */
    public final ArrayList<String> features = new ArrayList<>();
    /**
     * {@link File} object specifying the output directory.
     */
    public File outputDirectory;
    /**
     * {@link VariantsDictionary} to infer sequences from.
     */
    public VariantsDictionary inputVDict;

    /**
     * Constructor of the {@link CLIParametersInferSequences} class.
     * <p>
     * Used to parse and validate command line interface input. Parsed arguments are stored - and accessible by other
     * components - via class properties.
     *
     * @param args {@link String} {@link Array} containing the command line arguments.
     * @throws MusialIOException If any error occurs during input file validation.
     */
    public CLIParametersInferSequences(String[] args)
            throws MusialIOException {
        // Store original command line arguments.
        this.ARGUMENTS = args;
        // Initialize `Option` object with all parameters.
        Options options = new Options();
        options.addOption(Option.builder("I")
                .longOpt("inputVDict")
                .desc("Path to .vdict.json file used as input variants dictionary from which sequence information shall be inferred.")
                .hasArg()
                .required()
                .build()
        );
        options.addOption(Option.builder("O")
                .longOpt("outputDirectory")
                .desc("Path to output directory to store inferred sequences.")
                .hasArg()
                .required()
                .build()
        );
        options.addOption(Option.builder("MSA")
                .longOpt("outputMSA")
                .desc("If specified, sequences will be additionally written as multiple sequence alignment.")
                .build()
        );
        options.addOption(Option.builder("PS")
                .longOpt("includePS")
                .desc("If specified, protein sequences, if protein information is allocated to specified features, will be output.")
                .build()
        );
        options.addOption(Option.builder("S")
                .longOpt("samples")
                .desc("Comma separated internal sample names to match in the specified input .vdict.json. If none are specified all samples will be used.")
                .hasArg()
                .build()
        );
        options.addOption(Option.builder("F")
                .longOpt("features")
                .desc("Comma separated internal feature names to match in the sepcified input .vdict.json. If none are specified all samples will be used.")
                .hasArg()
                .build());
        // Instantiate a formatter for the help message and a default command line parser.
        HelpFormatter helpformatter = new HelpFormatter();
        CommandLineParser parser = new DefaultParser();
        Consumer<Options> checkHelp = (Options checkOptions) -> {
            for (String arg : args) {
                if (Objects.equals(arg, "-h") || Objects.equals(arg, "--help") || Objects.equals(arg, "h") || Objects.equals(arg, "help")) {
                    printHelp(options, helpformatter);
                    System.exit(0);
                }
            }
        };
        // If `-h`/`--help` was specified with other options apply the same behaviour.
        checkHelp.accept(options);
        // Else, command line interface arguments are parsed.
        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
            // Parse specified input variants dictionary.
            File inputVDictFile = new File(cmd.getOptionValue("I"));
            try {
                this.inputVDict = IO.readVariantsDictionary(inputVDictFile);
            } catch (IOException e) {
                throw new MusialIOException("Failed to parse specified variants dictionary at " + cmd.getOptionValue("I") + ".");
            }
            // Parse specified output directory.
            this.outputDirectory = new File(cmd.getOptionValue("O"));
            if (!Validation.isDirectory(this.outputDirectory)) {
                throw new MusialIOException("The specified output directory can not be accessed or is no valid directory.");
            }
            // Parse optional MSA and PS parameters.
            this.outputMSA = cmd.hasOption("MSA");
            this.includePS = cmd.hasOption("PS");
            // Parse and validate samples.
            if (cmd.hasOption("S")) {
                for (String s : cmd.getOptionValue("S").split(",")) {
                    if (this.inputVDict.samples.containsKey(s)) {
                        this.samples.add(s);
                    } else {
                        Logging.logWarning("Sample " + s + " does not exist in the specified variants dictionary and will be skipped.");
                    }
                }
            }
            if (this.samples.isEmpty()) {
                this.samples.addAll(this.inputVDict.samples.keySet());
            }
            // Parse and validate features.
            if (cmd.hasOption("F")) {
                for (String f : cmd.getOptionValue("F").split(",")) {
                    if (this.inputVDict.features.containsKey(f)) {
                        this.features.add(f);
                    } else {
                        Logging.logWarning("Feature " + f + " does not exist in the specified variants dictionary and will be skipped.");
                    }
                }
            }
            if (this.features.isEmpty()) {
                this.features.addAll(this.inputVDict.features.keySet());
            }
        } catch (ParseException e) {
            Logging.logError("An error occurred during command line interface parameter parsing: " + e.getMessage());
        }
    }

    /**
     * Print help information for command line interface arguments.
     *
     * @param options       The specified command line interface {@link Options}.
     * @param helpFormatter A {@link HelpFormatter} for displaying help information.
     */
    public void printHelp(Options options, HelpFormatter helpFormatter) {
        helpFormatter.printHelp(
                100,
                "java -jar MUSIAL-" + Musial.VERSION + ".jar inferSequences",
                IO.LINE_SEPARATOR + "Extract nucleotide and/or amino-acid sequences as plain sequences or multiple sequence alignment",
                options,
                "",
                true
        );
    }

}