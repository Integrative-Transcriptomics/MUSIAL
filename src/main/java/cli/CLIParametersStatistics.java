package cli;

import components.IO;
import components.Logging;
import components.Validation;
import datastructure.VariantsDictionary;
import exceptions.MusialException;
import main.Musial;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Objects;
import java.util.function.Consumer;

/**
 * Parses command line interface arguments for the `MUSIAL statistics` module.
 * <p>
 * Used to parse and validate passed command line options or print help information to the user. Parsed arguments are
 * stored to be accessible for other components of the tool.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public final class CLIParametersStatistics implements CLIParameters {

    /**
     * The unprocessed command line content passed by the user.
     */
    public String[] ARGUMENTS;
    /**
     * {@link ArrayList} of {@link String}s yielding information about for which samples variants should be annotated and
     * used for statistics computation.
     * <p>
     * If empty, all samples are considered.
     */
    public final ArrayList<String> samples = new ArrayList<>();
    /**
     * {@link ArrayList} of {@link String}s yielding information about for which features variants should be annotated and
     * used for statistics computation.
     * <p>
     * If empty, all features are considered.
     */
    public final ArrayList<String> features = new ArrayList<>();
    /**
     * {@link File} object specifying the reference genome sequence .fasta file.
     */
    public File referenceFASTA;
    /**
     * {@link File} object specifying the reference genome annotation .gff file.
     */
    public File referenceGFF;
    /**
     * {@link File} object specifying the output directory.
     */
    public File outputDirectory;
    /**
     * Input {@link VariantsDictionary} to compute statistics for.
     */
    public VariantsDictionary inputVDict;

    /**
     * Constructor of the {@link CLIParametersStatistics} class.
     * <p>
     * Used to parse and validate command line interface input. Parsed arguments are stored - and accessible by other
     * components - via class properties.
     *
     * @param args {@link String} {@link Array} containing the command line arguments.
     * @throws MusialException If IO file validation fails.
     */
    public CLIParametersStatistics(String[] args)
            throws MusialException {
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
        options.addOption(Option.builder("R")
                .longOpt("referenceFasta")
                .desc("Path to .fasta file containing the reference genome; Has to be the same file as the one used for constructing the input variants dictionary.")
                .hasArg()
                .required()
                .build()
        );
        options.addOption(Option.builder("A")
                .longOpt("referenceAnnotation")
                .desc("Path to .gff file containing the reference genome annotation; Has to be the same file as the one used for constructing the input variants dictionary.")
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
                throw new MusialException("(CLI Parameter Parsing/Initialization) Failed to parse specified variants dictionary at " + cmd.getOptionValue("I") + ".");
            }
            // Parse reference genome (.fasta)
            File r = new File(cmd.getOptionValue("R"));
            if (Validation.isFile(r)) {
                this.referenceFASTA = r;
            } else {
                throw new MusialException(
                        "(CLI Parameter Parsing/Initialization) Invalid `referenceFASTA` " + cmd.getOptionValue("R") + "; failed to read file.");
            }
            // Parse reference annotation (.gff/.gff3)
            File a = new File(cmd.getOptionValue("A"));
            if (Validation.isFile(a)) {
                this.referenceGFF = a;
            } else {
                throw new MusialException(
                        "(CLI Parameter Parsing/Initialization) Invalid `referenceGFF` " + cmd.getOptionValue("A") + "; failed to read file.");
            }
            // Parse specified output directory.
            this.outputDirectory = new File(cmd.getOptionValue("O"));
            if (!Validation.isDirectory(this.outputDirectory)) {
                throw new MusialException("(CLI Parameter Parsing/Initialization) The specified output directory can not be accessed or is no valid directory.");
            }
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
                "java -jar MUSIAL-" + Musial.VERSION + ".jar statistics",
                IO.LINE_SEPARATOR + "Annotate variants with SnpEff and infer various statistics",
                options,
                "",
                true
        );
    }

}