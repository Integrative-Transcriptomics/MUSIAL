package cli;

import exceptions.MusialException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.io.file.PathUtils;
import util.Logging;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.Set;
import java.util.function.Function;

/**
 * Handles the {@code sequence} task CLI parameters.
 * <p>
 * This class defines the command-line options and validation logic for the {@code sequence} task. See {@link #options()} for details.
 */
public class CLISequence implements CLI {

    /**
     * Defines the command-line options for the {@code sequence} task.
     * <p>
     * The options include:
     * <ul>
     *     <li><b>-I, --storage</b>: The path to the input storage file (required).</li>
     *     <li><b>-c, --content</b>: Whether to generate NUCLEOTIDE or AMINOACID sequences (default: NUCLEOTIDE).</li>
     *     <li><b>-l, --locations</b>: One or more feature identifiers or genomic ranges to generate sequence data for (optional).</li>
     *     <li><b>-s, --samples</b>: One or more sample identifiers to retrieve sequences for (optional).</li>
     *     <li><b>-m, --merge</b>: Whether to merge identical sequences (optional, default: false).</li>
     *     <li><b>-f, --split</b>: How to split output files (optional, default: FEATURE).</li>
     *     <li><b>-a, --align</b>: Whether to align sequences (optional, default: false).</li>
     *     <li><b>-v, --variable</b>: Whether to consider only variable sites (optional, default: false).</li>
     *     <li><b>-o, --output</b>: Path to write the output. By default, files will be created in the input's directory (optional).</li>
     * </ul>
     *
     * @return An {@link Options} object containing the defined command-line options.
     */
    public static Options options() {
        Options options = new Options();
        options.addOption(Option.builder("I")
                .longOpt("storage")
                .desc("Path to a .json(.gz) file generated with the build task of MUSIAL.")
                .hasArg()
                .required()
                .build());
        options.addOption(Option.builder("l")
                .longOpt("locations")
                .desc("One or multiple feature identifiers or genomic ranges (contig:start-end) to generate sequence data of. If none are" +
                        " provided, all features or full contig ranges will be considered.")
                .hasArgs()
                .build());
        options.addOption(Option.builder("c")
                .longOpt("content")
                .desc("Whether to generate NUCLEOTIDE or AMINOACID sequences (optional, case-insensitive, default: NUCLEOTIDE).")
                .hasArg()
                .build());
        options.addOption(Option.builder("s")
                .longOpt("samples")
                .desc("One or multiple sample identifiers to retrieve sequences for (optional).")
                .hasArgs()
                .build());
        options.addOption(Option.builder("m")
                .longOpt("merge")
                .desc("Whether to merge identical sequences (optional, default: false).")
                .build());
        options.addOption(Option.builder("f")
                .longOpt("split")
                .desc("Whether to split output files by FEATURE, SAMPLE, BOTH, or NONE (optional, case-insensitive, default: FEATURE).")
                .hasArg()
                .build());
        options.addOption(Option.builder("a")
                .longOpt("align")
                .desc("Whether to align sequences (optional, default: false).")
                .build());
        options.addOption(Option.builder("v")
                .longOpt("variable")
                .desc("Whether to only consider variable positions (optional, default: false).")
                .build());
        options.addOption(Option.builder("o")
                .longOpt("output")
                .desc("Path to write the output. If not provided, the directory of the input storage is used. If a directory is provided," +
                        " files are created there. If a file is provided, its parent directory is used.")
                .hasArg()
                .build());
        return options;
    }

    /**
     * The path to the input storage file.
     */
    public final Path input;

    /**
     * The content type to generate sequences of.
     */
    public enum Content {
        /**
         * Nucleotide/genomic sequences.
         */
        NUCLEOTIDE,
        /**
         * Protein sequences.
         */
        AMINOACID
    }

    /**
     * The mode to split output files by.
     */
    public enum Split {
        /**
         * Split output files by feature/genomic range.
         */
        FEATURE,
        /**
         * Split output files by sample (or sequence type) identifier.
         */
        SAMPLE,
        /**
         * Do not split output files.
         */
        NONE,
        /**
         * Split output files by both feature/genomic range and sample (or sequence type) identifier.
         */
        BOTH
    }

    /**
     * The content type specified by the user.
     */
    public final Content content;

    /**
     * The mode to split output files by.
     */
    public final Split split;

    /**
     * Generator for output paths per specified locus.
     */
    public final Function<String, String> outputGenerator;

    /**
     * The set of loci (features or genomic ranges) provided by the user.
     */
    public final Set<String> loci;

    /**
     * The set of sample identifiers provided by the user.
     */
    public final Set<String> samples;

    /**
     * Indicates whether to merge identical sequences.
     */
    public final boolean merge;

    /**
     * Indicates whether to align sequences.
     */
    public final boolean align;

    /**
     * Indicates whether only variable sites should be considered.
     */
    public final boolean variable;

    /**
     * Constructs a new {@code CLISequence} instance and initializes its fields based on the provided command-line arguments.
     * <p>
     * This constructor parses and validates the command-line arguments to initialize the fields required for the {@code sequence} task. It
     * sets the input storage file, content type, loci, samples, merge and align options, and the output path generator.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @throws MusialException If there is an error parsing the output path or validating the input storage file.
     */
    public CLISequence(CommandLine arguments) throws MusialException {
        // Parse and validate the input storage file path.
        this.input = Common.parseInputStorageFile(arguments);

        // Parse the content type (NUCLEOTIDE or AMINOACID) from the arguments.
        this.content = parseContent(arguments);

        // Parse the mode to split output files by from the arguments.
        this.split = parseSplit(arguments);

        // Parse the loci (features or genomic ranges) from the arguments.
        this.loci = parseLoci(arguments);

        // Parse the sample identifiers from the arguments.
        this.samples = parseSamples(arguments);

        // Check if the merge option is enabled in the arguments.
        this.merge = arguments.hasOption("m");

        // Validate that merging is not used with incompatible split modes.
        if (this.merge && (this.split.equals(Split.SAMPLE))) {
            throw new MusialException("Merging sequences is not compatible with splitting output files by SAMPLE.");
        }

        // Check if the align option is enabled in the arguments.
        this.align = arguments.hasOption("a");

        // Check if the variable option is enabled in the arguments.
        this.variable = arguments.hasOption("v");

        // Generate the output path function based on the arguments.
        this.outputGenerator = parseOutput(arguments);
    }

    /**
     * Parses the content type specified in the command-line arguments.
     * <p>
     * This method checks if the `-c` or `--content` option is provided in the command-line arguments. If the option is not provided, it
     * defaults to {@link Content#NUCLEOTIDE}. If the option is provided, it attempts to parse the value as a valid {@link Content} enum.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @return The {@link Content} type specified by the user, or {@link Content#NUCLEOTIDE} if not specified.
     * @throws IllegalArgumentException If the provided content type is not one of the valid {@link Content} values.
     */
    private Content parseContent(CommandLine arguments) {
        if (!arguments.hasOption("c")) {
            return Content.NUCLEOTIDE;
        } else {
            String content = arguments.getOptionValue("c").toUpperCase();
            try {
                return Content.valueOf(content);
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("The content to generate sequences must be one of NUCLEOTIDE or AMINOACID " +
                        "(case-insensitive).");
            }
        }
    }

    /**
     * Parses the mode to split output files by from the command-line arguments.
     * <p>
     * This method checks if the `-f` or `--split` option is provided in the command-line arguments. If the option is not provided, it
     * defaults to {@link Split#FEATURE}. If the option is provided, it attempts to parse the value as a valid {@link Split} enum.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @return The {@link Split} mode specified by the user, or {@link Split#FEATURE} if not specified.
     * @throws IllegalArgumentException If the provided split mode is not one of the valid {@link Split} values.
     */
    private Split parseSplit(CommandLine arguments) {
        if (!arguments.hasOption("f")) {
            return Split.FEATURE;
        } else {
            String split = arguments.getOptionValue("f").toUpperCase();
            try {
                return Split.valueOf(split);
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("The mode to split output files by must be one of FEATURE, SAMPLE, NONE, or BOTH " +
                        "(case-insensitive).");
            }
        }
    }

    /**
     * Parses the loci (features or genomic ranges) provided in the command-line arguments.
     * <p>
     * This method retrieves the values associated with the `-l` or `--locations` option from the command-line arguments. If no loci are
     * provided, an empty set is returned. The values are returned as a {@link Set} for further processing. The values are not validated at
     * this stage.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @return A {@link Set} of loci (features or genomic ranges) specified by the user.
     */
    private Set<String> parseLoci(CommandLine arguments) {
        String[] loci;
        if (arguments.hasOption("l")) {
            loci = arguments.getOptionValues("l");
        } else {
            loci = new String[0];
        }
        return Set.of(loci);
    }

    /**
     * Parses the sample identifiers provided in the command-line arguments.
     * <p>
     * This method retrieves the values associated with the `-s` or `--samples` option from the command-line arguments. If no sample
     * identifiers are provided, it returns an empty set. Otherwise, it converts the array of sample identifiers into a set. The values are
     * not validated at this stage.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @return A {@link Set} of sample identifiers. If no sample identifiers are provided, an empty set is returned.
     */
    private Set<String> parseSamples(CommandLine arguments) {
        String[] samples = arguments.getOptionValues("s");
        if (samples == null || samples.length == 0) {
            return Collections.emptySet();
        }
        return Set.of(samples);
    }

    /**
     * Generates a function to construct output file names.
     * <p>
     * This method constructs a suffix for the output file name based on the content type, alignment, and merge options. It validates and
     * creates the base output path, which can either be a directory or a file.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @return A {@link Function} that takes a locus identifier as input and returns the corresponding output file path.
     * @throws MusialException If the output path specified in the arguments is invalid or cannot be created.
     */
    private Function<String, String> parseOutput(CommandLine arguments) throws MusialException {
        // Construct the suffix for the output file name based on content type, alignment, and merge options.
        String suffix = (this.align ? "-aligned" : "")
                + (this.merge ? "-merged" : "")
                + (this.variable ? "-variants" : "")
                + (this.content == Content.NUCLEOTIDE ? ".fna" : ".faa");

        try {
            // Determine the base output path. Use the specified output path if provided, otherwise use the parent directory of the input
            // or specified file.
            Path basePath = arguments.hasOption("o")
                    ? Path.of(arguments.getOptionValue("o"))
                    : Path.of(arguments.getOptionValue("I"));

            // Get the parent directory if the base path is a file.
            if (PathUtils.isRegularFile(basePath)) {
                basePath = basePath.getParent();
            }

            // Ensure the parent directories for the base path exist.
            Files.createDirectories(basePath);

            Logging.logConfig("`output` will be generated at %s.".formatted(basePath));
            Path finalBasePath = basePath;
            return s -> finalBasePath.resolve("%s%s".formatted(s, suffix)).toAbsolutePath().toString();
        } catch (Exception e) {
            // Throw an exception if the output path is invalid or cannot be created.
            throw new MusialException("Failed to validate path %s specified for `output`.".formatted(arguments.getOptionValue("o")));
        }
    }
}
