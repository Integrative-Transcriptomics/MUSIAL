package cli;

import exceptions.MusialException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.io.FileUtils;
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
     *     <li><b>-I, --storage</b>: Specifies the path to the input storage file (required).</li>
     *     <li><b>-C, --content</b>: Specifies whether to generate NUCLEOTIDE or AMINOACID sequences (required).</li>
     *     <li><b>-L, --locations</b>: Specifies one or more feature identifiers or genomic ranges to generate sequence data for
     *     (required).</li>
     *     <li><b>-s, --samples</b>: Specifies one or more sample identifiers to retrieve sequences for (optional).</li>
     *     <li><b>-m, --merge</b>: Indicates whether to merge identical sequences (optional, default: false).</li>
     *     <li><b>-a, --align</b>: Indicates whether to align sequences (optional, default: false).</li>
     *     <li><b>-v, --variable</b>: Indicates whether to consider only variable sites (optional, default: false).</li>
     *     <li><b>-o, --output</b>: Specifies the path to write the output. If not provided, a default file will be created (optional).</li>
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
        options.addOption(Option.builder("L")
                .longOpt("locations")
                .desc("One or multiple feature identifiers or genomic ranges (contig:start-end) to generate sequence data of.")
                .hasArgs()
                .required()
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
                .desc("Path to write the output. If not provided, a file with default file name will be created next to the input file " +
                        " for each specified location (default). If a directory is provided, a respective file is created there. If a " +
                        "file is provided, all sequences will be written to the same file.")
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
     * The content type specified by the user.
     */
    public final Content content;

    /**
     * Generator for output paths per specified locus.
     */
    public final Function<String, String> outputGenerator;

    /**
     * Indicates whether to append to existing output files.
     */
    public boolean append;

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
     * </p>
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @throws MusialException If there is an error parsing the output path or validating the input storage file.
     */
    public CLISequence(CommandLine arguments) throws MusialException {
        // Parse and validate the input storage file path.
        this.input = Common.parseInputStorageFile(arguments);

        // Parse the content type (NUCLEOTIDE or AMINOACID) from the arguments.
        this.content = parseContent(arguments);

        // Parse the loci (features or genomic ranges) from the arguments.
        this.loci = parseLoci(arguments);

        // Parse the sample identifiers from the arguments.
        this.samples = parseSamples(arguments);

        // Check if the merge option is enabled in the arguments.
        this.merge = arguments.hasOption("m");

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
     * Parses the loci (features or genomic ranges) provided in the command-line arguments.
     * <p>
     * This method retrieves the values associated with the `-L` or `--locations` option from the command-line arguments. If no loci are
     * provided, it throws an {@link IllegalArgumentException}. The values are returned as a {@link Set} for further processing. The values
     * are not validated at this stage.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @return A {@link Set} of loci (features or genomic ranges) specified by the user.
     * @throws IllegalArgumentException If no loci are provided in the command-line arguments.
     * @throws MusialException          If an error occurs while parsing the loci.
     */
    private Set<String> parseLoci(CommandLine arguments) throws MusialException {
        String[] loci = arguments.getOptionValues("L");
        if (loci == null || loci.length == 0) {
            throw new MusialException("At least one location (feature ID or genomic range) must be specified.");
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
     * Generates a function to determine the output file path for each locus based on the command-line arguments.
     * <p>
     * This method constructs a suffix for the output file name based on the content type, alignment, and merge options. It validates and
     * creates the base output path, which can either be a directory or a file. If the base path is a directory, the output file path is
     * generated per locus. If the base path is a file, all sequences are written to the same file.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @return A {@link Function} that takes a locus identifier as input and returns the corresponding output file path.
     * @throws MusialException If the output path specified in the arguments is invalid or cannot be created.
     */
    private Function<String, String> parseOutput(CommandLine arguments) throws MusialException {
        // Construct the suffix for the output file name based on content, alignment, and merge options.
        String suffix =
                (this.content == Content.NUCLEOTIDE ? "n" : "a")
                        + (this.align ? "a" : "s")
                        + (this.merge ? "m" : "")
                        + (this.variable ? "v" : "");

        try {
            // Determine the base output path. Use the specified output path if provided, otherwise use the parent directory of the input
            // file.
            Path basePath = arguments.hasOption("o")
                    ? Path.of(arguments.getOptionValue("o"))
                    : Path.of(arguments.getOptionValue("I")).getParent();

            // Ensure the parent directories for the base path exist.
            FileUtils.createParentDirectories(basePath.toFile());

            // If the base path is a directory, generate output paths per locus.
            if (Files.isDirectory(basePath)) {
                Logging.logConfig("`output` will be generated per locus at %s.".formatted(basePath));
                this.append = false;
                return s -> basePath.resolve("musial-%s-%s-%s.fasta".formatted(s, suffix, Logging.getDate())).toAbsolutePath().toString();
            } else {
                // If the base path is a file, use it as the output path for all sequences.
                Logging.logConfig("`output` set to file %s.".formatted(basePath));
                this.append = true;
                return s -> basePath.toAbsolutePath().toString();
            }
        } catch (Exception e) {
            // Throw an exception if the output path is invalid or cannot be created.
            throw new MusialException("Failed to validate path %s specified for `output`.".formatted(arguments.getOptionValue("o")));
        }
    }
}
