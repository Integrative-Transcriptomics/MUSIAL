package cli;

import exceptions.MusialException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

import java.nio.file.Path;
import java.util.Set;

/**
 * Handles the {@code profile} task CLI parameters.
 * <p>
 * This class defines the command-line options and validation logic for the {@code profile} task. See {@link #options()} for details.
 */
public class CLIProfile implements CLI {

    /**
     * Defines the command-line options for the {@code view} task.
     * <p>
     * The options include:
     * <ul>
     *     <li><b>-I, --storage</b>: Path to the input storage file (required).</li>
     *     <li><b>-C, --content</b>: The content type to profile (VARIANTS, ALLELES, PROTEOFORMS) (required).</li>
     *     <li><b>-q, --query</b>: One or more identifiers or genomic ranges to query (optional).</li>
     *     <li><b>-x, --reduced</b>: Represent entries in a reduced format, i.e., sequence types as numbers with 0 as the reference or
     *     synonymous sequence and variants without detailed call information.</li>
     *     <li><b>-o, --output</b>: Path to the output file or a special value ("print" or "stdout") (optional).</li>
     * </ul>
     *
     * @return An {@link Options} object containing the defined CLI options.
     */
    public static Options options() {
        Options options = new Options();
        options.addOption(Option.builder("I")
                .longOpt("storage")
                .desc("Path to a .json(.gz) file generated with the build task of MUSIAL.")
                .hasArg()
                .required()
                .build());
        options.addOption(Option.builder("C")
                .longOpt("content")
                .desc("The content to view. One of VARIANTS, ALLELES, PROTEOFORMS (case-insensitive).")
                .hasArg()
                .required()
                .build());
        options.addOption(Option.builder("q")
                .longOpt("query")
                .desc("One or multiple identifiers or genomic ranges (contig:start-end) to consider.")
                .hasArgs()
                .build());
        options.addOption(Option.builder("x")
                .longOpt("reduced")
                .desc("Represent entries in a reduced format, i.e., sequence types as numbers with 0 as the reference or synonymous " +
                        "sequence " +
                        "and variants without detailed call information.")
                .hasArg(false)
                .build());
        options.addOption(Option.builder("o")
                .longOpt("output")
                .desc("Path to write the output file. If not provided, a default file will be created based on the input file (default). " +
                        "If `print` or `stdout` is specified, the output will be printed to the console.")
                .hasArg()
                .build());
        return options;
    }

    /**
     * The path to the input storage file.
     */
    public final Path input;

    /**
     * The content type to profile.
     */
    public enum Content {
        // Per sample variants.
        VARIANTS,
        // Per sample alleles of features.
        ALLELES,
        // Per sample proteoforms of features.
        PROTEOFORMS
    }

    /**
     * The content type specified by the user.
     */
    public final Content content;

    /**
     * The path to the output destination.
     */
    public final Path output;

    /**
     * The set of query filters provided by the user.
     */
    public final Set<String> query;

    /**
     * Indicates whether to represent entries in a reduced format.
     * <p>
     * If {@code true}, sequence types are represented as numbers with 0 as the reference or synonymous sequence, and variants are shown
     * without detailed call information.
     */
    public final boolean reduced;

    /**
     * Constructs a CLIView instance by parsing the provided command-line arguments.
     * <p>
     * This constructor validates and extracts the input storage file path, content type, output destination, and query filters from the
     * provided {@link CommandLine} arguments.
     *
     * @param arguments The {@link CommandLine} object containing the parsed arguments.
     * @throws MusialException If an error occurs while parsing the arguments.
     */
    public CLIProfile(CommandLine arguments) throws MusialException {
        this.input = Common.parseInputStorageFile(arguments);
        this.content = parseContent(arguments);
        String suffix = switch (this.content) {
            case VARIANTS -> "variants-profile";
            case ALLELES -> "allele-profile";
            case PROTEOFORMS -> "proteoform-profile";
        };
        this.reduced = arguments.hasOption("x");
        this.output = Common.parseOutputFile(arguments, suffix, "tsv");
        this.query = Common.parseQuery(arguments);
    }

    /**
     * Parses the content type from the command-line arguments.
     * <p>
     * This method retrieves the value of the "C" option, converts it to uppercase, and validates it against the {@link Content} enum. If
     * the value is invalid, an {@link IllegalArgumentException} is thrown.
     *
     * @param arguments The {@link CommandLine} object containing the parsed arguments.
     * @return The {@link Content} enum value corresponding to the specified content type.
     * @throws IllegalArgumentException If the content type is invalid.
     */
    private Content parseContent(CommandLine arguments) {
        String content = arguments.getOptionValue("C").toUpperCase();
        try {
            return Content.valueOf(content);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("The content to view must be one of VARIANTS, ALLELES, PROTEOFORMS (case-insensitive).");
        }
    }
}
