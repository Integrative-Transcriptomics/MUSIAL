package cli;

import exceptions.MusialCLException;
import java.io.File;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Parses command line interface arguments for Musial generateMSA module.
 * <p>
 * Used to parse and validate passed command line options or print help information to the user. Parsed arguments are
 * stored to be accessible for other components of the tool.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public final class CLIParametersGenerateMSA implements CLIParameters {

  /**
   * The unprocessed command line content passed by the user.
   */
  final String[] arguments;
  /**
   *
   */
  public final File inputDB;
  /**
   * {@link File} object specifying the output directory.
   */
  public final File outputDirectory;
  /**
   *
   */
  public final ArrayList<String> features = new ArrayList<>();
  /**
   * Boolean value to indicate if debugging mode is enabled.
   */
  public final boolean debug;

  /**
   * Constructor of the {@link CLIParametersGenerateMSA} class.
   * <p>
   * Used to parse and validate the user command line input. Parsed arguments are stored - and accessible by other
   * components - via class properties.
   *
   * @param args {@link String} {@link Array} containing the command line arguments.
   * @throws ParseException     If any parsing error occurs.
   */
  public CLIParametersGenerateMSA(String[] args)
      throws ParseException, MusialCLException {
    // Store original command line arguments.
    arguments = args;
    // Add option to print help message.
    Options helpOptions = new Options();
    helpOptions.addOption("h", "help", false, "Display help information.");
    // Add non-help options.
    Options options = new Options();
    options.addOption("h", "help", false, "Display help information.");
    options.addOption(Option.builder("i")
        .longOpt("inputDB")
        .desc("TODO.")
        .hasArg()
        .required()
        .build());
    options.addOption(Option.builder("f")
        .longOpt("features")
        .desc("TODO.")
        .hasArg()
        .required()
        .build());
    options.addOption(Option.builder("o")
        .longOpt("output")
        .desc("Path to the directory at which results files shall be generated.")
        .required()
        .hasArg()
        .build());
    // Add options used to specify input data.
    options.addOption("v", "verbose", false, "Verbose mode on, i.e. print stacktrace on errors. [false]");
    // Instantiate a formatter for the help message and a default command line parser.
    HelpFormatter helpformatter = new HelpFormatter();
    CommandLineParser parser = new DefaultParser();
    CommandLine cmd;
    // First it is checked if the user specified only the -h option.
    try {
      cmd = parser.parse(helpOptions, args);
      if (cmd.hasOption('h')) {
        printHelp(options, helpformatter);
        System.exit(0);
      }
    } catch (ParseException ignored) {
    }
    // If the user has specified the -h option together with other options the same behaviour is applied.
    try {
      cmd = parser.parse(options, args);
      if (cmd.hasOption('h')) {
        printHelp(options, helpformatter);
        System.exit(0);
      }
    } catch (ParseException ignored) {
    }
    // Else it is tried to parse and validate command line options.
    cmd = parser.parse(options, args);
    outputDirectory = new File(cmd.getOptionValue('o'));
    boolean outputDirectoryExists = outputDirectory.exists();
    if (!outputDirectoryExists) {
      outputDirectoryExists = outputDirectory.mkdirs();
    }
    if (!outputDirectory.isDirectory()) {
      throw new MusialCLException("`-o` The specified output directory is no directory.");
    }
    if (!outputDirectoryExists) {
      throw new MusialCLException("`-o` The specified output directory does not exist and could not be created.");
    }
    // In the following each possible command line argument is validated, if any faulty input is detected the
    // application will exit with an exception.
    debug = cmd.hasOption("v");
    // Validation of data input parameters:
    inputDB = new File(cmd.getOptionValue("i"));
    // Validate optional behaviour options:
    features.addAll(Arrays.asList(cmd.getOptionValue("f").split(",")));
  }
}