package cli;

import java.lang.reflect.Array;
import main.Musial;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

/**
 * Interface for command line interface argument parsers.
 * <p>
 * Each MUSIAL module uses a distinct command line interface parser extending this interface.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public interface CLIParameters {

  /**
   * Print help information for command line interface arguments.
   *
   * @param options       The specified command line interface {@link Options}.
   * @param helpFormatter A {@link HelpFormatter} for displaying help information.
   */
  default void printHelp(Options options, HelpFormatter helpFormatter) {
    helpFormatter.printHelp("java -jar MUSIAL-" + Musial.VERSION + ".jar generateVDict", options);
  }

}
