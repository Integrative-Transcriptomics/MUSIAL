package cli;

import org.apache.commons.cli.Options;

/**
 * Interface representing a Command Line Interface (CLI) component.
 */
public interface CLI {

    /**
     * Returns the command line options for this CLI component.
     *
     * @return An {@link Options} object containing the command line options.
     */
    static Options options() {
        return new Options();
    }

}
