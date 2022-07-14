package exceptions;

/**
 * Exception thrown if any error occurs during parsing and validating the user specified command line arguments.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class MusialCLException extends Exception {

    /**
     * Constructor that accepts a String message.
     *
     * @param message The message comprised by the exception.
     */
    public MusialCLException(String message) {
        super(message);
    }
}
