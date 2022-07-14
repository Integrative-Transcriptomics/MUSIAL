package exceptions;

/**
 * Exception thrown to indicate an error during I/O operations of MUSIAL internal methods.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class MusialIOException extends Exception {

    /**
     * Constructor that accepts a String message.
     *
     * @param message The message comprised by the exception.
     */
    public MusialIOException(String message) {
        super(message);
    }
}
