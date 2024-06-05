package exceptions;

/**
 * Exception thrown if any data generation or storage process is incorrect in the context of the software, e.g. gene
 * features with negative lengths.
 *
 * @author Simon Hackl
 */
public final class MusialException extends Exception {

    /**
     * Constructor that accepts a string message.
     *
     * @param message The message comprised by the exception.
     */
    public MusialException(String message) {
        super(message);
    }
}
