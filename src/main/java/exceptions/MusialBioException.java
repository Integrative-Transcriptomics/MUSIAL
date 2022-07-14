package exceptions;

/**
 * Exception thrown if any data generation or storage process is incorrect in the context of the software, e.g. gene
 * features with negative lengths.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class MusialBioException extends Exception {

    /**
     * Constructor that accepts a String message.
     *
     * @param message The message comprised by the exception.
     */
    public MusialBioException(String message) {
        super(message);
    }
}
