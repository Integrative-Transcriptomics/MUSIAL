package exceptions;

/**
 * Exception thrown if any data manipulation procedure is not logical in the context of the software, for example the
 * assignment of a protein resource to a whole genome reference feature.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public final class MusialIntegrityException extends Exception {

    /**
     * Constructor that accepts a String message.
     *
     * @param message The message comprised by the exception.
     */
    public MusialIntegrityException(String message) {
        super(message);
    }
}
