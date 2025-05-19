package exceptions;

/**
 * A custom exception class for handling application-specific errors.
 * <p>
 * This final class extends {@link Exception} and is used to represent
 * errors specific to the Musial application. It provides a constructor
 * for initializing the exception with a custom error message.
 * </p>
 */
public final class MusialException extends Exception {

    /**
     * Constructs a new {@link MusialException} with the specified detail message.
     * <p>
     * This constructor initializes the exception with a custom error message
     * that can be retrieved later using the {@link Throwable#getMessage()} method.
     * </p>
     *
     * @param message The detail message for the exception.
     */
    public MusialException(String message) {
        super(message);
    }

}
