package exceptions;

/**
 * IO exception of MUSIAL.
 * <p>
 * This exception is thrown if any generation of
 *
 * @author Simon Hackl
 * @version 2.0
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
