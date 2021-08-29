package exceptions;

/**
 * Command line argument exception of MUSIAL.
 * <p>
 * This exception is thrown if any error occurs during parsing and validating the user specified command line arguments.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class MusialCLAException extends Exception {

  /**
   * Constructor that accepts a String message.
   *
   * @param message The message comprised by the exception.
   */
  public MusialCLAException(String message) {
    super(message);
  }
}
