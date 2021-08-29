package exceptions;

/**
 * Faulty data exception of MUSIAL.
 * <p>
 * This exception is thrown if any data generation or storage is incorrect in the context of the software, i.e. gene
 * features with negative locations.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class MusialFaultyDataException extends Exception {

  /**
   * Constructor that accepts a String message.
   *
   * @param message The message comprised by the exception.
   */
  public MusialFaultyDataException(String message) {
    super(message);
  }
}
