package exceptions;

/**
 * Duplication exception of MUSIAL.
 * <p>
 * This exception is thrown if any duplication in key assignments, i.e. multiple samples with the same name or
 * multiple entries within one sample for the same genetic location, occur.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class MusialDuplicationException extends Exception {

  /**
   * Constructor that accepts a String message.
   *
   * @param message The message comprised by the exception.
   */
  public MusialDuplicationException(String message) {
    super(message);
  }
}
