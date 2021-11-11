package datastructure;

/**
 * Represents and stores information about a single genomic position of one sample relative to a reference (genetic
 * sequence).
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class VariablePosition {

  /**
   * Indicates reference base call, no variation in this position.
   */
  public static final char REFERENCE = '.';
  /**
   * Indicates reference base call, but at least one filter criteria was not met and the call was discarded.
   */
  public static final char REFERENCE_DISCARDED = ':';
  /**
   * Represents a deletion (relative to the reference).
   */
  public static final char DELETION = '-';
  /**
   * Indicates alternative base call as A.
   */
  public static final char ALT_A = 'A';
  /**
   * Indicates alternative base call as A, but at least one filter criteria was not met and the call was discarded.
   */
  public static final char ALT_A_DISCARDED = 'a';
  /**
   * Indicates alternative base call as C.
   */
  public static final char ALT_C = 'C';
  /**
   * Indicates alternative base call as C, but at least one filter criteria was not met and the call was discarded.
   */
  public static final char ALT_C_DISCARDED = 'c';
  /**
   * Indicates alternative base call as G.
   */
  public static final char ALT_G = 'G';
  /**
   * Indicates alternative base call as G, but at least one filter criteria was not met and the call was discarded.
   */
  public static final char ALT_G_DISCARDED = 'g';
  /**
   * Indicates alternative base call as T.
   */
  public static final char ALT_T = 'T';
  /**
   * Indicates alternative base call as T, but at least one filter criteria was not met and the call was discarded.
   */
  public static final char ALT_T_DISCARDED = 't';
  /**
   * Indicates that no base call was possible.
   */
  public static final char NO_CALL = 'N';

  /**
   * The content (base symbol) of the variable position.
   */
  public final Character content;
  /**
   * The annotation string of the variable position.
   */
  public final String annotation;

  /**
   * Constructs a {@link VariablePosition} with meta-information.
   *
   * @param content    {@link Character} representing the variants base change.
   * @param annotation Any annotations for the alternative.
   */
  public VariablePosition(char content, String annotation) {
    this.content = content;
    this.annotation = annotation;
  }

}