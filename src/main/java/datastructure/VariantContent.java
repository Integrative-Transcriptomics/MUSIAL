package datastructure;

import java.util.HashMap;

/**
 * Stores information about the content of a single genomic position of one sample relative to a reference (genetic
 * sequence).
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class VariantContent {

  /**
   * Character to indicate a reference call, i.e. no variant at the position.
   */
  public static final char REFERENCE = '.';
  /**
   * Character to indicate a deletion relative to the reference.
   */
  public static final char DELETION = '-';
  /**
   * Character to indicate an alternative A nucleotide.
   */
  public static final char ALT_A = 'A';
  /**
   * Character to indicate an alternative C nucleotide.
   */
  public static final char ALT_C = 'C';
  /**
   * Character to indicate an alternative G nucleotide.
   */
  public static final char ALT_G = 'G';
  /**
   * Character to indicate an alternative T nucleotide.
   */
  public static final char ALT_T = 'T';
  /**
   * Character to indicate a no call position.
   */
  public static final char NO_CALL = 'N';
  /**
   * ...
   */
  public static final char INSERTION_DUMMY = '~';
  /**
   * Separator character used for generating annotation strings.
   */
  public static final char SEPARATOR = '@';
  /**
   * The alternate content (base symbol) of the variant position.
   */
  public final char content;
  /**
   * The Phred-scaled quality score of the position.
   */
  public final double quality;
  /**
   * The depth of coverage of the position
   */
  public final double coverage;
  /**
   * The frequency of the alternate content of the variant position.
   */
  public final double frequency;
  /**
   * Additional annotation key, value pairs associated with the position.
   */
  public final HashMap<String, String> annotations = new HashMap<>();
  /**
   * Indicates if the content is of the most frequent allele.
   */
  public final boolean MFA;

  /**
   * Constructor of {@link VariantContent}.
   */
  public VariantContent(char content, double quality, double coverage,double frequency, boolean isMFA) {
    this.content = content;
    this.quality = quality;
    this.coverage = coverage;
    this.frequency = frequency;
    this.MFA = isMFA;
  }

  /**
   * Adds a new annotation value accessible via the specified key to the annotations map of this object.
   *
   * @param key {@link String} representing the key of the annotation.
   * @param value {@link String} representing the value of the annotation.
   */
  public void addAnnotation( String key, String value ) {
    this.annotations.put( key, value );
  }
}