package datastructure;

import exceptions.MusialBioException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Represents and stores information about a single genomic position of one sample relative to a reference (genetic
 * sequence).
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class VariantPosition {

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
   * Annotation tag to indicate low coverage.
   */
  public static final char ANN_SEP = '@';
  /**
   * The alternate content (base symbol) of the variant position.
   */
  private final char content;
  /**
   * The Phred-scaled quality score of the position.
   */
  private final double positionQuality;
  /**
   * The depth of coverage of the position
   */
  private final double positionCoverage;
  /**
   * The frequency of the alternate content of the variant position.
   */
  private final double frequency;
  /**
   * Additional annotation tags associated with the position.
   */
  private final ArrayList<String> annotationTags = new ArrayList<>();
  /**
   * Additional annotation key, value pairs associated with the position.
   */
  private final HashMap<String, String> annotationMap = new HashMap<>();
  /**
   * Whether the alleles reflect an insertion
   */
  private final ArrayList<Boolean> isInsertion = new ArrayList<>();
  /**
   * Whether the alleles reflect a deletion.
   */
  private final ArrayList<Boolean> isDeletion = new ArrayList<>();
  /**
   * Whether the instance reflects not a possible variant of a sample, but the content of the reference.
   */
  private final boolean isReference;

  /**
   * Constructor of {@link VariantPosition}.
   */
  public VariantPosition(ArrayList<String> contents, double quality, double coverage, ArrayList<Double> frequencies,
                         boolean isReference) throws MusialBioException {
    if ( contents.size() != frequencies.size() ) {
      throw new MusialBioException( "Construction of variant position failed: The number of contents and frequencies did " +
          "not match" );
    }
    if ( contents.size() == 0 ) {
      throw new MusialBioException( "Construction of variant position failed: An empty content list was passed." );
    }
    this.isReference = isReference;
    this.positionQuality = quality;
    this.positionCoverage = coverage;
    for (int i = 0; i < contents.size(); i++) {
      String content = contents.get( i );
      double frequency = frequencies.get( i );
      this.contentList.add( content );
      this.frequencyList.add( frequency );
    }

  }
}