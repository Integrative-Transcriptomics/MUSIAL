package datastructure;

import exceptions.MusialFaultyDataException;

/**
 * Represents and stores information about a gene feature parsed from a `.gff` file.
 * <p>
 * Internal use to store information about one gene of a reference sequence.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 * @TODO: 29.08.2021: Currently only the representation of features on the sense strand is implemented.
 */
public final class GeneFeature {

  /**
   * The name of the reference sequence the feature is located at.
   */
  public final String seqName;
  /**
   * The name of the feature.
   */
  public final String featureName;
  /**
   * The start position of the feature.
   */
  public final int startPosition;
  /**
   * The end position of the feature.
   */
  public final int endPosition;

  /**
   * Constructor of {@link GeneFeature}.
   *
   * @param seqName       {@link String} represents the reference sequence name the feature is located on.
   * @param featureName   {@link String} represents
   * @param startPosition {@link Integer} 1-based indexed starting position of the feature on the reference sequence.
   * @param endPosition   {@link Integer} 1-based indexed end position of the feature on the reference sequence.
   * @throws MusialFaultyDataException If a feature with faulty start and end position is tried to generate, i.e. if
   *                                   the start position is negative, both positions are equal or the start position is lower than the end position.
   */
  public GeneFeature(String seqName, String featureName, int startPosition, int endPosition) throws
      MusialFaultyDataException {
    if (startPosition < 1 || (startPosition > endPosition) || (endPosition - startPosition) == 0) {
      throw new MusialFaultyDataException(
          "It was tried to generate a `GeneFeature` with faulty position data:\t" + startPosition + ", " + endPosition);
    } else {
      this.seqName = seqName;
      this.featureName = featureName;
      this.startPosition = startPosition;
      this.endPosition = endPosition;
    }
  }

}
