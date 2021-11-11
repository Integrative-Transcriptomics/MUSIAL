package datastructure;

import exceptions.MusialFaultyDataException;

/**
 * Internal representation of a reference sequence location that is subject to analysis. May represent the full
 * genome, a single gene, contigs or plasmids and chromosomes.
 * <p>
 * Pairs of instances of this class and {@link SampleAnalysisEntry} instances are used internally as so called 'Run
 * entries' to specify a set of analysis tasks, i.e. which sample is analyzed with respect to which specified
 * reference feature.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class FeatureAnalysisEntry {

  /**
   * The internal name of the entry.
   */
  public final String identifier;
  /**
   * The reference DNA sequence of the entry.
   */
  private String referenceSequence;
  /**
   *
   */
  private int referenceSequenceLength;
  /**
   * The location of the entry on the reference, i.e. for genes the contig or chromosome the feature is found on.
   */
  public final String referenceSequenceLocation;
  /**
   * The 1-based indexed starting position of the feature.
   */
  public final int locationStart;
  /**
   * The 1-based indexed end position of the feature.
   */
  public final int locationEnd;
  /**
   * Indicates if the feature is located on the sense strand.
   */
  public final boolean isSense;

  /**
   * Constructor of {@link FeatureAnalysisEntry}.
   *
   * @param entryIdentifier {@link String} representing the internal name of the reference feature to analyze.
   * @param entryLocation   {@link String} the name of the reference location (contig, chromosome, plasmid) the feature
   *                        is located on.
   * @param entryStart      {@link Integer} The 1-based indexed starting position of the feature on the reference.
   * @param entryEnd        {@link Integer} The 1-based indexed end position of the feature on the reference.
   * @throws MusialFaultyDataException If the specified locus is ambiguous.
   */
  public FeatureAnalysisEntry(String entryIdentifier, String entryLocation, int entryStart,
                              int entryEnd) throws MusialFaultyDataException {
    this.identifier = entryIdentifier;
    this.referenceSequenceLocation = entryLocation;
    if (entryStart > 0 && entryEnd > 0 && (entryEnd - entryStart) != 0) {
      // CASE: Feature is on sense strand.
      this.isSense = true;
      this.locationStart = entryStart + 1;
      this.locationEnd = entryEnd;
    } else if (entryStart < 0 && entryEnd < 0 && (entryEnd - entryStart) != 0) {
      // CASE: Feature is on anti-sense strand.
      this.isSense = false;
      this.locationStart = -entryStart + 1;
      this.locationEnd = -entryEnd;
    } else {
      throw new MusialFaultyDataException("It was tried to generate a gene feature with faulty position " +
          "data:\t" + entryStart + ", " + entryEnd);
    }
  }

  /**
   * @return The stored reference sequence.
   */
  public String getReferenceSequence() {
    return referenceSequence;
  }

  /**
   * Sets the reference sequence.
   *
   * @param referenceSequence A {@link String} representing the reference sequence.
   */
  public void setReferenceSequence(String referenceSequence) {
    this.referenceSequence = referenceSequence;
  }

  /**
   * @return The total length of the parent sequence (if a gene feature is stored) or stored sequence. Used to
   * compute positions on the reverse complement.
   */
  public int getReferenceSequenceLength() {
    return referenceSequenceLength;
  }

  /**
   * Sets the total length of the parent sequence (if a gene feature is stored) or stored sequence. Used to
   * compute positions on the reverse complement.
   *
   * @param referenceSequenceLength A {@link Integer} representing the sequence length to store.
   */
  public void setReferenceSequenceLength(int referenceSequenceLength) {
    this.referenceSequenceLength = referenceSequenceLength;
  }
}