package datastructure;

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
 * @TODO: 29.08.2021: This class is highly redundant compared to {@link GeneFeature} class and may extend it.
 * @since 2.0
 */
public final class ReferenceAnalysisEntry {

  /**
   * The internal name of the entry.
   */
  public final String analysisIdentifier;
  /**
   * The reference DNA sequence of the entry.
   */
  public final String analysisSequence;
  /**
   * The location of the entry on the reference, i.e. for genes the contig or chromosome the feature is found on.
   */
  public final String analysisSequenceLocation;
  /**
   * The 1-based indexed starting position of the feature.
   */
  public final int analysisSequenceStart;
  /**
   * The 1-based indexed end position of the feature.
   */
  public final int analysisSequenceEnd;
  /**
   * Indicates if the feature is located on the sense strand.
   */
  public final boolean isSense;

  /**
   * Constructor of {@link ReferenceAnalysisEntry}.
   *
   * @param entryIdentifier {@link String} representing the internal name of the reference feature to analyze.
   * @param entrySequence   {@link String} the full length sequence of the reference feature.
   * @param entryLocation   {@link String} the name of the reference location (contig, chromosome, plasmid) the feature
   *                        is located on.
   * @param entryStart      {@link Integer} The 1-based indexed starting position of the feature on the reference.
   * @param entryEnd        {@link Integer} The 1-based indexed end position of the feature on the reference.
   * @param isSense         {@link Boolean} indicating the strandedness of the feature.
   */
  public ReferenceAnalysisEntry(String entryIdentifier, String entrySequence, String entryLocation, int entryStart,
                                int entryEnd, boolean isSense) {
    this.analysisIdentifier = entryIdentifier;
    this.analysisSequence = entrySequence;
    this.analysisSequenceLocation = entryLocation;
    this.analysisSequenceStart = entryStart;
    this.analysisSequenceEnd = entryEnd;
    this.isSense = isSense;
  }

}