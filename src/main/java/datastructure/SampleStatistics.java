package datastructure;

import utility.IO;

/**
 * Container class to store counts regarding different call types of one analyzed sample.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public class SampleStatistics {

  /**
   * Count of hom. calls with respect to the sample.
   */
  public int homCalls = 0;
  /**
   * Count of het. calls with respect to the sample.
   */
  public int hetCalls = 0;
  /**
   * Count of rejected calls with respect to the sample.
   */
  public int rejectedCalls = 0;
  /**
   * Count of deleted positions with respect to the sample.
   */
  public int deletedPositions = 0;
  /**
   * Count of inserted positions with respect to the sample.
   */
  public int insertedPositions = 0;
  /**
   * Count of called single nucleotide variants with respect to the sample.
   */
  public int SNV = 0;

  /**
   * Constructor of {@link SampleStatistics}.
   */
  public SampleStatistics() { }

  /**
   * Returns a {@link String} yielding the tab-separated names of the single counts.
   * <p>
   * Used for the generation of output files.
   *
   * @return {@link String} yielding tab-separated count names.
   */
  public static String getHeaderString() {
    return "Sample\tHomozygous_Calls\tHeterozygous_Calls\tRejected_Calls\tDeleted_Positions\tInserted_Positions" +
        "\tSNVs" + IO.LINE_SEPARATOR;
  }

  /**
   * Returns a {@link String} yielding the tab-separated values of the single counts.
   * <p>
   * Used for the generation of output files.
   *
   * @return {@link String} yielding tab-separated count values.
   */
  public String getContentString() {
    return homCalls + "\t" + hetCalls + "\t" + rejectedCalls + "\t" +
        deletedPositions + "\t" + insertedPositions +
        "\t" + SNV + IO.LINE_SEPARATOR;
  }

}
