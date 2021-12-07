package datastructure;

import utility.IO;

/**
 * Container class to store counts regarding different call types of one analyzed position.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public class PositionStatistics {

  /**
   * Count of alternative A calls at the position.
   */
  public int ACalls = 0;
  /**
   * Count of alternative C calls at the position.
   */
  public int CCalls = 0;
  /**
   * Count of alternative G calls at the position.
   */
  public int GCalls = 0;
  /**
   * Count of alternative T calls at the position.
   */
  public int TCalls = 0;
  /**
   * Count of deletion calls at the position.
   */
  public int deletions = 0;
  /**
   * Count of no-calls at the position.
   */
  public int noCalls = 0;
  /**
   * Count of rejected calls at the position.
   */
  public int rejectedCalls = 0;

  /**
   * Constructor of {@link PositionStatistics}.
   */
  public PositionStatistics() {
  }

  /**
   * Returns a {@link String} yielding the tab-separated names of the single counts.
   * <p>
   * Used for the generation of output files.
   *
   * @return {@link String} yielding tab-separated count names.
   */
  public static String getHeaderString() {
    return "Position\tA_Calls\tT_Calls\tC_Calls\tG_Calls\tDeletions\tNo_Calls\tRejected_Calls" + IO.LINE_SEPARATOR;
  }

  /**
   * Returns a {@link String} yielding the tab-separated values of the single counts.
   * <p>
   * Used for the generation of output files.
   *
   * @return {@link String} yielding tab-separated count values.
   */
  public String getContentString() {
    return ACalls + "\t" + TCalls + "\t" + CCalls + "\t" + GCalls + "\t" +
        deletions + "\t" + noCalls + "\t" + rejectedCalls + IO.LINE_SEPARATOR;
  }

}
