package utility;

import java.util.Comparator;

/**
 * Implementation of the {@link Comparator<String>} interface to sort {@link String} representations of genome
 * positions that may contain "I" characters to represent insertions according to their integer value and number of
 * contained "I"s.
 * <p>
 * For example the strings "3","5","2","4","4I","4II","4III" are sorted as "2","3","4","4I","4II","4III","5" and
 * would represent an insertion of length three at position 4.
 *
 * @author Alexander Seitz
 * @author Simon Hackl
 * @version 2.0
 */
public final class VariablePositionTableComparator implements Comparator<String> {

  /**
   * Implementation of the interfaces compare method.
   *
   * @param vp1 The first {@link String} to compare.
   * @param vp2 The second {@link String} to compare.
   * @return {@link Integer}: 0 if vp1 and vp2 are equal, 1 if vp1 is greater than vp2 and -1 if vp1 is lower than vp2.
   */
  @Override
  public int compare(String vp1, String vp2) {
    if (vp1.equals(vp2)) {
      return 0;
    } else {
      if (vp1.contains("I") && vp2.contains("I")) {
        int vp1Pos = getIntegerValue(vp1);
        int vp2Pos = getIntegerValue(vp2);
        if (vp1Pos == vp2Pos) {
          int vp1Is = vp1.substring(vp1.indexOf("I")).length();
          int vp2Is = vp2.substring(vp2.indexOf("I")).length();
          return Integer.compare(vp1Is, vp2Is);
        } else {
          return Integer.compare(vp1Pos, vp2Pos);
        }
      } else if (vp1.contains("I") && !vp2.contains("I")) {
        int vp1Pos = getIntegerValue(vp1);
        int vp2Pos = Integer.parseInt(vp2);
        return Integer.compare(vp1Pos, vp2Pos);
      } else if (!vp1.contains("I") && vp2.contains("I")) {
        int vp1Pos = Integer.parseInt(vp1);
        int vp2Pos = getIntegerValue(vp2);
        return Integer.compare(vp1Pos, vp2Pos);
      } else {
        int vp1Pos = Integer.parseInt(vp1);
        int vp2Pos = Integer.parseInt(vp2);
        return Integer.compare(vp1Pos, vp2Pos);
      }
    }
  }

  /**
   * Parses the integer value of a genome position with a "I"s suffix.
   *
   * @param vp {@link String} of the form {0,1,2,3,4,5,6,7,8,9}*{I}*.
   * @return The {@link Integer} value of the passed string.
   */
  private int getIntegerValue(String vp) {
    return Integer.parseInt(vp.substring(0, vp.indexOf("I")));
  }
}