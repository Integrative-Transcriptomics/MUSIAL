package datastructure;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.TreeSet;
import utility.Bio;

/**
 *
 */
public class ProteoformEntry {

  /**
   *
   */
  public final String name;
  /**
   *
   */
  public final HashMap<String, String> annotations = new HashMap<>();
  /**
   *
   */
  public final TreeSet<String> samples = new TreeSet<>();

  /**
   *
   */
  public ProteoformEntry(String name, String sId, String proteoformVSwab) {
    this.name = name;
    this.annotations.put("vSwab", proteoformVSwab);
    this.samples.add(sId);
  }

}
