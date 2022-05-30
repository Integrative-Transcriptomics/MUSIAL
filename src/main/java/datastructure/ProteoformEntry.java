package datastructure;

import java.util.HashMap;
import java.util.TreeSet;

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
