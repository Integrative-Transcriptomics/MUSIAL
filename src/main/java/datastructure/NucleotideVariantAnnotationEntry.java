package datastructure;

import java.util.HashMap;
import java.util.concurrent.ConcurrentSkipListMap;

public class NucleotideVariantAnnotationEntry {

  public final HashMap<String, String> annotations = new HashMap<>();

  public final ConcurrentSkipListMap<String, String> occurrence = new ConcurrentSkipListMap<>();

  public NucleotideVariantAnnotationEntry() {

  }

  public static String constructSampleSpecificAnnotation(boolean rejected, boolean primary, double quality,
                                                         double frequency, double coverage) {
    return rejected + "|" + primary + "|" + quality + "|" + frequency + "|" + coverage;
  }

}
