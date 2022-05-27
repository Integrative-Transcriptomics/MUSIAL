package datastructure;

import java.util.HashMap;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;

public class AminoacidVariantAnnotationEntry {

  public final HashMap<String, String> annotations = new HashMap<>();

  public final ConcurrentSkipListSet<String> occurrence = new ConcurrentSkipListSet<>();

  public AminoacidVariantAnnotationEntry() {

  }

}
