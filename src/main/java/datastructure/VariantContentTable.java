package datastructure;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;
import utility.VariantPositionComparator;

/**
 * Major data structure to store the results of Musial.
 * <p>
 * Implements a variety of thread safe data structures to store the information parsed from `.vcf` files. In addition
 * accessor methods to retrieve the stored data either as statistics or in a {@link String} format for the generation
 * of human readable output are implemented.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class VariantContentTable {

  /**
   * Thread safe table structure that stores all information parsed from `.vcf` files in the form of
   * {@link VariantContent} instances.
   * <p>
   * The first key represents the analyzed reference feature, e.g. in the case of a full genome analysis the internal
   * name of the genome or in the case of single genes the gene names.
   * <p>
   * The second key represents the name of the analyzed sample.
   * <p>
   * The third key represents the position on the reference as a {@link String}. In the most cases this will simply
   * be the {@link String} representation of the {@link Integer} position, however, to represent an insertion in any
   * sample `+i` is appended: I.e. if an insertion at position 10 of length 2 was analyzed, the keys 10+1
   * and 10+2 will be used. The third key set is sorted by the {@link VariantPositionComparator}. This system
   * allows for an easy insertion of `-` gap symbols in those samples in which no insertion occurred.
   * <p>
   * The final values stored are {@link ConcurrentLinkedQueue} instances that contain {@link VariantContent} objects,
   * i.e. multiple {@link VariantContent} objects are stored if a sample had a het. call at any position.
   * The concurrent linked queue implementation allows the {@link VariantContent} objects to stay sorted according to
   * their frequency.
   */
  private final ConcurrentHashMap<String, ConcurrentHashMap<String, ConcurrentSkipListMap<String,
      ConcurrentLinkedQueue<VariantContent>>>>
      table;
  /**
   * Thread safe structure mapping reference feature names to a list of positions at which a variant was called.
   */
  private final ConcurrentHashMap<String, ConcurrentSkipListSet<String>> variantPositions;
  /**
   * Thread safe structure mapping reference feature names to a list of {@link PositionStatistics} instances.
   */
  private final ConcurrentHashMap<String, ConcurrentHashMap<String, PositionStatistics>> positionStatistics;
  /**
   * Thread safe structure mapping reference feature names to a list of {@link SampleStatistics} instances.
   */
  private final ConcurrentHashMap<String, ConcurrentHashMap<String, SampleStatistics>> sampleStatistics;
  /**
   * Thread safe structure storing per reference feature names per sample the starting positions of consecutive
   * deletions.
   */
  private final ConcurrentHashMap<String, ConcurrentHashMap<String, ConcurrentHashMap<String, ArrayList<String>>>>
      deletions;
  /**
   * Number of samples, used for counting.
   */
  private final int noSamples;

  /**
   * Constructor of {@link VariantContentTable}.
   *
   * @param numberOfSamples {@link Integer} indicating the number of samples, used for counting.
   */
  public VariantContentTable(int numberOfSamples) {
    this.table = new ConcurrentHashMap<>();
    this.variantPositions = new ConcurrentHashMap<>();
    this.positionStatistics = new ConcurrentHashMap<>();
    this.sampleStatistics = new ConcurrentHashMap<>();
    this.deletions = new ConcurrentHashMap<>();
    this.noSamples = numberOfSamples;
  }

  /**
   * Inserts a new {@link VariantContent} instance into the table structure.
   *
   * @param referenceFeatureName {@link String} the name of the reference feature to access.
   * @param sampleName           {@link String} representing the analyzed sample.
   * @param position             {@link String} representing the analyzed position.
   * @param content              {@link Character} representing the nucleotide content to insert.
   * @param quality              {@link Double} representing the quality of the analyzed position.
   * @param coverage             {@link Double} representing the coverage of the analyzed position.
   * @param frequency            {@link Double} representing the frequency of the analyzed content.
   * @param isMfa                {@link Boolean} indicating if the content is from the most frequently observed allele.
   * @param additionalAnnotation {@link HashMap} containing key, value pairs of {@link String} that reflect
   *                             annotations associated with the variant content.
   */
  public void putVariantPosition(String referenceFeatureName, String sampleName, String position,
                                 char content, double quality, double coverage, double frequency, boolean isMfa,
                                 HashMap<String, String> additionalAnnotation) {
    initializeInternalContentStructure(referenceFeatureName, sampleName);
    initializeInternalSampleStatisticsStructure(referenceFeatureName, sampleName);
    initializeInternalPositionStatisticsStructure(referenceFeatureName, position);
    VariantContent variantContent = new VariantContent(content, quality, coverage, frequency, isMfa);
    if (!this.table.get(referenceFeatureName).get(sampleName).containsKey(position)) {
      this.table.get(referenceFeatureName).get(sampleName).put(position, new ConcurrentLinkedQueue<>());
    }
    this.table.get(referenceFeatureName).get(sampleName).get(position).add(variantContent);
    if (additionalAnnotation.size() > 0) {
      for (Map.Entry<String, String> entry : additionalAnnotation.entrySet()) {
        variantContent.addAnnotation(entry.getKey(), entry.getValue());
      }
    }
    // If the content of the passed variable position is no reference symbol and it does not represent the reference
    // sequence directly, the position is added as variant position.
    if (content != VariantContent.REFERENCE && content != VariantContent.NO_CALL && !sampleName.equals(
        "Reference")) {
      this.variantPositions.get(referenceFeatureName).add(position);
    }
  }

  /**
   * Increments the rejected call counter by one for the specified reference feature (i.e. accessed via
   * `referenceFeatureName`), position and sample.
   *
   * @param referenceFeatureName {@link String} specifying the reference feature for which the counter shall be
   *                             incremented.
   * @param sampleName           {@link String} specifying the sample for which the count shall be increased.
   * @param position             {@link String } specifying the position for which the count shall be increased.
   */
  public void countRejectedCall(String referenceFeatureName, String sampleName, String position) {
    initializeInternalSampleStatisticsStructure(referenceFeatureName, sampleName);
    initializeInternalPositionStatisticsStructure(referenceFeatureName, position);
    // Count no call for position statistics.
    this.positionStatistics.get(referenceFeatureName).get(position).rejectedCalls += 1;
    this.positionStatistics.get(referenceFeatureName).get(position).noCalls -= 1;
    // Count no call for sample statistics.
    this.sampleStatistics.get(referenceFeatureName).get(sampleName).rejectedCalls += 1;
  }

  /**
   * Increments the respective call type counter by one for the specified reference feature (i.e. accessed via
   * `referenceFeatureName`) and sample, dependent on the value of `isHet`.
   *
   * @param referenceFeatureName {@link String} specifying the reference feature for which the counter shall be
   *                             incremented.
   * @param sampleName           {@link String} specifying the sample for which the count shall be increased.
   * @param isHet                {@link Boolean } specifying the call type, i.e. false for hom. and true for het. calls.
   */
  public void countCallType(String referenceFeatureName, String sampleName, boolean isHet) {
    initializeInternalSampleStatisticsStructure(referenceFeatureName, sampleName);
    if (isHet) {
      this.sampleStatistics.get(referenceFeatureName).get(sampleName).hetCalls += 1;
    } else {
      this.sampleStatistics.get(referenceFeatureName).get(sampleName).homCalls += 1;
    }
  }

  /**
   * Increments the variant call counter by one for the specified reference feature (i.e. accessed via
   * `referenceFeatureName`), position and sample. Which counter is increased depends on the specified `content` and
   * `position`.
   *
   * @param referenceFeatureName {@link String} specifying the reference feature for which the counter shall be
   *                             incremented.
   * @param sampleName           {@link String} specifying the sample for which the count shall be increased.
   * @param position             {@link String } specifying the position for which the count shall be increased.
   * @param content              {@link Character} specifying the variants nucleotide content.
   */
  public void countAlternateCall(String referenceFeatureName, String sampleName, String position, char content) {
    initializeInternalSampleStatisticsStructure(referenceFeatureName, sampleName);
    initializeInternalPositionStatisticsStructure(referenceFeatureName, position);
    // Count no call for position statistics.
    switch (content) {
      case VariantContent.ALT_A -> {
        this.positionStatistics.get(referenceFeatureName).get(position).ACalls += 1;
        this.positionStatistics.get(referenceFeatureName).get(position).noCalls -= 1;
        if (!position.contains("+")) {
          this.sampleStatistics.get(referenceFeatureName).get(sampleName).SNV += 1;
        }
      }
      case VariantContent.ALT_T -> {
        this.positionStatistics.get(referenceFeatureName).get(position).TCalls += 1;
        this.positionStatistics.get(referenceFeatureName).get(position).noCalls -= 1;
        if (!position.contains("+")) {
          this.sampleStatistics.get(referenceFeatureName).get(sampleName).SNV += 1;
        }
      }
      case VariantContent.ALT_C -> {
        this.positionStatistics.get(referenceFeatureName).get(position).CCalls += 1;
        this.positionStatistics.get(referenceFeatureName).get(position).noCalls -= 1;
        if (!position.contains("+")) {
          this.sampleStatistics.get(referenceFeatureName).get(sampleName).SNV += 1;
        }
      }
      case VariantContent.ALT_G -> {
        this.positionStatistics.get(referenceFeatureName).get(position).GCalls += 1;
        this.positionStatistics.get(referenceFeatureName).get(position).noCalls -= 1;
        if (!position.contains("+")) {
          this.sampleStatistics.get(referenceFeatureName).get(sampleName).SNV += 1;
        }
      }
      case VariantContent.DELETION -> {
        this.positionStatistics.get(referenceFeatureName).get(position).deletions += 1;
        this.positionStatistics.get(referenceFeatureName).get(position).noCalls -= 1;
        this.sampleStatistics.get(referenceFeatureName).get(sampleName).deletedPositions += 1;
      }
    }
    if (position.contains("+")) {
      this.sampleStatistics.get(referenceFeatureName).get(sampleName).insertedPositions += 1;
    }
  }

  /**
   * Inserts empty {@link ConcurrentHashMap} and {@link ConcurrentSkipListMap} structures into the table if none are
   * present for a specified reference analysis entry and sample name.
   * <p>
   * This method is intended to be used for samples with only confident reference calls as these will not be
   * considered for the output files otherwise.
   *
   * @param referenceFeatureName {@link String} the name of the reference feature to access.
   * @param sampleName           {@link String} representing the name of the sample.
   */
  public void addSampleToReference(String referenceFeatureName, String sampleName) {
    // Insert new ConcurrentHashMap for reference feature if none is present.
    if (!this.table.containsKey(referenceFeatureName)) {
      this.table.put(referenceFeatureName, new ConcurrentHashMap<>());
    }
    // Insert new ConcurrentSkipListMap for sample name if none is present.
    if (!this.table.get(referenceFeatureName).containsKey(sampleName)) {
      this.table.get(referenceFeatureName).put(sampleName,
          new ConcurrentSkipListMap<>(new VariantPositionComparator()));
    }
    // Initialize empty sample statistics container.
    initializeInternalSampleStatisticsStructure(referenceFeatureName, sampleName);
  }

  /**
   * Returns an {@link Iterator} of {@link VariantContent} objects, if stored in the table.
   * <p>
   * Returns null if no entry exists.
   *
   * @param referenceFeatureName {@link String} representing the reference analysis identifier / feature identifier to
   *                             access.
   * @param sampleName           {@link String} representing the sample identifier to access.
   * @param position             {@link String} representing the position to access.
   * @return {@link Iterator} of {@link VariantContent} objects, if stored in the table; else null.
   */
  public Iterator<VariantContent> getContent(String referenceFeatureName, String sampleName, String position) {
    if (this.table.get(referenceFeatureName).get(sampleName).containsKey(position)) {
      return this.table.get(referenceFeatureName).get(sampleName).get(position).iterator();
    } else {
      return null;
    }
  }

  /**
   * Returns all stored sample names for one reference location name.
   *
   * @param referenceFeatureName {@link String} the name of the reference location.
   * @return {@link HashSet<String>} containing all sample names of the specified reference location.
   */
  public HashSet<String> getSampleNames(String referenceFeatureName) {
    HashSet<String> sampleNames = new HashSet<>();
    Iterator<String> keyIterator = this.table.get(referenceFeatureName).keys().asIterator();
    while (keyIterator.hasNext()) {
      String sampleName = keyIterator.next();
      if (!sampleName.equals("Reference")) {
        sampleNames.add(sampleName);
      }
    }
    return sampleNames;
  }

  /**
   * Returns the {@link SampleStatistics} instance stored for the specified reference location and sample name.
   *
   * @param referenceFeatureName {@link String} the name of the reference location.
   * @param sampleName           {@link String} the name of the sample.
   * @return {@link SampleStatistics} instance.
   */
  public SampleStatistics getSampleStatistics(String referenceFeatureName, String sampleName) {
    return this.sampleStatistics.get(referenceFeatureName).get(sampleName);
  }

  /**
   * Returns the {@link PositionStatistics} instance stored for the specified reference location and position.
   *
   * @param referenceFeatureName {@link String} the name of the reference location.
   * @param position             {@link String} the identifier of the position.
   * @return {@link PositionStatistics} instance.
   */
  public PositionStatistics getPositionStatistics(String referenceFeatureName, String position) {
    return this.positionStatistics.get(referenceFeatureName).get(position);
  }

  /**
   * Returns a navigable set of all (including non-variant) positions of one reference feature.
   *
   * @param referenceFeatureName {@link String} the reference feature to access.
   * @return {@link TreeSet<String>} containing all positions of the specified reference feature.
   */
  public NavigableSet<String> getAllPositionsSet(String referenceFeatureName) {
    ConcurrentSkipListSet<String> positions = new ConcurrentSkipListSet<>();
    for (String position : this.table.get(referenceFeatureName).get("Reference").navigableKeySet()) {
      //noinspection UseBulkOperation
      positions.add(position);
    }
    ConcurrentSkipListSet<String> variantPositions = this.getVariantPositionsSet(referenceFeatureName);
    for (String p : variantPositions) {
      //noinspection UseBulkOperation
      positions.add(p);
    }
    return positions;
  }

  /**
   * Returns an iterator of all (including non-variant) positions of one reference feature.
   *
   * @param referenceFeatureName {@link String} the reference feature to access.
   * @return {@link Iterator<String>} over all positions of the specified reference feature.
   */
  @SuppressWarnings("unused")
  public Iterator<String> getAllPositionsIterator(String referenceFeatureName) {
    NavigableSet<String> allPositions = getAllPositionsSet(referenceFeatureName);
    return allPositions.iterator();
  }

  /**
   * Returns a navigable set of all variant positions of one reference feature.
   *
   * @param referenceFeatureName {@link String} the reference feature to access.
   * @return {@link ConcurrentSkipListSet<String>} containing all variant positions of the specified reference feature.
   */
  public ConcurrentSkipListSet<String> getVariantPositionsSet(String referenceFeatureName) {
    if (this.variantPositions.containsKey(referenceFeatureName)) {
      return this.variantPositions.get(referenceFeatureName);
    } else {
      return new ConcurrentSkipListSet<>();
    }
  }

  /**
   * Returns an iterator over all variant positions of a reference feature.
   *
   * @param referenceFeatureName {@link String} the reference feature to access.
   * @return {@link Iterator<String>} over all variant positions of the specified reference feature.
   */
  public Iterator<String> getVariantPositions(String referenceFeatureName) {
    if (this.variantPositions.containsKey(referenceFeatureName)) {
      return this.variantPositions.get(referenceFeatureName).iterator();
    } else {
      return Collections.emptyIterator();
    }
  }

  /**
   * Inserts a new block of consecutively deleted positions with respect to the specified sample and reference
   * location starting exclusively at the specified position.
   *
   * @param referenceFeatureName          {@link String} the name of the reference location.
   * @param sampleName                    {@link String} the name of the sample.
   * @param position                      {@link String} the identifier of the position.
   * @param consecutivelyDeletedPositions {@link ArrayList} of {@link String}s representing consecutively deleted
   *                                      positions in the specified sample and reference location.
   */
  public void putDeletion(String referenceFeatureName, String sampleName, String position,
                          ArrayList<String> consecutivelyDeletedPositions) {
    if (!this.deletions.containsKey(referenceFeatureName)) {
      this.deletions.put(referenceFeatureName, new ConcurrentHashMap<>());
    }
    if (!this.deletions.get(referenceFeatureName).containsKey(sampleName)) {
      this.deletions.get(referenceFeatureName).put(sampleName, new ConcurrentHashMap<>());
    }
    this.deletions.get(referenceFeatureName).get(sampleName).put(position, consecutivelyDeletedPositions);
  }

  /**
   * Returns a block of consecutively deleted positions with respect to the specified sample and reference
   * location starting exclusively at the specified position or null if no such entry is present.
   *
   * @param referenceFeatureName {@link String} the name of the reference location.
   * @param sampleName           {@link String} the name of the sample.
   * @param position             {@link String} the identifier of the position.
   * @return {@link ArrayList} if the position is a starting point of a block of consecutively deleted positions or
   * null.
   */
  public ArrayList<String> getDeletion(String referenceFeatureName, String sampleName, String position) {
    if (!this.deletions.containsKey(referenceFeatureName) ||
        !this.deletions.get(referenceFeatureName).containsKey(sampleName) ||
        !this.deletions.get(referenceFeatureName).get(sampleName).containsKey(position)) {
      return null;
    } else {
      return this.deletions.get(referenceFeatureName).get(sampleName).get(position);
    }
  }

  /**
   * Internal method to initialize objects used for sample statistic counting.
   *
   * @param referenceFeatureName {@link String} specifying the reference feature for which internal objects are
   *                             generated.
   * @param sampleName           {@link String} specifying the sample for which internal objects are generated.
   */
  private void initializeInternalSampleStatisticsStructure(String referenceFeatureName, String sampleName) {
    // Insert new ConcurrentHashMap for sample statistics of reference, if none is present.
    if (!this.sampleStatistics.containsKey(referenceFeatureName)) {
      this.sampleStatistics.put(referenceFeatureName, new ConcurrentHashMap<>());
    }
    // Insert new SampleStatistics instance for sample of reference, if none is present.
    if (!this.sampleStatistics.get(referenceFeatureName).containsKey(sampleName)) {
      SampleStatistics sampleStatistics = new SampleStatistics();
      this.sampleStatistics.get(referenceFeatureName).put(sampleName, sampleStatistics);
    }
  }

  /**
   * Internal method to initialize objects used for position statistic counting.
   *
   * @param referenceFeatureName {@link String} specifying the reference feature for which internal objects are
   *                             generated.
   * @param position             {@link String} specifying the position for which internal objects are generated.
   */
  private void initializeInternalPositionStatisticsStructure(String referenceFeatureName, String position) {
    // Insert new ConcurrentHashMap for position statistics of reference, if none is present.
    if (!this.positionStatistics.containsKey(referenceFeatureName)) {
      this.positionStatistics.put(referenceFeatureName, new ConcurrentHashMap<>());
    }
    // Insert new PositionStatistics instance for position of reference, if none is present.
    if (!this.positionStatistics.get(referenceFeatureName).containsKey(position)) {
      PositionStatistics positionStatistics = new PositionStatistics();
      positionStatistics.noCalls = this.noSamples;
      this.positionStatistics.get(referenceFeatureName).put(position, positionStatistics);
    }
  }

  /**
   * Internal method to initialize objects used for variant content storage.
   *
   * @param referenceFeatureName {@link String} specifying the reference feature for which internal objects are
   *                             generated.
   * @param sampleName           {@link String} specifying the sample for which internal objects are generated.
   */
  private void initializeInternalContentStructure(String referenceFeatureName, String sampleName) {
    // Insert new ConcurrentHashMaps for reference feature if none is present.
    if (!this.table.containsKey(referenceFeatureName)) {
      this.table.put(referenceFeatureName, new ConcurrentHashMap<>());
    }
    // Insert new ConcurrentSkipListMap for sample name if none is present.
    if (!this.table.get(referenceFeatureName).containsKey(sampleName)) {
      this.table.get(referenceFeatureName).put(sampleName,
          new ConcurrentSkipListMap<>(new VariantPositionComparator()));
    }
    // Insert new ConcurrentSkipListSet to store variant positions if none is present.
    if (!this.variantPositions.containsKey(referenceFeatureName)) {
      this.variantPositions
          .put(referenceFeatureName, new ConcurrentSkipListSet<>(new VariantPositionComparator()));
    }
  }

}