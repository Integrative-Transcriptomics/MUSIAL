package datastructure;

import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NavigableSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;
import utility.VariablePositionTableComparator;

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
public final class VariablePositionsTable {

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
   * sample the character `I` is appended: I.e. if an insertion at position 10 of length 2 was analyzed, the keys 10I
   * and 10II will be used. The third key set is sorted by the {@link VariablePositionTableComparator}. This system
   * allows for an easy insertion of `-` gap symbols in those samples in which no insertion occurred.
   */
  private final ConcurrentHashMap<String, ConcurrentHashMap<String, ConcurrentSkipListMap<String, VariantContent>>>
      table;
  /**
   * Thread safe structure mapping reference locations to a list of positions at which a variant was called.
   */
  private final ConcurrentHashMap<String, ConcurrentSkipListSet<String>> variantPositions;
  /**
   * Thread safe structure mapping reference locations to a list of {@link PositionStatistics} instances.
   */
  private final ConcurrentHashMap<String, ConcurrentHashMap<String, PositionStatistics>> positionStatistics;
  /**
   * Thread safe structure mapping reference locations to a list of {@link SampleStatistics} instances.
   */
  private final ConcurrentHashMap<String, ConcurrentHashMap<String, SampleStatistics>> sampleStatistics;

  /**
   * Constructor of {@link VariablePositionsTable}.
   */
  public VariablePositionsTable() {
    this.table = new ConcurrentHashMap<>();
    this.variantPositions = new ConcurrentHashMap<>();
    this.positionStatistics = new ConcurrentHashMap<>();
    this.sampleStatistics = new ConcurrentHashMap<>();
  }

  /**
   * Inserts a new {@link VariantContent} instance into the table structure.
   *
   * @param referenceAnalysisId {@link String} the name of the reference feature to access.
   * @param genomePosition      {@link String} representing the analyzed reference feature.
   * @param sampleName          {@link String} representing the name of the sample.
   * @param variablePosition    {@link VariantContent} to insert.
   */
  public void putVariablePosition(String referenceAnalysisId, String sampleName, String genomePosition,
                                  VariantContent variablePosition) {
    // Insert new ConcurrentHashMaps for reference feature if none is present.
    if (!this.table.containsKey(referenceAnalysisId)) {
      this.table.put(referenceAnalysisId, new ConcurrentHashMap<>());
    }
    if (!this.positionStatistics.containsKey(referenceAnalysisId)) {
      this.positionStatistics.put(referenceAnalysisId, new ConcurrentHashMap<>());
    }
    if (!this.sampleStatistics.containsKey(referenceAnalysisId)) {
      this.sampleStatistics.put(referenceAnalysisId, new ConcurrentHashMap<>());
    }
    // Insert new ConcurrentSkipListMap for sample name if none is present.
    if (!this.table.get(referenceAnalysisId).containsKey(sampleName)) {
      this.table.get(referenceAnalysisId).put(sampleName,
          new ConcurrentSkipListMap<>(new VariablePositionTableComparator()));
    }
    // Insert new ConcurrentHashMap for sample statistics if none is present.
    if (!this.sampleStatistics.get(referenceAnalysisId).containsKey(sampleName)) {
      SampleStatistics sampleStatistics = new SampleStatistics();
      this.sampleStatistics.get(referenceAnalysisId).put(sampleName, sampleStatistics);
    }
    // Insert new ConcurrentHashMap for position statistics if none is present.
    if (!this.positionStatistics.get(referenceAnalysisId).containsKey(genomePosition)) {
      PositionStatistics positionStatistics = new PositionStatistics();
      if (genomePosition.contains("I")) {
        positionStatistics.isInsertion = true;
        this.sampleStatistics.get(referenceAnalysisId).get(sampleName).insertions += 1;
      }
      this.positionStatistics.get(referenceAnalysisId).put(genomePosition, positionStatistics);
    }
    // Insert new ConcurrentSkipListSet to store variant positions if none is present.
    if (!this.variantPositions.containsKey(referenceAnalysisId)) {
      this.variantPositions
          .put(referenceAnalysisId, new ConcurrentSkipListSet<>(new VariablePositionTableComparator()));
    }
    // Insert the variable position information into the table.
    this.table.get(referenceAnalysisId).get(sampleName).put(genomePosition, variablePosition);
    // If the content of the passed variable position is no reference symbol and it does not represent the reference
    // sequence directly, the position is added as variant position. Thereby, the respective statistics entries are 
    // updated.
    boolean isHetCall = false;
    char content;
    for (int i = 0; i < variablePosition.numberOfAlleles; i++) {
      //noinspection ConstantConditions
      content = variablePosition.getContentOf(i);
      if (content != VariantContent.REFERENCE_REJECTED && content != VariantContent.DELETION_REJECTED &&
          content != VariantContent.ALT_A_REJECTED && content != VariantContent.ALT_C_REJECTED &&
          content != VariantContent.ALT_G_REJECTED && content != VariantContent.ALT_T_REJECTED &&
          content != VariantContent.REFERENCE && content != VariantContent.NO_CALL && !sampleName.equals(
          "Reference")) {
        this.variantPositions.get(referenceAnalysisId).add(genomePosition);
        if (i > 0) {
          isHetCall = true;
        }
      }
      switch (content) {
        case VariantContent.REFERENCE -> {
          this.positionStatistics.get(referenceAnalysisId).get(genomePosition).referenceCalls += 1;
          this.sampleStatistics.get(referenceAnalysisId).get(sampleName).referenceCalls += 1;
        }
        case VariantContent.DELETION -> {
          this.positionStatistics.get(referenceAnalysisId).get(genomePosition).deletions += 1;
          this.sampleStatistics.get(referenceAnalysisId).get(sampleName).deletions += 1;
        }
        case VariantContent.NO_CALL -> {
          this.positionStatistics.get(referenceAnalysisId).get(genomePosition).noCalls += 1;
          this.sampleStatistics.get(referenceAnalysisId).get(sampleName).noCalls += 1;
        }
        case VariantContent.ALT_A -> {
          this.positionStatistics.get(referenceAnalysisId).get(genomePosition).altACalls += 1;
          this.sampleStatistics.get(referenceAnalysisId).get(sampleName).altACalls += 1;
        }
        case VariantContent.ALT_C -> {
          this.positionStatistics.get(referenceAnalysisId).get(genomePosition).altCCalls += 1;
          this.sampleStatistics.get(referenceAnalysisId).get(sampleName).altCCalls += 1;
        }
        case VariantContent.ALT_G -> {
          this.positionStatistics.get(referenceAnalysisId).get(genomePosition).altGCalls += 1;
          this.sampleStatistics.get(referenceAnalysisId).get(sampleName).altGCalls += 1;
        }
        case VariantContent.ALT_T -> {
          this.positionStatistics.get(referenceAnalysisId).get(genomePosition).altTCalls += 1;
          this.sampleStatistics.get(referenceAnalysisId).get(sampleName).altTCalls += 1;
        }
        default -> {
          this.positionStatistics.get(referenceAnalysisId).get(genomePosition).rejectedCalls += 1;
          this.sampleStatistics.get(referenceAnalysisId).get(sampleName).rejectedCalls += 1;
        }
      }
    }
    if (isHetCall) {
      this.positionStatistics.get(referenceAnalysisId).get(genomePosition).hetCalls += 1;
      this.sampleStatistics.get(referenceAnalysisId).get(sampleName).hetCalls += 1;
    } else {
      this.positionStatistics.get(referenceAnalysisId).get(genomePosition).homCalls += 1;
      this.sampleStatistics.get(referenceAnalysisId).get(sampleName).homCalls += 1;
    }
  }

  /**
   * Inserts empty {@link ConcurrentHashMap} and {@link ConcurrentSkipListMap} structures into the table if none are
   * present for a specified reference analysis entry and sample name.
   * <p>
   * This method is intended to be used for samples with only confident reference calls as these will not be
   * considered for the output files otherwise.
   *
   * @param referenceAnalysisId {@link String} the name of the reference feature to access.
   * @param sampleName          {@link String} representing the name of the sample.
   */
  public void putSampleInReference(String referenceAnalysisId, String sampleName) {
    // Insert new ConcurrentHashMap for reference feature if none is present.
    if (!this.table.containsKey(referenceAnalysisId)) {
      this.table.put(referenceAnalysisId, new ConcurrentHashMap<>());
    }
    // Insert new ConcurrentSkipListMap for sample name if none is present.
    if (!this.table.get(referenceAnalysisId).containsKey(sampleName)) {
      this.table.get(referenceAnalysisId).put(sampleName,
          new ConcurrentSkipListMap<>(new VariablePositionTableComparator()));
    }
  }

  /**
   * Returns whether the table contains an entry for a specified reference location and sample name.
   *
   * @param referenceAnalysisId {@link String} the name of the reference feature to access.
   * @param sampleName          {@link String} representing the name of the sample.
   * @return Whether the table contains an entry for the passed referenceAnalysisId and sampleName.
   */
  public boolean referenceHasSample(String referenceAnalysisId, String sampleName) {
    if (this.table.containsKey(referenceAnalysisId)) {
      return this.table.get(referenceAnalysisId).containsKey(sampleName);
    } else {
      return false;
    }
  }

  /**
   * Returns a {@link VariantContent} stored in the table.
   * <p>
   * May return null if no entry exists.
   *
   * @param referenceAnalysisId {@link String} representing the reference analysis identifier / feature identifier to
   *                            access.
   * @param sampleName          {@link String} representing the sample identifier to access.
   * @param position            {@link String} representing the position to access.
   * @return {@link VariantContent} from the table.
   */
  public VariantContent getVariablePosition(String referenceAnalysisId, String sampleName, String position) {
    return this.table.get(referenceAnalysisId).get(sampleName).get(position);
  }

  /**
   * Returns all stored sample names for one reference location name.
   *
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @return {@link HashSet<String>} containing all sample names of the specified reference location.
   */
  public HashSet<String> getSampleNamesOf(String referenceAnalysisId) {
    HashSet<String> sampleNames = new HashSet<>();
    Iterator<String> keyIterator = this.table.get(referenceAnalysisId).keys().asIterator();
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
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @param sampleName          {@link String} the name of the sample.
   * @return {@link SampleStatistics} instance.
   */
  public SampleStatistics getSampleStatisticsOf(String referenceAnalysisId, String sampleName) {
    return this.sampleStatistics.get(referenceAnalysisId).get(sampleName);
  }

  /**
   * Returns the {@link PositionStatistics} instance stored for the specified reference location and position.
   *
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @param position            {@link String} the identifier of the position.
   * @return {@link PositionStatistics} instance.
   */
  public PositionStatistics getPositionStatisticsOf(String referenceAnalysisId, String position) {
    return this.positionStatistics.get(referenceAnalysisId).get(position);
  }

  /**
   * Returns a navigable set of all (including non-variant) positions of one reference feature.
   *
   * @param referenceAnalysisId {@link String} the reference feature to access.
   * @return A {@link TreeSet<String>} containing all positions of the specified reference feature.
   */
  public TreeSet<String> getAllPositions(String referenceAnalysisId) {
    TreeSet<String> allPositions = new TreeSet<>(new VariablePositionTableComparator());
    NavigableSet<String> referencePositions =
        this.table.get(referenceAnalysisId).get("Reference").navigableKeySet();
    allPositions.addAll(referencePositions);
    ConcurrentSkipListSet<String> variantPositions = this.getVariantPositionsSet(referenceAnalysisId);
    allPositions.addAll(variantPositions);
    return allPositions;
  }

  /**
   * Returns an iterator over all variant positions of a reference feature.
   *
   * @param referenceAnalysisId {@link String} the reference feature to access.
   * @return A {@link Iterator<String>} containing all variant positions of the specified reference feature.
   */
  public Iterator<String> getVariantPositions(String referenceAnalysisId) {
    if (this.variantPositions.containsKey(referenceAnalysisId)) {
      return this.variantPositions.get(referenceAnalysisId).iterator();
    } else {
      return Collections.emptyIterator();
    }
  }

  /**
   * Returns a navigable set of all variant positions of one reference feature.
   *
   * @param referenceAnalysisId {@link String} the reference feature to access.
   * @return A {@link ConcurrentSkipListSet<String>} containing all variant positions of the specified reference feature.
   */
  public ConcurrentSkipListSet<String> getVariantPositionsSet(String referenceAnalysisId) {
    if (this.variantPositions.containsKey(referenceAnalysisId)) {
      return this.variantPositions.get(referenceAnalysisId);
    } else {
      return new ConcurrentSkipListSet<>();
    }
  }

}