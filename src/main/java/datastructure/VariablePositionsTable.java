package datastructure;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;
import runnables.SampleAnalyserRunnable;
import utility.IO;
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
   * {@link VariablePosition} instances.
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
  private final ConcurrentHashMap<String, ConcurrentHashMap<String, ConcurrentSkipListMap<String, VariablePosition>>>
      table;
  /**
   * Thread safe structure mapping reference locations to a list of positions at which a variant was called.
   */
  private final ConcurrentHashMap<String, ConcurrentSkipListSet<String>> variantPositions;
  /**
   * Thread safe structure mapping reference locations to a list of positions at which an insertion is present.
   */
  private final ConcurrentHashMap<String, ConcurrentSkipListSet<String>> insertionPositions;

  /**
   * Constructor of {@link VariablePositionsTable}.
   */
  public VariablePositionsTable() {
    this.table = new ConcurrentHashMap<>();
    this.variantPositions = new ConcurrentHashMap<>();
    this.insertionPositions = new ConcurrentHashMap<>();
  }

  /**
   * Inserts a new {@link VariablePosition} instance into the table structure.
   *
   * @param referenceAnalysisId {@link String} the name of the reference feature to access.
   * @param genomePosition      {@link String} representing the analyzed reference feature.
   * @param sampleName          {@link String} representing the name of the sample.
   * @param variablePosition    {@link VariablePosition} to insert.
   * @param alleleIndex         {@link Integer} representing the index of a possible alternate allele.
   */
  public void putVariablePosition(String referenceAnalysisId, String sampleName, String genomePosition,
                                  VariablePosition variablePosition, int alleleIndex) {
    // If the alleleIndex is not zero, the sampleName is adjusted in order to represent all 2nd, 3rd, etc.
    // alternative alleles.
    if (alleleIndex != 0) {
      sampleName = sampleName + "~" + alleleIndex + "Alternative";
    }
    // Insert new ConcurrentHashMap for reference feature if none is present.
    if (!this.table.containsKey(referenceAnalysisId)) {
      this.table.put(referenceAnalysisId, new ConcurrentHashMap<>());
    }
    // Insert new ConcurrentSkipListMap for sample name if none is present.
    if (!this.table.get(referenceAnalysisId).containsKey(sampleName)) {
      this.table.get(referenceAnalysisId).put(sampleName,
          new ConcurrentSkipListMap<>(new VariablePositionTableComparator()));
    }
    // Insert the variable position information into the table.
    this.table.get(referenceAnalysisId).get(sampleName).put(genomePosition, variablePosition);
    // Insert new ConcurrentSkipListSet to store variant positions if none is present.
    if (!this.variantPositions.containsKey(referenceAnalysisId)) {
      this.variantPositions
          .put(referenceAnalysisId, new ConcurrentSkipListSet<>(new VariablePositionTableComparator()));
    }
    // Insert new ConcurrentSkipListSet to store insertion positions if none is present.
    if (genomePosition.contains("I")) {
      this.insertionPositions.put(referenceAnalysisId,
          new ConcurrentSkipListSet<>(new VariablePositionTableComparator()));
      this.insertionPositions.get(referenceAnalysisId).add(genomePosition);
    }
    // If the content of the passed variable position is no reference symbol and it does not represent the reference
    // sequence directly, the position is added as variant position.
    char vPContent = variablePosition.content;
    if (vPContent != VariablePosition.REFERENCE && vPContent != VariablePosition.REFERENCE_DISCARDED &&
        vPContent != VariablePosition.NO_CALL && !sampleName.equals(
        "Reference")) {
      this.variantPositions.get(referenceAnalysisId).add(genomePosition);
    }
  }

  /**
   * Inserts empty {@link ConcurrentHashMap} and {@link ConcurrentSkipListMap} structures into the table if none are
   * present for a specified reference analysis entry and sample name.
   * <p>
   * This method is intended to be used for samples with only confident reference calls are these will not be
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
   * Returns a {@link VariablePosition} stored in the table.
   * <p>
   * May return null of no entry exists.
   *
   * @param referenceAnalysisId {@link String} representing the reference analysis identifier / feature identifier to
   *                            access.
   * @param sampleName          {@link String} representing the sample identifier to access.
   * @param position            {@link String} representing the position to access.
   * @return {@link VariablePosition} from the table.
   */
  public VariablePosition getVariablePosition(String referenceAnalysisId, String sampleName, String position) {
    return this.table.get(referenceAnalysisId).get(sampleName).get(position);
  }

  /**
   * Returns all stored feature location names as list.
   *
   * @return {@link List<String>} containing all feature location names.
   */
  public List<String> getFeatureIdentifiers() {
    ArrayList<String> featureIds = new ArrayList<>();
    Iterator<String> keyIterator = this.table.keys().asIterator();
    while (keyIterator.hasNext()) {
      featureIds.add(keyIterator.next());
    }
    return featureIds;
  }

  /**
   * Returns all stored sample names for one reference location name.
   *
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @return {@link List<String>} containing all sample names of the specified reference location.
   */
  public List<String> getSampleNamesOf(String referenceAnalysisId) {
    ArrayList<String> sampleNames = new ArrayList<>();
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
   * Returns the concatenation of all variant position contents for one sample.
   * <p>
   * The contents are appended in ascending order with respect to the position value. Currently, if multiple alleles
   * for one position are present, the one with the highest frequency is used.
   *
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @param sampleName          {@link String} the name of the sample.
   * @param fastaFormat         {@link Boolean} whether the sequence should be returned in Fasta format.
   * @return {@link String} representing the concatenation of all variant position contents.
   */
  public String getSampleVariantsSequence(String referenceAnalysisId, String sampleName, boolean fastaFormat) {
    StringBuilder variantSequenceBuilder = new StringBuilder();
    // Append `.fasta` header to variant sequence builder.
    if (fastaFormat) {
      variantSequenceBuilder.append("> ").append(sampleName).append(IO.LINE_SEPARATOR);
    }
    // Append base symbol of each single position.
    int lineBreak = 0;
    for (String variantPosition : this.getVariantPositionsSet(referenceAnalysisId)) {
      VariablePosition variablePosition = this.table.get(referenceAnalysisId).get(sampleName).get(variantPosition);
      if (variablePosition == null) {
        if (variantPosition.contains("I")) {
          // CASE: The position represents an insertion in any sample.
          variantSequenceBuilder.append("-");
        } else if (sampleName.endsWith("Alternative")) {
          // CASE: The sample is a dummy to represent an alternative allele of the original sample, but the position
          // is a hom. call.
          variantSequenceBuilder.append("=");
        } else {
          // CASE: The position in the sample equals the reference sequence (and is not contained in the variant
          // positions table).
          variantSequenceBuilder.append(".");
        }
      } else {
        // Append the content.
        variantSequenceBuilder.append(variablePosition.content);
      }
      // Add a line break after 80 symbols to match maximal line length of the `.fasta` format.
      if (fastaFormat) {
        lineBreak += 1;
      }
      if (lineBreak == 80) {
        variantSequenceBuilder.append(IO.LINE_SEPARATOR);
        lineBreak = 0;
      }
    }
    return variantSequenceBuilder.append(IO.LINE_SEPARATOR).toString();
  }

  /**
   * Returns the concatenation of all position contents for one sample.
   * <p>
   * The contents are appended in ascending order with respect to the position value. Currently, if multiple alleles
   * for one position are present, the one with the highest frequency is used.
   *
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @param sampleName          {@link String} the name of the sample.
   * @param fastaFormat         {@link Boolean} whether the sequence should be returned in Fasta format.
   * @return {@link String} representing the concatenation of all position contents.
   */
  public String getSampleFullSequence(String referenceAnalysisId, String sampleName, boolean fastaFormat) {
    StringBuilder fullSequenceBuilder = new StringBuilder();
    // Append `.fasta` header to full sequence builder.
    if (fastaFormat) {
      fullSequenceBuilder.append("> ").append(sampleName).append(IO.LINE_SEPARATOR);
    }
    // Append base symbol of each single position.
    int lineBreak = 0;
    for (String position : this.getAllPositions(referenceAnalysisId)) {
      VariablePosition variablePosition = this.table.get(referenceAnalysisId).get(sampleName).get(position);
      if (variablePosition == null) {
        if (position.contains("I")) {
          // CASE: The position represents an insertion in any sample.
          fullSequenceBuilder.append("-");
        } else if (sampleName.endsWith("Alternative")) {
          // CASE: The sample is a dummy to represent an alternative allele of the original sample, but the position
          // is a hom. call.
          fullSequenceBuilder.append("=");
        } else {
          // CASE: The position in the sample equals the reference sequence (and is not contained in the variant
          // positions table).
          fullSequenceBuilder.append(".");
        }
      } else {
        // Append the content.
        fullSequenceBuilder.append(variablePosition.content);
      }
      // Add a line break after 80 symbols to match maximal line length of the `.fasta` format.
      if (fastaFormat) {
        lineBreak += 1;
      }
      if (lineBreak == 80) {
        fullSequenceBuilder.append(IO.LINE_SEPARATOR);
        lineBreak = 0;
      }
    }
    if (fastaFormat) {
      fullSequenceBuilder.append(IO.LINE_SEPARATOR);
    }
    return fullSequenceBuilder.toString();
  }

  /**
   * Returns the concatenation of all sample contents for one position in tab delimited format of one reference
   * feature location.
   * <p>
   * The contents are appended in the order given by the sample names.
   *
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @param position            {@link String} the position to access.
   * @param shift               {@link Integer} representing a position shift, i.e. the sum of all insertions that occurred up to
   *                            and including the passed position.
   * @return {@link String} representing the concatenation of all sample contents on the specified position.
   */
  public String getPositionContentTabDelimited(String referenceAnalysisId, String position, int shift) {
    StringBuilder positionStringBuilder = new StringBuilder();
    List<String> sampleNames = this.getSampleNamesOf(referenceAnalysisId);
    // If the passed position contains an I it represents an insertion position and the integer value of the string
    // without any Is is parsed. Note: The shift has to be passed to the method, including the current position if it
    // is an insertion position.
    if (position.contains("I")) {
      positionStringBuilder.append(Integer.parseInt(position.substring(0, position.indexOf("I"))) + shift);
    } else {
      positionStringBuilder.append(Integer.parseInt(position) + shift);
    }
    positionStringBuilder.append("\t");
    // The reference is accessed and added at the first
    // position of the tab delimited format. Note: For this the reference information has to be stored in the
    // variable positions table with the key "Reference" as sample name.
    VariablePosition referenceVariablePosition =
        this.table.get(referenceAnalysisId).get("Reference").get(position);
    if (referenceVariablePosition == null) {
      positionStringBuilder.append("-");
    } else {
      positionStringBuilder.append(referenceVariablePosition.content);
    }
    positionStringBuilder.append("\t");
    for (String sampleName : sampleNames) {
      VariablePosition sampleVariablePosition =
          this.table.get(referenceAnalysisId).get(sampleName).get(position);
      if (sampleVariablePosition == null) {
        if (position.contains("I")) {
          // CASE: The position represents an insertion in any sample.
          positionStringBuilder.append("-");
        } else if (sampleName.endsWith("Alternative")) {
          // CASE: The sample is a dummy to represent an alternative allele of the original sample, but the position
          // is a hom. call.
          positionStringBuilder.append("=");
        } else {
          // CASE: The position in the sample equals the reference sequence (and is not contained in the variant
          // positions table).
          positionStringBuilder.append(".");
        }
      } else {
        // CASE: The position represents a variant.
        positionStringBuilder.append(sampleVariablePosition.content);
      }
      positionStringBuilder.append("\t");
    }
    String positionVariationsTabDelimited = positionStringBuilder.toString();
    return positionVariationsTabDelimited.substring(0, positionVariationsTabDelimited.length() - 1);
  }

  /**
   * Returns the concatenation of all sample annotations for one position in tab delimited format of one reference
   * feature location.
   * <p>
   * The annotations are appended in the order given by the sample names.
   *
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @param position            {@link String} the position to access.
   * @param shift               {@link Integer} representing a position shift, i.e. the sum of all insertions that occurred up to
   *                            and including the passed position.
   * @return {@link String} representing the concatenation of all sample annotations on the specified position.
   */
  public String getPositionAnnotationTabDelimited(String referenceAnalysisId, String position, int shift) {
    StringBuilder positionStringBuilder = new StringBuilder();
    List<String> sampleNames = this.getSampleNamesOf(referenceAnalysisId);
    // If the passed position contains an I it represents an insertion position and the integer value of the string
    // without any Is is parsed. Note: The shift has to be passed to the method, including the current position if it
    // is an insertion position.
    if (position.contains("I")) {
      positionStringBuilder.append(Integer.parseInt(position.substring(0, position.indexOf("I"))) + shift);
    } else {
      positionStringBuilder.append(Integer.parseInt(position) + shift);
    }
    positionStringBuilder.append("\t");
    // The reference is accessed and added at the first
    // position of the tab delimited format. Note: For this the reference information has to be stored in the
    // variable positions table with the key "Reference" as sample name.
    VariablePosition referenceVariablePosition =
        this.table.get(referenceAnalysisId).get("Reference").get(position);
    if (referenceVariablePosition != null) {
      positionStringBuilder.append(referenceVariablePosition.content);
    }
    positionStringBuilder.append("\t");
    for (String sampleName : sampleNames) {
      VariablePosition sampleVariablePosition =
          this.table.get(referenceAnalysisId).get(sampleName).get(position);
      if (sampleVariablePosition != null) {
        positionStringBuilder.append(sampleVariablePosition.annotation);
      }
      positionStringBuilder.append("\t");
    }
    String positionVariationsTabDelimited = positionStringBuilder.toString();
    return positionVariationsTabDelimited.substring(0, positionVariationsTabDelimited.length() - 1);
  }

  /**
   * Returns a header line for tab delimited output in which columns are possible variant content identifiers.
   *
   * @param firstColumnName {@link String} representing the name of the first column.
   * @return {@link String} representing a tab delimited header line for output files.
   */
  public String getStatisticsHeaderTabDelimited(String firstColumnName) {
    StringBuilder headerStringBuilder = new StringBuilder();
    headerStringBuilder.append(firstColumnName).append("\tRef\tRef_LowQuality\tRef_LowCoverage\tRef_LowFrequency" +
        "\tDeletion\tNo_Call\t");
    char[] bases = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < bases.length; i++) {
      char base = bases[i];
      headerStringBuilder.append("Alt_").append(base).append("\tAlt_").append(base).append("_LowQuality\tAlt_")
          .append(base).append("_LowCoverage\tAlt_").append(base).append("_LowFrequency");
      if (i != (bases.length - 1)) {
        headerStringBuilder.append("\t");
      }
    }
    return headerStringBuilder.toString();
  }

  /**
   * Computes and returns the number of each possible contents at one position for all samples of one reference
   * feature location.
   *
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @param position            {@link String} the position to access.
   * @param shift               {@link Integer} representing a position shift, i.e. the sum of all insertions that occurred up to
   *                            and including the passed position.
   * @return {@link String} representing the counts of each possible content at the specified position.
   */
  public String getPositionStatisticsTabDelimited(String referenceAnalysisId, String position, int shift) {
    // First, initialize a counter for each possible variant content. See the documentation of the VariablePosition
    // class for more information.
    int ref = 0;
    int refLowQual = 0;
    int refLowCov = 0;
    int refLowFreq = 0;
    int del = 0;
    int altA = 0;
    int altALowQual = 0;
    int altALowCov = 0;
    int altALowFreq = 0;
    int altC = 0;
    int altCLowQual = 0;
    int altCLowCov = 0;
    int altCLowFreq = 0;
    int altG = 0;
    int altGLowQual = 0;
    int altGLowCov = 0;
    int altGLowFreq = 0;
    int altT = 0;
    int altTLowQual = 0;
    int altTLowCov = 0;
    int altTLowFreq = 0;
    int noCall = 0;
    // Iterate over each sample of the specified reference feature location and position.
    List<String> sampleNames = this.getSampleNamesOf(referenceAnalysisId);
    for (String sampleName : sampleNames) {
      VariablePosition sampleVariablePosition =
          this.table.get(referenceAnalysisId).get(sampleName).get(position);
      if (sampleVariablePosition == null) {
        if (position.contains("I")) {
          // CASE: The position represents an insertion in any sample.
          del += 1;
        } else {
          // CASE: The position in the sample equals the reference sequence (and is not contained in the variant
          // positions table).
          ref += 1;
        }
      } else {
        // CASE: The position in the sample is an alternative base.
        char sampleVariablePositionContent = sampleVariablePosition.content;
        List<String> sampleVariablePositionAnnotations;
        String[] splitAnnotation =
            sampleVariablePosition.annotation.split(String.valueOf(SampleAnalyserRunnable.ANN_SEP));
        switch (sampleVariablePositionContent) {
          case VariablePosition.REFERENCE -> ref += 1;
          case VariablePosition.REFERENCE_DISCARDED -> {
            sampleVariablePositionAnnotations =
                Arrays.asList(splitAnnotation);
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_QUALITY)) {
              refLowQual += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_COVERAGE)) {
              refLowCov += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_FREQUENCY)) {
              refLowFreq += 1;
            }
          }
          case VariablePosition.ALT_A -> altA += 1;
          case VariablePosition.ALT_A_DISCARDED -> {
            sampleVariablePositionAnnotations =
                Arrays.asList(splitAnnotation);
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_QUALITY)) {
              altALowQual += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_COVERAGE)) {
              altALowCov += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_FREQUENCY)) {
              altALowFreq += 1;
            }
          }
          case VariablePosition.ALT_C -> altC += 1;
          case VariablePosition.ALT_C_DISCARDED -> {
            sampleVariablePositionAnnotations =
                Arrays.asList(splitAnnotation);
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_QUALITY)) {
              altCLowQual += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_COVERAGE)) {
              altCLowCov += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_FREQUENCY)) {
              altCLowFreq += 1;
            }
          }
          case VariablePosition.ALT_G -> altG += 1;
          case VariablePosition.ALT_G_DISCARDED -> {
            sampleVariablePositionAnnotations =
                Arrays.asList(splitAnnotation);
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_QUALITY)) {
              altGLowQual += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_COVERAGE)) {
              altGLowCov += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_FREQUENCY)) {
              altGLowFreq += 1;
            }
          }
          case VariablePosition.ALT_T -> altT += 1;
          case VariablePosition.ALT_T_DISCARDED -> {
            sampleVariablePositionAnnotations =
                Arrays.asList(splitAnnotation);
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_QUALITY)) {
              altTLowQual += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_COVERAGE)) {
              altTLowCov += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_FREQUENCY)) {
              altTLowFreq += 1;
            }
          }
          case VariablePosition.NO_CALL -> noCall += 1;
        }
      }
    }
    // Build a tab delimited string from the conducted counting.
    StringBuilder positionStringBuilder = new StringBuilder();
    if (position.contains("I")) {
      positionStringBuilder.append(Integer.parseInt(position.substring(0, position.indexOf("I"))) + shift);
    } else {
      positionStringBuilder.append(Integer.parseInt(position) + shift);
    }
    positionStringBuilder
        .append("\t").append(ref)
        .append("\t").append(refLowQual).append("\t").append(refLowCov).append("\t").append(refLowFreq)
        .append("\t").append(del).append("\t").append(noCall)
        .append("\t")
        .append(altA).append("\t").append(altALowQual).append("\t").append(altALowCov).append("\t").append(altALowFreq)
        .append("\t")
        .append(altC).append("\t").append(altCLowQual).append("\t").append(altCLowCov).append("\t").append(altCLowFreq)
        .append("\t")
        .append(altG).append("\t").append(altGLowQual).append("\t").append(altGLowCov).append("\t").append(altGLowFreq)
        .append("\t")
        .append(altT).append("\t").append(altTLowQual).append("\t").append(altTLowCov).append("\t").append(altTLowFreq);
    return positionStringBuilder.toString();
  }

  /**
   * Computes and returns the number of each possible contents of one sample for all positions of one reference
   * feature location.
   *
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @param sampleName          {@link String} the sample to access.
   * @return {@link String} representing the counts of each possible content at the specified position.
   */
  public String getSampleStatisticsTabDelimited(String referenceAnalysisId, String sampleName) {
    // First, initialize a counter for each possible variant content. See the documentation of the VariablePosition
    // class for more information.
    int ref = 0;
    int refLowQual = 0;
    int refLowCov = 0;
    int refLowFreq = 0;
    int del = 0;
    int altA = 0;
    int altALowQual = 0;
    int altALowCov = 0;
    int altALowFreq = 0;
    int altC = 0;
    int altCLowQual = 0;
    int altCLowCov = 0;
    int altCLowFreq = 0;
    int altG = 0;
    int altGLowQual = 0;
    int altGLowCov = 0;
    int altGLowFreq = 0;
    int altT = 0;
    int altTLowQual = 0;
    int altTLowCov = 0;
    int altTLowFreq = 0;
    int noCall = 0;
    // Iterate over each sample of the specified reference feature location and position.
    for (Iterator<String> it = this.getVariantPositions(referenceAnalysisId); it.hasNext(); ) {
      String variantPosition = it.next();
      VariablePosition sampleVariablePosition =
          this.table.get(referenceAnalysisId).get(sampleName).get(variantPosition);
      if (sampleVariablePosition == null) {
        if (variantPosition.contains("I")) {
          // CASE: The position represents an insertion in any other sample.
          del += 1;
        } else {
          // CASE: The position in the sample equals the reference sequence (and is not contained in the variant
          // positions table).
          ref += 1;
        }
      } else {
        // CASE: The position in the sample is an alternative base.
        char sampleVariablePositionContent = sampleVariablePosition.content;
        List<String> sampleVariablePositionAnnotations;
        var split = sampleVariablePosition.annotation.split(String.valueOf(SampleAnalyserRunnable.ANN_SEP));
        switch (sampleVariablePositionContent) {
          case VariablePosition.REFERENCE -> ref += 1;
          case VariablePosition.REFERENCE_DISCARDED -> {
            sampleVariablePositionAnnotations =
                Arrays.asList(split);
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_QUALITY)) {
              refLowQual += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_COVERAGE)) {
              refLowCov += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_FREQUENCY)) {
              refLowFreq += 1;
            }
          }
          case VariablePosition.ALT_A -> altA += 1;
          case VariablePosition.ALT_A_DISCARDED -> {
            sampleVariablePositionAnnotations =
                Arrays.asList(split);
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_QUALITY)) {
              altALowQual += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_COVERAGE)) {
              altALowCov += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_FREQUENCY)) {
              altALowFreq += 1;
            }
          }
          case VariablePosition.ALT_C -> altC += 1;
          case VariablePosition.ALT_C_DISCARDED -> {
            sampleVariablePositionAnnotations =
                Arrays.asList(split);
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_QUALITY)) {
              altCLowQual += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_COVERAGE)) {
              altCLowCov += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_FREQUENCY)) {
              altCLowFreq += 1;
            }
          }
          case VariablePosition.ALT_G -> altG += 1;
          case VariablePosition.ALT_G_DISCARDED -> {
            sampleVariablePositionAnnotations =
                Arrays.asList(split);
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_QUALITY)) {
              altGLowQual += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_COVERAGE)) {
              altGLowCov += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_FREQUENCY)) {
              altGLowFreq += 1;
            }
          }
          case VariablePosition.ALT_T -> altT += 1;
          case VariablePosition.ALT_T_DISCARDED -> {
            sampleVariablePositionAnnotations =
                Arrays.asList(split);
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_QUALITY)) {
              altTLowQual += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_COVERAGE)) {
              altTLowCov += 1;
            }
            if (sampleVariablePositionAnnotations.contains(SampleAnalyserRunnable.LOW_FREQUENCY)) {
              altTLowFreq += 1;
            }
          }
          case VariablePosition.NO_CALL -> noCall += 1;
        }
      }
    }
    // Build a tab delimited string from the conducted counting.
    return sampleName +
        "\t" + ref +
        "\t" + refLowQual + "\t" + refLowCov + "\t" + refLowFreq +
        "\t" + del + "\t" + noCall +
        "\t" +
        altA + "\t" + altALowQual + "\t" + altALowCov + "\t" + altALowFreq +
        "\t" +
        altC + "\t" + altCLowQual + "\t" + altCLowCov + "\t" + altCLowFreq +
        "\t" +
        altG + "\t" + altGLowQual + "\t" + altGLowCov + "\t" + altGLowFreq +
        "\t" +
        altT + "\t" + altTLowQual + "\t" + altTLowCov + "\t" + altTLowFreq;
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

  /**
   * Returns a navigable set of all insertion positions of one reference feature.
   *
   * @param referenceAnalysisId {@link String} the reference feature to access.
   * @return A {@link TreeSet<String>} containing all insertion positions of the specified reference feature.
   */
  public Iterator<String> getInsertionPositions(String referenceAnalysisId) {
    if (this.insertionPositions.containsKey(referenceAnalysisId)) {
      return this.insertionPositions.get(referenceAnalysisId).iterator();
    } else {
      return Collections.emptyIterator();
    }
  }

  /**
   * Returns a header line for tab delimited output in which columns are sample names and each row starts with a
   * position.
   * <p>
   * If addReference is set to true, the second field will be "Reference".
   * <p>
   * The output string will have the format: Position\tReference(optional)\tSAMPLE1_NAME\tSAMPLE2_NAME\t...
   *
   * @param referenceAnalysisId {@link String} the name of the reference location.
   * @return {@link String} representing a tab delimited header line for output files.
   */
  public String getSampleNamesHeaderTabDelimited(String referenceAnalysisId) {
    StringBuilder positionsStringBuilder = new StringBuilder();
    List<String> sampleNames = this.getSampleNamesOf(referenceAnalysisId);
    positionsStringBuilder.append("Position");
    positionsStringBuilder.append("\t");
    positionsStringBuilder.append("Reference");
    positionsStringBuilder.append("\t");
    for (String sampleName : sampleNames) {
      positionsStringBuilder.append(sampleName);
      positionsStringBuilder.append("\t");
    }
    String headerString = positionsStringBuilder.toString();
    return headerString.substring(0, headerString.length() - 1);
  }

}