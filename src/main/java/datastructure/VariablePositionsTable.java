package datastructure;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NavigableSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;
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
   * Thread safe structure mapping reference locations to a list of positions, represented as for the table property.
   */
  private final ConcurrentHashMap<String, ConcurrentSkipListSet<String>> variantPositions;

  /**
   * Constructor of {@link VariablePositionsTable}.
   */
  public VariablePositionsTable() {
    this.table = new ConcurrentHashMap<>();
    this.variantPositions = new ConcurrentHashMap<>();
  }

  /**
   * Inserts a new {@link VariablePosition} instance into the table structure.
   *
   * @param genomePosition   {@link String} representing the analyzed reference feature.
   * @param sampleName       {@link String} representing the name of the sample.
   * @param variablePosition {@link VariablePosition} to insert.
   */
  public void putVariablePosition(String referenceAnalysisId, String sampleName, String genomePosition,
                                  VariablePosition variablePosition) {
    // Insert new ConcurrentHashMap for reference feature if none is present.
    if (!this.table.containsKey(referenceAnalysisId)) {
      this.table.put(referenceAnalysisId, new ConcurrentHashMap<>());
    }
    // Insert new ConcurrentSkipListMap for sample name if none is present.
    if (!this.table.get(referenceAnalysisId).containsKey(sampleName)) {
      this.table.get(referenceAnalysisId).put(sampleName,
          new ConcurrentSkipListMap<>(new VariablePositionTableComparator()));
    }
    // If a VariablePosition entry for the specified reference feature, sample and position already exists, a new
    // allele entry is added.
    if (this.table.get(referenceAnalysisId).get(sampleName).containsKey(genomePosition)) {
      this.table.get(referenceAnalysisId).get(sampleName).get(genomePosition).addAllele(
          variablePosition.content.get(0),
          variablePosition.annotation.get(0),
          variablePosition.frequency.get(0)
      );
      // Else a totally new VariablePosition is inserted.
    } else {
      this.table.get(referenceAnalysisId).get(sampleName).put(genomePosition, variablePosition);
    }
    // Insert new ConcurrentSkipListSet to store variant positions if none is present.
    if (!this.variantPositions.containsKey(referenceAnalysisId)) {
      this.variantPositions
          .put(referenceAnalysisId, new ConcurrentSkipListSet<>(new VariablePositionTableComparator()));
    }
    // If the content of the passed variable position is no reference symbol and it does not represent the reference
    // sequence directly, the position is added as variant position.
    char vPContent = variablePosition.content.get(0);
    if (vPContent != VariablePosition.REFERENCE && vPContent != VariablePosition.REFERENCE_LOW_COV &&
        vPContent != VariablePosition.REFERENCE_LOW_QUAL && !sampleName.equals("Reference")) {
      this.variantPositions.get(referenceAnalysisId).add(genomePosition);
    }
  }

  /**
   * Returns all stored reference location names as list.
   *
   * @return {@link List<String>} containing all reference location names.
   */
  public List<String> getReferenceAnalysisIds() {
    ArrayList<String> referenceIds = new ArrayList<>();
    Iterator<String> keyIterator = this.table.keys().asIterator();
    while (keyIterator.hasNext()) {
      referenceIds.add(keyIterator.next());
    }
    return referenceIds;
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
   * @return {@link String} representing the concatenation of all variant position contents.
   */
  public String getSampleVariantsSequence(String referenceAnalysisId, String sampleName) {
    StringBuilder variantSequenceBuilder = new StringBuilder();
    // Append `.fasta` header to variant sequence builder.
    variantSequenceBuilder.append("> ").append(sampleName).append(IO.LINE_SEPARATOR);
    // Append base symbol of each single position.
    int lineBreak = 0;
    for (String variantPosition : this.variantPositions.get(referenceAnalysisId)) {
      VariablePosition variablePosition = this.table.get(referenceAnalysisId).get(sampleName).get(variantPosition);
      if (variablePosition == null) {
        if (variantPosition.contains("I")) {
          // CASE: The position represents an insertion in any sample.
          variantSequenceBuilder.append("-");
        } else {
          // CASE: The position in the sample equals the reference sequence (and is not contained in the variant
          // positions table).
          variantSequenceBuilder.append(".");
        }
      } else {
        // Choose the allele content with the highest frequency.
        int maxIndex = -1;
        double maxValue = 0.0;
        for (int i = 0; i < variablePosition.frequency.size(); i++) {
          if (variablePosition.frequency.get(i) > maxValue) {
            maxValue = variablePosition.frequency.get(i);
            maxIndex = i;
          }
        }
        variantSequenceBuilder.append(variablePosition.content.get(maxIndex));
      }
      // Add a line break after 80 symbols to match maximal line length of the `.fasta` format.
      lineBreak += 1;
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
   * @return {@link String} representing the concatenation of all position contents.
   */
  public String getSampleFullSequence(String referenceAnalysisId, String sampleName) {
    StringBuilder fullSequenceBuilder = new StringBuilder();
    // Append `.fasta` header to full sequence builder.
    fullSequenceBuilder.append("> ").append(sampleName).append(IO.LINE_SEPARATOR);
    // Append base symbol of each single position.
    int lineBreak = 0;
    for (String position : this.getAllPositions(referenceAnalysisId)) {
      VariablePosition variablePosition = this.table.get(referenceAnalysisId).get(sampleName).get(position);
      if (variablePosition == null) {
        if (position.contains("I")) {
          // CASE: The position represents an insertion in any sample.
          fullSequenceBuilder.append("-");
        } else {
          // CASE: The position in the sample equals the reference sequence (and is not contained in the variant
          // positions table).
          fullSequenceBuilder.append(".");
        }
      } else {
        // Choose the allele content with the highest frequency.
        int maxIndex = -1;
        double maxValue = 0.0;
        for (int i = 0; i < variablePosition.frequency.size(); i++) {
          if (variablePosition.frequency.get(i) > maxValue) {
            maxValue = variablePosition.frequency.get(i);
            maxIndex = i;
          }
        }
        fullSequenceBuilder.append(variablePosition.content.get(maxIndex));
      }
      // Add a line break after 80 symbols to match maximal line length of the `.fasta` format.
      lineBreak += 1;
      if (lineBreak == 80) {
        fullSequenceBuilder.append(IO.LINE_SEPARATOR);
        lineBreak = 0;
      }
    }
    return fullSequenceBuilder.append(IO.LINE_SEPARATOR).toString();
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
   * @param addReference        {@link Boolean} whether to include reference information.
   * @return {@link String} representing the concatenation of all sample contents on the specified position.
   */
  public String getPositionContentTabDelimited(String referenceAnalysisId, String position, int shift,
                                               boolean addReference) {
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
    // If reference information should be added, the content of the reference is accessed and added at the first
    // position of the tab delimited format. Note: For this the reference information has to be stored in the
    // variable positions table with the key "Reference" as sample name.
    if (addReference) {
      VariablePosition referenceVariablePosition =
          this.table.get(referenceAnalysisId).get("Reference").get(position);
      if (referenceVariablePosition == null) {
        positionStringBuilder.append("-");
      } else {
        positionStringBuilder.append(referenceVariablePosition.content.get(0));
      }
      positionStringBuilder.append("\t");
    }
    for (String sampleName : sampleNames) {
      VariablePosition sampleVariablePosition =
          this.table.get(referenceAnalysisId).get(sampleName).get(position);
      if (sampleVariablePosition == null) {
        if (position.contains("I")) {
          // CASE: The position represents an insertion in any sample.
          positionStringBuilder.append("-");
        } else {
          // CASE: The position in the sample equals the reference sequence (and is not contained in the variant
          // positions table).
          positionStringBuilder.append(".");
        }
      } else {
        // If multiple alleles are present, each allele content is separated by the "|" symbol.
        StringBuilder allelesStringBuilder = new StringBuilder();
        for (Character alleleBase : sampleVariablePosition.content) {
          allelesStringBuilder.append(alleleBase);
          allelesStringBuilder.append("|");
        }
        String allelesString = allelesStringBuilder.toString();
        positionStringBuilder.append(allelesString, 0, allelesString.length() - 1);
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
    for (String sampleName : sampleNames) {
      VariablePosition sampleVariablePosition =
          this.table.get(referenceAnalysisId).get(sampleName).get(position);
      if (sampleVariablePosition == null) {
        if (position.contains("I")) {
          // CASE: The position represents an insertion in any sample.
          positionStringBuilder.append("-");
        } else {
          // CASE: The position in the sample equals the reference sequence (and is not contained in the variant
          // positions table).
          positionStringBuilder.append(".");
        }
      } else {
        // If multiple alleles are present, each allele content is separated by the "|" symbol.
        StringBuilder allelesStringBuilder = new StringBuilder();
        for (Double alleleFrequency : sampleVariablePosition.frequency) {
          allelesStringBuilder.append("COV=").append(sampleVariablePosition.coverage).append(";").append("QUAL=")
              .append(sampleVariablePosition.quality).append(";").append("FREQ=").append(alleleFrequency).append("|");
        }
        String allelesString = allelesStringBuilder.toString();
        positionStringBuilder.append(allelesString, 0, allelesString.length() - 1);
      }
      positionStringBuilder.append("\t");
    }
    String positionVariationsTabDelimited = positionStringBuilder.toString();
    return positionVariationsTabDelimited.substring(0, positionVariationsTabDelimited.length() - 1);
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
        // Currently only the allele with the highest frequency is taken into account.
        int maxIndex = -1;
        double maxValue = -1;
        for (int i = 0; i < sampleVariablePosition.frequency.size(); i++) {
          if (sampleVariablePosition.frequency.get(i) > maxValue) {
            maxValue = sampleVariablePosition.frequency.get(i);
            maxIndex = i;
          }
        }
        char sampleVariablePositionContent = sampleVariablePosition.content.get(maxIndex);
        switch (sampleVariablePositionContent) {
          case VariablePosition.REFERENCE_LOW_QUAL:
            refLowQual += 1;
            break;
          case VariablePosition.REFERENCE_LOW_COV:
            refLowCov += 1;
            break;
          case VariablePosition.ALT_A:
            altA += 1;
            break;
          case VariablePosition.ALT_A_LOW_QUAL:
            altALowQual += 1;
            break;
          case VariablePosition.ALT_A_LOW_COV:
            altALowCov += 1;
            break;
          case VariablePosition.ALT_A_LOW_FREQ:
            altALowFreq += 1;
            break;
          case VariablePosition.ALT_C:
            altC += 1;
            break;
          case VariablePosition.ALT_C_LOW_QUAL:
            altCLowQual += 1;
            break;
          case VariablePosition.ALT_C_LOW_COV:
            altCLowCov += 1;
            break;
          case VariablePosition.ALT_C_LOW_FREQ:
            altCLowFreq += 1;
            break;
          case VariablePosition.ALT_G:
            altG += 1;
            break;
          case VariablePosition.ALT_G_LOW_QUAL:
            altGLowQual += 1;
            break;
          case VariablePosition.ALT_G_LOW_COV:
            altGLowCov += 1;
            break;
          case VariablePosition.ALT_G_LOW_FREQ:
            altGLowFreq += 1;
            break;
          case VariablePosition.ALT_T:
            altT += 1;
            break;
          case VariablePosition.ALT_T_LOW_QUAL:
            altTLowQual += 1;
            break;
          case VariablePosition.ALT_T_LOW_COV:
            altTLowCov += 1;
            break;
          case VariablePosition.ALT_T_LOW_FREQ:
            altTLowFreq += 1;
            break;
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
        .append("\t").append(refLowQual).append("\t").append(refLowCov)
        .append("\t").append(del)
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
   * Returns a navigable set of all (including non-variant) positions of one reference feature.
   *
   * @param referenceAnalysisId {@link String} the reference feature to access.
   * @return A {@link TreeSet<String>} containing all positions of the specified reference feature.
   */
  public TreeSet<String> getAllPositions(String referenceAnalysisId) {
    TreeSet<String> allPositions = new TreeSet<>(new VariablePositionTableComparator());
    NavigableSet<String> referencePositions =
        this.table.get(referenceAnalysisId).get("Reference").navigableKeySet();
    ConcurrentSkipListSet<String> variantPositions = this.variantPositions.get(referenceAnalysisId);
    allPositions.addAll(referencePositions);
    allPositions.addAll(variantPositions);
    return allPositions;
  }

  /**
   * Returns a navigable set of all variant positions of one reference feature.
   *
   * @param referenceAnalysisId {@link String} the reference feature to access.
   * @return A {@link TreeSet<String>} containing all positions of the specified reference feature.
   */
  public Iterator<String> getVariantPositions(String referenceAnalysisId) {
    return this.variantPositions.get(referenceAnalysisId).iterator();
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
   * @param addReference        {@link Boolean} whether a field "Reference" should be added.
   * @return {@link String} representing a tab delimited header line for output files.
   */
  public String getSampleNamesHeaderTabDelimited(String referenceAnalysisId, boolean addReference) {
    StringBuilder positionsStringBuilder = new StringBuilder();
    List<String> sampleNames = this.getSampleNamesOf(referenceAnalysisId);
    positionsStringBuilder.append("Position");
    positionsStringBuilder.append("\t");
    if (addReference) {
      positionsStringBuilder.append("Reference");
      positionsStringBuilder.append("\t");
    }
    for (String sampleName : sampleNames) {
      positionsStringBuilder.append(sampleName);
      positionsStringBuilder.append("\t");
    }
    String headerString = positionsStringBuilder.toString();
    return headerString.substring(0, headerString.length() - 1);
  }

  /**
   * Returns a header line for tab delimited output in which columns are possible variant content identifiers and
   * each row starts with a position.
   *
   * @return {@link String} representing a tab delimited header line for output files.
   */
  public String getPositionStatisticsHeaderTabDelimited() {
    StringBuilder headerStringBuilder = new StringBuilder();
    headerStringBuilder.append("Position\tRef\tRef_LowQuality\tRef_LowCoverage\tDeletion\t");
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

}