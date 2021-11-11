package runnables;

import datastructure.FeatureAnalysisEntry;
import datastructure.VariablePosition;
import datastructure.VariablePositionsTable;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.ArrayList;
import java.util.List;
import me.tongfei.progressbar.ProgressBar;
import utility.Bio;

/**
 * Implementation of the {@link Runnable} interface to analyze a sample.
 * <p>
 * Runs a sample analysis in a single thread.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class SampleAnalyserRunnable implements Runnable {

  /**
   * Annotation tag to indicate low coverage.
   */
  public static final String LOW_COVERAGE = "LOW_COVERAGE";
  /**
   * Annotation tag to indicate low quality.
   */
  public static final String LOW_QUALITY = "LOW_QUALITY";
  /**
   * Annotation tag to indicate low frequency.
   */
  public static final String LOW_FREQUENCY = "LOW_FREQUENCY";
  /**
   * Annotation tag to indicate no call.
   */
  public static final String NO_CALL = "NO_CALL";
  /**
   * Annotation tag separator.
   */
  public static final char ANN_SEP = '@';
  /**
   * The internal name of the sample to analyze.
   */
  private final String sampleName;
  /**
   * An {@link FeatureAnalysisEntry} instance storing information about the reference with regard to which the
   * sample is analyzed.
   */
  private final FeatureAnalysisEntry featureAnalysisEntry;
  /**
   * An {@link VariantContext} iterator storing information about single entries (rows) of the `.vcf` file of the
   * respective sample.
   */
  private final CloseableIterator<VariantContext> variantInformation;
  /**
   * An {@link VariablePositionsTable} instance to which the parsed information about variants is passed.
   */
  private final VariablePositionsTable variablePositionsTable;
  /**
   * The minimum coverage to accept a variant specified by the user.
   */
  private final double minCoverage;
  /**
   * The minimum hom. call frequency to accept a variant specified by the user.
   */
  private final double minHomFrequency;
  /**
   * The minimum het. call frequency to accept a variant specified by the user.
   */
  private final double minHetFrequency;
  /**
   * The maximum het. call frequency to accept a variant specified by the user.
   */
  private final double maxHetFrequency;
  /**
   * The minimum frequency to accept a variant specified by the user.
   */
  private final double minQuality;
  /**
   * An {@link ProgressBar} instance to output user information during runtime.
   */
  private final ProgressBar progress;

  /**
   * Constructor of {@link SampleAnalyserRunnable}.
   *
   * @param sampleName             {@link String} the samples name.
   * @param featureAnalysisEntry   {@link FeatureAnalysisEntry} information about the reference to which respect
   *                               the sample is analyzed.
   * @param variantInformation     {@link CloseableIterator<VariantContext>} returned by
   *                               a {@link htsjdk.variant.vcf.VCFFileReader}
   *                               giving access to the single input `.vcf`
   *                               files entries.
   * @param variablePositionsTable {@link VariablePositionsTable} to store the processed information.
   * @param minCoverage            {@link Double} minimum coverage specified by the user.
   * @param minHomFrequency        {@link Double} minimum hom. call frequency (wrt. read support) specified by the
   *                               user.
   * @param minHetFrequency        {@link Double} minimum het. call frequency (wrt. read support) specified by the
   *                               user.
   * @param maxHetFrequency        {@link Double} maximum het. call frequency (wrt. read support) specified by the
   *                               user.
   * @param minQuality             {@link Double} minimum quality specified by the user.
   * @param progress               {@link ProgressBar} to indicate runtime information.
   */
  public SampleAnalyserRunnable(String sampleName, FeatureAnalysisEntry featureAnalysisEntry,
                                CloseableIterator<VariantContext> variantInformation,
                                VariablePositionsTable variablePositionsTable,
                                double minCoverage, double minHomFrequency, double minHetFrequency,
                                double maxHetFrequency, double minQuality,
                                ProgressBar progress) {
    this.sampleName = sampleName;
    this.featureAnalysisEntry = featureAnalysisEntry;
    this.variantInformation = variantInformation;
    this.variablePositionsTable = variablePositionsTable;
    this.minCoverage = minCoverage;
    this.minHomFrequency = minHomFrequency;
    this.minHetFrequency = minHetFrequency;
    this.maxHetFrequency = maxHetFrequency;
    this.minQuality = minQuality;
    this.progress = progress;
  }

  /**
   * Runs the threaded analysis of the data specified with the instances properties.
   */
  @Override
  public void run() {
    try {
      // For each variant information contained within the `.vcf` file of the respective sample.
      variablePositionsTable.putSampleInReference( featureAnalysisEntry.identifier, sampleName );
      while (variantInformation.hasNext()) {
        VariantContext variantContext = variantInformation.next();
        String variantContig = variantContext.getContig();
        // If the variant is on a contig/location other than the one specified in the reference analysis entry, the
        // variant can be skipped.
        if (!variantContig.equals(featureAnalysisEntry.referenceSequenceLocation)) {
          continue;
        }
        int variantPosition = variantContext.getStart();
        GenotypeType genotype = variantContext.getGenotypes().get(0).getType();
        if (genotype.equals(GenotypeType.MIXED) || genotype.equals(GenotypeType.UNAVAILABLE) ||
            genotype.equals(GenotypeType.NO_CALL)) {
          processNoCall(String.valueOf(variantPosition));
        } else {
          Allele referenceAllele = variantContext.getReference();
          List<Allele> alternateAlleles = variantContext.getAlternateAlleles();
          double variantCoverage = variantContext.getAttributeAsDouble("DP", 0.0);
          double variantQuality = variantContext.getPhredScaledQual();
          List<Double> allelesReadCounts = new ArrayList<>();
          if (alternateAlleles.size() == 0) {
            double referenceAlleleReadCount = variantContext.getAttributeAsDouble("DP", 0.0);
            // If the number of alternate alleles is zero, the variant is a no call, i.e. the reference allele is called.
            processReference(variantPosition, variantCoverage, variantQuality, referenceAlleleReadCount,
                referenceAllele);
          } else {
            String additionalAnnotation = "";
            // If a SnpEff annotation attribute (key is EFF) is present, add its content as additional annotation.
            if (variantContext.hasAttribute("EFF")) {
              additionalAnnotation = "EFF=" + variantContext.getAttributeAsString("EFF", "");
            }
            // Assign allele read counts.
            int[] alleleReadCounts = variantContext.getGenotype(0).getAD();
            for (int i = 1; i < alleleReadCounts.length; i++) { // Skip index 0 as it represents reference genotype.
              allelesReadCounts.add((double) alleleReadCounts[i]);
            }
            // For each alternate allele:
            for (int i = 0; i < alternateAlleles.size(); i++) {
              processVariation(variantPosition, variantQuality, variantCoverage, allelesReadCounts.get(i),
                  alternateAlleles.get(i), referenceAllele, additionalAnnotation, i, alternateAlleles.size() == 1);
            }
          }
        }
      }
      progress.step();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Internal method to process a position at which no call was possible.
   *
   * @param position {@link String} the position.
   */
  private void processNoCall(String position) {
    variablePositionsTable.putVariablePosition(
        featureAnalysisEntry.identifier,
        sampleName,
        position,
        new VariablePosition(
            VariablePosition.NO_CALL,
            NO_CALL
        ),
        0
    );
  }

  /**
   * Internal method to process a reference call.
   *
   * @param position        {@link String} the position.
   * @param coverage        {@link Double} the coverage of the position.
   * @param quality         {@link Double} the quality of the position.
   * @param readCount       {@link Double} the number of reads supporting the reference call.
   * @param referenceAllele {@link Allele} the reference allele information.
   */
  private void processReference(int position, double coverage, double quality, double readCount,
                                Allele referenceAllele) {
    // Access the reference alleles bases as character array. It should always have a length of one!
    char[] basesArray = referenceAllele.getBaseString().toCharArray();
    char base = VariablePosition.REFERENCE;
    String genomePosition;
    StringBuilder annotationStringBuilder = new StringBuilder();
    // Assume that the reference call fails on some filters, thus no further processing is necessary.
    boolean process = false;
    for (int i = 0; i < basesArray.length; i++) {
      if (coverage < minCoverage) {
        base = VariablePosition.REFERENCE_DISCARDED;
        annotationStringBuilder.append(LOW_COVERAGE + ANN_SEP);
        process = true;
      }
      if (quality < minQuality) {
        base = VariablePosition.REFERENCE_DISCARDED;
        annotationStringBuilder.append(LOW_QUALITY + ANN_SEP);
        process = true;
      }
      if ((readCount / coverage) < minHomFrequency) {
        base = VariablePosition.REFERENCE_DISCARDED;
        annotationStringBuilder.append(LOW_FREQUENCY + ANN_SEP);
        process = true;
      }
      if (process) {
        genomePosition = String.valueOf(position + i);
        annotationStringBuilder.append("COVERAGE=").append(coverage).append(ANN_SEP);
        annotationStringBuilder.append("QUALITY=").append(quality).append(ANN_SEP);
        annotationStringBuilder.append("FREQUENCY=").append((readCount / coverage)).append(ANN_SEP);
        String annotationString = annotationStringBuilder.toString();
        annotationString = annotationString.substring(0, annotationString.length() - 1);
        // Build and pass a variablePosition instance to the variablePositionsTable.
        variablePositionsTable.putVariablePosition(
            featureAnalysisEntry.identifier,
            sampleName,
            genomePosition,
            new VariablePosition(
                base,
                annotationString
            ),
            0
        );
      }
    }
  }

  /**
   * Internal method to process a alternate base call.
   *
   * @param variantPosition      {@link String} the position of the variant.
   * @param variantCoverage      {@link Double} the coverage of the variant.
   * @param variantQuality       {@link Double} the quality of the variant.
   * @param variantReadCount     {@link Double} the number of reads supporting the variant.
   * @param alternateAllele      {@link Allele} the alternative allele information.
   * @param referenceAllele      {@link Allele} the reference allele information.
   * @param additionalAnnotation {@link String} an additional annotation string that is added to the entry.
   * @param alleleIndex          {@link Integer} representing the index of the allele.
   * @param homCall              {@link Boolean} to indicate if a hom. or het. call is made.
   */
  private void processVariation(int variantPosition, double variantQuality, double variantCoverage,
                                double variantReadCount, Allele alternateAllele, Allele referenceAllele,
                                String additionalAnnotation, int alleleIndex, boolean homCall) {
    char[] alternateBasesArray = alternateAllele.getBaseString().toCharArray();
    char[] referenceBasesArray = referenceAllele.getBaseString().toCharArray();
    char base;
    String position;
    double variantFrequency = variantReadCount / variantCoverage;
    if (alternateBasesArray.length == referenceBasesArray.length) {
      // CASE: If the number of characters in the alternative and reference allele are equal, a simple base
      // substitution is observed.
      for (int i = 0; i < alternateBasesArray.length; i++) {
        position = String.valueOf(variantPosition + i);
        base = assignBase(alternateBasesArray[i], variantQuality, variantCoverage, variantFrequency);
        StringBuilder annotationStringBuilder = new StringBuilder();
        if (variantCoverage < minCoverage) {
          annotationStringBuilder.append(LOW_COVERAGE + ANN_SEP);
        }
        if (variantQuality < minQuality) {
          annotationStringBuilder.append(LOW_QUALITY + ANN_SEP);
        }
        if (homCall) {
          if (variantFrequency < minHomFrequency) {
            annotationStringBuilder.append(LOW_FREQUENCY + ANN_SEP);
          }
        } else {
          if (variantFrequency < minHetFrequency || maxHetFrequency < variantFrequency) {
            annotationStringBuilder.append(LOW_FREQUENCY + ANN_SEP);
          }
        }
        annotationStringBuilder.append("COVERAGE=").append(variantCoverage).append(ANN_SEP);
        annotationStringBuilder.append("QUALITY=").append(variantQuality).append(ANN_SEP);
        annotationStringBuilder.append("FREQUENCY=").append(variantFrequency).append(ANN_SEP);
        if (!additionalAnnotation.equals("")) {
          annotationStringBuilder.append(additionalAnnotation).append(ANN_SEP);
        }
        String annotationString = annotationStringBuilder.toString();
        annotationString = annotationString.substring(0, annotationString.length() - 1);
        variablePositionsTable.putVariablePosition(
            featureAnalysisEntry.identifier,
            sampleName,
            position,
            new VariablePosition(
                base,
                annotationString
            ),
            alleleIndex
        );
      }
    } else if (alternateBasesArray.length > referenceBasesArray.length) {
      // CASE: If the number of characters in the alternative allele is greater than in the reference allele, an
      // insertion is observed. NOTE: As .vcf files may contain ambiguous indel entries all indel cases are treated
      // as if as long as the alternate and reference bases match no indel has started.
      StringBuilder annotationStringBuilder = new StringBuilder();
      if (variantCoverage < minCoverage) {
        annotationStringBuilder.append(LOW_COVERAGE + ANN_SEP);
      }
      if (variantQuality < minQuality) {
        annotationStringBuilder.append(LOW_QUALITY + ANN_SEP);
      }
      if (homCall) {
        if (variantFrequency < minHomFrequency) {
          annotationStringBuilder.append(LOW_FREQUENCY + ANN_SEP);
        }
      } else {
        if (variantFrequency < minHetFrequency || maxHetFrequency < variantFrequency) {
          annotationStringBuilder.append(LOW_FREQUENCY + ANN_SEP);
        }
      }
      annotationStringBuilder.append("COVERAGE=").append(variantCoverage).append(ANN_SEP);
      annotationStringBuilder.append("QUALITY=").append(variantQuality).append(ANN_SEP);
      annotationStringBuilder.append("FREQUENCY=").append(variantFrequency).append(ANN_SEP);
      if (!additionalAnnotation.equals("")) {
        annotationStringBuilder.append(additionalAnnotation).append(ANN_SEP);
      }
      String annotationString = annotationStringBuilder.toString();
      annotationString = annotationString.substring(0, annotationString.length() - 1);
      for (int i = 0; i < alternateBasesArray.length; i++) {
        // If the reference analysis entry is on the negative sense strand the genome position and base have to be
        // converted.
        if (i < referenceBasesArray.length) {
          continue;
        }
        position = (variantPosition + i) + "I".repeat(i);
        base = assignBase(alternateBasesArray[i], variantQuality, variantCoverage, variantFrequency);
        variablePositionsTable.putVariablePosition(
            featureAnalysisEntry.identifier,
            sampleName,
            position,
            new VariablePosition(
                base,
                annotationString
            ),
            alleleIndex
        );
      }
    } else {
      // CASE: If the number of characters in the reference allele is greater than in the alternative allele, an
      // deletion is observed. NOTE: As .vcf files may contain ambiguous indel entries all indel cases are treated
      // as if as long as the alternate and reference bases match no indel has started.
      StringBuilder annotationStringBuilder = new StringBuilder();
      if (variantCoverage < minCoverage) {
        annotationStringBuilder.append(LOW_COVERAGE + ANN_SEP);
      }
      if (variantQuality < minQuality) {
        annotationStringBuilder.append(LOW_QUALITY + ANN_SEP);
      }
      if (homCall) {
        if (variantFrequency < minHomFrequency) {
          annotationStringBuilder.append(LOW_FREQUENCY + ANN_SEP);
        }
      } else {
        if (variantFrequency < minHetFrequency || maxHetFrequency < variantFrequency) {
          annotationStringBuilder.append(LOW_FREQUENCY + ANN_SEP);
        }
      }
      annotationStringBuilder.append("COVERAGE=").append(variantCoverage).append(ANN_SEP);
      annotationStringBuilder.append("QUALITY=").append(variantQuality).append(ANN_SEP);
      annotationStringBuilder.append("FREQUENCY=").append(variantFrequency).append(ANN_SEP);
      if (!additionalAnnotation.equals("")) {
        annotationStringBuilder.append(additionalAnnotation).append(ANN_SEP);
      }
      String annotationString = annotationStringBuilder.toString();
      annotationString = annotationString.substring(0, annotationString.length() - 1);
      for (int i = 0; i < referenceBasesArray.length; i++) {
        if (i < alternateBasesArray.length) {
          continue;
        }
        base = VariablePosition.DELETION;
        position = String.valueOf(variantPosition + i);
        variablePositionsTable.putVariablePosition(
            featureAnalysisEntry.identifier,
            sampleName,
            position,
            new VariablePosition(
                base,
                annotationString
            ),
            alleleIndex
        );
      }
    }
  }

  /**
   * Internal method to assign the correct content symbol based on the called base, quality, coverage and frequency.
   * <p>
   * See the implementation of {@link VariablePosition} for all possible contents.
   *
   * @param alternativeBase {@link Character} the called alternative base.
   * @param quality         {@link Double} the quality of the base call.
   * @param coverage        {@link Double} the coverage of the base call.
   * @param frequency       {@link Double} the frequency of the base call.
   * @return A {@link Character} representing the variant content.
   */
  private char assignBase(char alternativeBase, double quality, double coverage,
                          double frequency) {
    if ((coverage < minCoverage) || (quality < minQuality) || (frequency < minHomFrequency)) {
      switch (alternativeBase) {
        case VariablePosition.ALT_A:
          return VariablePosition.ALT_A_DISCARDED;
        case VariablePosition.ALT_C:
          return VariablePosition.ALT_C_DISCARDED;
        case VariablePosition.ALT_G:
          return VariablePosition.ALT_G_DISCARDED;
        case VariablePosition.ALT_T:
          return VariablePosition.ALT_T_DISCARDED;
      }
    } else {
      return switch (alternativeBase) {
        case VariablePosition.ALT_A -> VariablePosition.ALT_A;
        case VariablePosition.ALT_C -> VariablePosition.ALT_C;
        case VariablePosition.ALT_G -> VariablePosition.ALT_G;
        case VariablePosition.ALT_T -> VariablePosition.ALT_T;
        default -> '?'; // TODO: Temporary solution to indicate ambiguous bases.
      };
    }
    return '?';
  }

}
