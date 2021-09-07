package runnables;

import datastructure.ReferenceAnalysisEntry;
import datastructure.VariablePosition;
import datastructure.VariablePositionsTable;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
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
 * @TODO: 29.08.2021: The processing of SnpEff annotations is currently not implemented.
 * @TODO: 29.08.2021: The representation of indels in .vcf files is highly ambiguous. Fix implementation.
 */
public final class SampleAnalyserRunnable implements Runnable {

  /**
   * The internal name of the sample to analyze.
   */
  private final String sampleName;
  /**
   * An {@link ReferenceAnalysisEntry} instance storing information about the reference with regard to which the
   * sample is analyzed.
   */
  private final ReferenceAnalysisEntry referenceAnalysisEntry;
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
   * The minimum frequency to accept a variant specified by the user.
   */
  private final double minFrequency;
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
   * @param referenceAnalysisEntry {@link ReferenceAnalysisEntry} information about the reference to which respect
   *                               the sample is analyzed.
   * @param variantInformation     {@link CloseableIterator<VariantContext>} returned by
   *                               a {@link htsjdk.variant.vcf.VCFFileReader}
   *                               giving access to the single input `.vcf`
   *                               files entries.
   * @param variablePositionsTable {@link VariablePositionsTable} to store the processed information.
   * @param minCoverage            {@link Double} minimum coverage specified by the user.
   * @param minFrequency           {@link Double} minimum frequency specified by the user. This may only be relevant if
   *                               multiple alleles are possible.
   * @param minQuality             {@link Double} minimum quality specified by the user.
   * @param progress               {@link ProgressBar} to indicate runtime information.
   */
  public SampleAnalyserRunnable(String sampleName, ReferenceAnalysisEntry referenceAnalysisEntry,
                                CloseableIterator<VariantContext> variantInformation,
                                VariablePositionsTable variablePositionsTable,
                                double minCoverage, double minFrequency, double minQuality, ProgressBar progress) {
    this.sampleName = sampleName;
    this.referenceAnalysisEntry = referenceAnalysisEntry;
    this.variantInformation = variantInformation;
    this.variablePositionsTable = variablePositionsTable;
    this.minCoverage = minCoverage;
    this.minFrequency = minFrequency;
    this.minQuality = minQuality;
    this.progress = progress;
  }

  /**
   * Runs the threaded analysis of the data specified with the instances properties.
   */
  @Override
  public void run() {
    // For each variant information contained within the `.vcf` file of the respective sample.
    while (variantInformation.hasNext()) {
      VariantContext variantContext = variantInformation.next();
      String variantContig = variantContext.getContig();
      // If the variant is on a contig/location other than the one specified in the reference analysis entry, the
      // variant can be skipped.
      if (!variantContig.equals(referenceAnalysisEntry.analysisSequenceLocation)) {
        continue;
      }
      int variantPosition = variantContext.getStart();
      Allele referenceAllele = variantContext.getReference();
      List<Allele> alternateAlleles = variantContext.getAlternateAlleles();
      double variantCoverage = variantContext.getAttributeAsDouble("DP", 0.0);
      double variantQuality = variantContext.getPhredScaledQual();
      if (alternateAlleles.size() == 0) {
        // If the number of alternate alleles is zero, the variant is a no call, i.e. the reference allele is called.
        processReference(variantPosition, variantCoverage, variantQuality, referenceAllele);
      } else {
        // Else, for each alternate allele the variant is further processed.
        List<Double> alternateAllelesFrequency = variantContext.getAttributeAsDoubleList("AF", 0.0);
        for (int i = 0; i < alternateAlleles.size(); i++) {
          processVariation(variantPosition, variantQuality, variantCoverage, alternateAllelesFrequency.get(i),
              alternateAlleles.get(i), referenceAllele);
        }
      }
    }
    progress.step();
  }

  /**
   * Internal method to process a reference call.
   *
   * @param position        {@link String} the position of the variant.
   * @param coverage        {@link Double} the coverage of the variant.
   * @param quality         {@link Double} the quality of the variant.
   * @param referenceAllele {@link Allele} the reference allele information.
   */
  private void processReference(int position, double coverage, double quality,
                                Allele referenceAllele) {
    // Access the reference alleles bases as character array. It may contain more than one entry, if some alleles are
    // reference calls and other alleles represent deletions wrt. the reference.
    char[] basesArray = referenceAllele.getBaseString().toCharArray();
    char base;
    String genomePosition;
    // Assume that the reference call passes all filters, thus no further processing is necessary.
    boolean process;
    for (int i = 0; i < basesArray.length; i++) {
      if (coverage < minCoverage) {
        // If the reference call has low coverage, this information is stored in the variablePositionsTable.
        base = VariablePosition.REFERENCE_LOW_COV;
        process = true;
      } else if (quality < minQuality) {
        // If the reference call has low quality, this information is stored in the variablePositionsTable.
        base = VariablePosition.REFERENCE_LOW_QUAL;
        process = true;
      } else {
        // If the reference call passes all filters, no further processing is applied.
        base = VariablePosition.REFERENCE;
        process = false;
      }
      if (process) {
        // If the reference analysis entry is on the negative sense strand the genome position has to be converted.
        if ( this.referenceAnalysisEntry.isSense ) {
          genomePosition = String.valueOf( position + i );
        } else {
          genomePosition = String.valueOf(
              Bio.getPositionOnReverseComplement(
                  position + i,
                  referenceAnalysisEntry.analysisSequenceStart,
                  referenceAnalysisEntry.analysisSequenceEnd
              )
          );
        }
        // Build and pass a variablePosition instance to the variablePositionsTable.
        variablePositionsTable.putVariablePosition(
            referenceAnalysisEntry.analysisIdentifier,
            sampleName,
            genomePosition,
            new VariablePosition(
                base,
                coverage,
                quality,
                "",
                1.0
            )
        );
      }
    }
  }

  /**
   * Internal method to process a alternate base call.
   *
   * @param variantPosition  {@link String} the position of the variant.
   * @param variantCoverage  {@link Double} the coverage of the variant.
   * @param variantQuality   {@link Double} the quality of the variant.
   * @param variantFrequency {@link Double} the frequency of the variant. Should only be different from one if
   *                         multiple alleles are present.
   * @param alternateAllele  {@link Allele} the alternative allele information.
   * @param referenceAllele  {@link Allele} the reference allele information.
   */
  private void processVariation(int variantPosition, double variantQuality, double variantCoverage,
                                double variantFrequency, Allele alternateAllele, Allele referenceAllele) {
    char[] alternateBasesArray = alternateAllele.getBaseString().toCharArray();
    char[] referenceBasesArray = referenceAllele.getBaseString().toCharArray();
    char base;
    String position;
    if (alternateBasesArray.length == referenceBasesArray.length) {
      // CASE: If the number of characters in the alternative and reference allele are equal, a simple base
      // substitution is observed.
      for (int i = 0; i < alternateBasesArray.length; i++) {
        // If the reference analysis entry is on the negative sense strand the genome position and base have to be
        // converted.
        if ( this.referenceAnalysisEntry.isSense ) {
          position = String.valueOf( variantPosition + i );
          base = assignBase(alternateBasesArray[i], variantQuality, variantCoverage, variantFrequency);
        } else {
          position = String.valueOf(
              Bio.getPositionOnReverseComplement(
                  variantPosition + i,
                  referenceAnalysisEntry.analysisSequenceStart,
                  referenceAnalysisEntry.analysisSequenceEnd
              )
          );
          base = assignBase(Bio.invertBase(alternateBasesArray[i]), variantQuality, variantCoverage, variantFrequency);
        }
        variablePositionsTable.putVariablePosition(
            referenceAnalysisEntry.analysisIdentifier,
            sampleName,
            position,
            new VariablePosition(
                base,
                variantCoverage,
                variantQuality,
                "",
                variantFrequency
            )
        );
      }
    } else if (alternateBasesArray.length > referenceBasesArray.length) {
      // CASE: If the number of characters in the alternative allele is greater than in the reference allele, an
      // insertion is observed. Currently only the simple case is implemented, i.e. the first base in both alleles is
      // equal and all following bases in the alternative allele are "inserted".
      for (int i = 1; i < alternateBasesArray.length; i++) {
        // If the reference analysis entry is on the negative sense strand the genome position and base have to be
        // converted.
        if ( this.referenceAnalysisEntry.isSense ) {
          position = (variantPosition + i) + "I".repeat(i);
          base = assignBase(alternateBasesArray[i], variantQuality, variantCoverage, variantFrequency);
        } else {
          position = Bio.getPositionOnReverseComplement(
              variantPosition + i,
              referenceAnalysisEntry.analysisSequenceStart,
              referenceAnalysisEntry.analysisSequenceEnd
          ) + "I".repeat(alternateBasesArray.length - i);
          base = assignBase(Bio.invertBase(alternateBasesArray[i]), variantQuality, variantCoverage, variantFrequency);
        }
        variablePositionsTable.putVariablePosition(
            referenceAnalysisEntry.analysisIdentifier,
            sampleName,
            position,
            new VariablePosition(
                base,
                variantCoverage,
                variantQuality,
                "",
                variantFrequency
            )
        );
      }
    } else {
      // CASE: If the number of characters in the reference allele is greater than in the alternative allele, an
      // deletion is observed. Currently only the simple case is implemented, i.e. the first base in both alleles is
      // equal and all following bases in the reference allele were deleted in the alternative allele.
      for (int i = 1; i < referenceBasesArray.length; i++) {
        base = VariablePosition.DELETION;
        // If the reference analysis entry is on the negative sense strand the genome position and base have to be
        // converted.
        if ( this.referenceAnalysisEntry.isSense ) {
          position = String.valueOf( variantPosition + i );
        } else {
          position = String.valueOf(
              Bio.getPositionOnReverseComplement(
                  variantPosition + i,
                  referenceAnalysisEntry.analysisSequenceStart,
                  referenceAnalysisEntry.analysisSequenceEnd
              )
          );
        }
        variablePositionsTable.putVariablePosition(
            referenceAnalysisEntry.analysisIdentifier,
            sampleName,
            position,
            new VariablePosition(
                base,
                variantCoverage,
                variantQuality,
                "",
                variantFrequency
            )
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
   * @param quality {@link Double} the quality of the base call.
   * @param coverage {@link Double} the coverage of the base call.
   * @param frequency {@link Double} the frequency of the base call.
   * @return A {@link Character} representing the variant content.
   */
  private char assignBase(char alternativeBase, double quality, double coverage,
                          double frequency) {
    if (coverage < minCoverage) {
      switch (alternativeBase) {
        case VariablePosition.ALT_A:
          return VariablePosition.ALT_A_LOW_COV;
        case VariablePosition.ALT_C:
          return VariablePosition.ALT_C_LOW_COV;
        case VariablePosition.ALT_G:
          return VariablePosition.ALT_G_LOW_COV;
        case VariablePosition.ALT_T:
          return VariablePosition.ALT_T_LOW_COV;
      }
    } else if (quality < minQuality) {
      switch (alternativeBase) {
        case VariablePosition.ALT_A:
          return VariablePosition.ALT_A_LOW_QUAL;
        case VariablePosition.ALT_C:
          return VariablePosition.ALT_C_LOW_QUAL;
        case VariablePosition.ALT_G:
          return VariablePosition.ALT_G_LOW_QUAL;
        case VariablePosition.ALT_T:
          return VariablePosition.ALT_T_LOW_QUAL;
      }
    } else if (frequency < minFrequency) {
      switch (alternativeBase) {
        case VariablePosition.ALT_A:
          return VariablePosition.ALT_A_LOW_FREQ;
        case VariablePosition.ALT_C:
          return VariablePosition.ALT_C_LOW_FREQ;
        case VariablePosition.ALT_G:
          return VariablePosition.ALT_G_LOW_FREQ;
        case VariablePosition.ALT_T:
          return VariablePosition.ALT_T_LOW_FREQ;
      }
    } else {
      switch (alternativeBase) {
        case VariablePosition.ALT_A:
          return VariablePosition.ALT_A;
        case VariablePosition.ALT_C:
          return VariablePosition.ALT_C;
        case VariablePosition.ALT_G:
          return VariablePosition.ALT_G;
        case VariablePosition.ALT_T:
          return VariablePosition.ALT_T;
      }
    }
    return '?'; // TODO: Temporary solution to indicate ambiguous bases.
  }

}
