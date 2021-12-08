package runnables;

import static org.forester.util.ForesterUtil.round;

import components.SnpEffAnnotator;
import datastructure.FeatureAnalysisEntry;
import datastructure.VariantContent;
import datastructure.VariantPositionsTable;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import main.Musial;
import me.tongfei.progressbar.ProgressBar;
import utility.IO;
import utility.Logging;

/**
 * Implementation of the {@link Runnable} interface to analyze a sample.
 * <p>
 * Runs a sample analysis in a single thread.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
@SuppressWarnings("ClassCanBeRecord")
public final class SampleAnalyserRunnable implements Runnable {

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
   * An {@link VariantPositionsTable} instance to which the parsed information about variants is passed.
   */
  private final VariantPositionsTable variantPositionsTable;
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
   * @param sampleName            {@link String} the samples name.
   * @param featureAnalysisEntry  {@link FeatureAnalysisEntry} information about the reference to which respect
   *                              the sample is analyzed.
   * @param variantInformation    {@link CloseableIterator<VariantContext>} returned by
   *                              a {@link htsjdk.variant.vcf.VCFFileReader}
   *                              giving access to the single input `.vcf`
   *                              files entries.
   * @param variantPositionsTable {@link VariantPositionsTable} to store the processed information.
   * @param minCoverage           {@link Double} minimum coverage specified by the user.
   * @param minHomFrequency       {@link Double} minimum hom. call frequency (wrt. read support) specified by the
   *                              user.
   * @param minHetFrequency       {@link Double} minimum het. call frequency (wrt. read support) specified by the
   *                              user.
   * @param maxHetFrequency       {@link Double} maximum het. call frequency (wrt. read support) specified by the
   *                              user.
   * @param minQuality            {@link Double} minimum quality specified by the user.
   * @param progress              {@link ProgressBar} to indicate runtime information.
   */
  public SampleAnalyserRunnable(String sampleName, FeatureAnalysisEntry featureAnalysisEntry,
                                CloseableIterator<VariantContext> variantInformation,
                                VariantPositionsTable variantPositionsTable,
                                double minCoverage, double minHomFrequency, double minHetFrequency,
                                double maxHetFrequency, double minQuality,
                                ProgressBar progress) {
    this.sampleName = sampleName;
    this.featureAnalysisEntry = featureAnalysisEntry;
    this.variantInformation = variantInformation;
    this.variantPositionsTable = variantPositionsTable;
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
      while (variantInformation.hasNext()) {
        VariantContext variantContext = variantInformation.next();
        String variantContig = variantContext.getContig();
        if (!variantContig.equals(featureAnalysisEntry.referenceSequenceLocation)) {
          // If the variant is on a contig/location other than the one specified in the reference analysis entry, the
          // variant can be skipped.
          continue;
        }
        int variantPosition = variantContext.getStart();
        double variantQuality = round(variantContext.getPhredScaledQual(), 2);
        double variantCoverage = variantContext.getAttributeAsDouble("DP", 0.0);
        double variantFrequency;
        // Check if the position is a no call position; true if a quality value of -10.0 is returned.
        boolean isNoCall = variantQuality == -10.0;
        if (isNoCall) {
          // CASE: The variant context is a no call.
          // FIXME: Possibly extend implementation if no-call information is provided.
          //noinspection UnnecessaryContinue
          continue;
        } else {
          Allele referenceAllele = variantContext.getReference();
          List<Allele> alternateAlleles = variantContext.getAlternateAlleles();
          Allele variantAllele;
          if (alternateAlleles.size() != 0) {
            // CASE: The variant context reflects a variant call.
            // Check for additional annotations.
            HashMap<String, String> additionalAnnotation = new HashMap<>();
            // FIXME: Currently only the possibility of an additional SnpEff annotation exists.
            // If a SnpEff annotation attribute (key is ANN, LOF and NMD are currently ignored) is present, add its
            // content as additional annotation if it affects the gene feature to analyze.
            String annotationString;
            StringBuilder filteredAnnotationStringBuilder = new StringBuilder();
            List<String> annotations;
            if (variantContext.hasAttribute(SnpEffAnnotator.snpEffAnnotationAttributeKey)) {
              annotationString = variantContext.getAttributeAsString(SnpEffAnnotator.snpEffAnnotationAttributeKey, "");
              annotationString = annotationString.replace("[", "");
              annotationString = annotationString.replace("]", "");
              if (featureAnalysisEntry.isGene) {
                annotations = Arrays.asList(annotationString.split(","));
                for (String annotation : annotations) {
                  if (annotation.split("\\|")[3].equals(featureAnalysisEntry.identifier)) {
                    filteredAnnotationStringBuilder.append(annotation).append(",");
                  }
                }
                if (filteredAnnotationStringBuilder.length() != 0) {
                  filteredAnnotationStringBuilder.setLength(filteredAnnotationStringBuilder.length() - 1);
                  additionalAnnotation
                      .put(SnpEffAnnotator.snpEffAnnotationAttributeKey, filteredAnnotationStringBuilder.toString());
                }
              } else {
                additionalAnnotation.put(SnpEffAnnotator.snpEffAnnotationAttributeKey, annotationString);
              }
            }
            // Initiate sample call type counting.
            variantPositionsTable
                .countCallType(featureAnalysisEntry.name, sampleName, alternateAlleles.size() > 1);
            // The allelic depths of coverage are extracted to order the alleles with respect to their decreasing
            // frequency.
            // FIXME: Each genotype should contain the full list of allele depths.
            try {
              int[] allelesReadCounts = variantContext.getGenotype(0).getAD();
              TreeMap<Double, Allele> sortedAlternateAlleles = new TreeMap<>(Collections.reverseOrder());
              for (int i = 0; i < alternateAlleles.size(); i++) {
                sortedAlternateAlleles.put(round(allelesReadCounts[i + 1] / variantCoverage, 2),
                    alternateAlleles.get(i));
              }
              boolean isMostFrequent = true;
              for (Map.Entry<Double, Allele> entry : sortedAlternateAlleles.entrySet()) {
                variantFrequency = entry.getKey();
                variantAllele = entry.getValue();
                processVariantCall(variantPosition, variantQuality, variantCoverage, variantFrequency, variantAllele,
                    referenceAllele, alternateAlleles.size() > 1, isMostFrequent, additionalAnnotation);
                isMostFrequent = false;
              }
            } catch (TribbleException e) {
              Logging.logWarning(
                  "An error occurred during the analysis of sample " + sampleName + " at position " + variantPosition +
                      " within feature " + featureAnalysisEntry.name + IO.LINE_SEPARATOR + "> " + e.getMessage());
            }
          } else {
            // CASE: The variant context is a reference call.
            // FIXME: Possibly extend implementation if reference-call information is provided.
            //noinspection UnnecessaryContinue
            continue;
          }
        }
      }
    } catch (Exception e) {
      if (Musial.debug) {
        e.printStackTrace();
      } else {
        Logging.logWarning( "An error occurred during sample analysis: " + e.getMessage() );
      }
    } finally {
      progress.step();
    }
  }

  /**
   * Internal method to process a variant-call.
   *
   * @param variantPosition      {@link String} the position of the variant.
   * @param variantQuality       {@link Double} the quality of the variant.
   * @param variantCoverage      {@link Double} the coverage of the variant.
   * @param variantFrequency     {@link Double} the frequency of the variant.
   * @param variantAllele        {@link Allele} the alternative allele information.
   * @param referenceAllele      {@link Allele} the reference allele information.
   * @param isHetCall            {@link Boolean} to indicate if the variant originates from a het. call; necessary to
   *                             apply the correct filter method.
   * @param isMFA                {@link Boolean} to indicate if the variant originates from the most frequent allele.
   * @param additionalAnnotation {@link String} an additional annotation string that is added to the entry.
   */
  private void processVariantCall(int variantPosition, double variantQuality, double variantCoverage,
                                  double variantFrequency, Allele variantAllele, Allele referenceAllele,
                                  boolean isHetCall, boolean isMFA, HashMap<String, String> additionalAnnotation) {

    char[] variantContentArray = variantAllele.getBaseString().toCharArray();
    char[] referenceContentArray = referenceAllele.getBaseString().toCharArray();
    if (variantQuality < minQuality || variantCoverage < minCoverage ||
        (!isHetCall && variantFrequency < minHomFrequency) ||
        (isHetCall && (variantFrequency < minHetFrequency || variantFrequency > maxHetFrequency))) {
      variantPositionsTable.countRejectedCall(featureAnalysisEntry.name, sampleName,
          String.valueOf(variantPosition));
      return;
    }
    String position;
    char content;
    if (variantContentArray.length == referenceContentArray.length && referenceContentArray.length == 1) {
      // CASE: Nucleotide substitution.
      position = String.valueOf(variantPosition);
      content = variantAllele.getBaseString().toCharArray()[0];
      variantPositionsTable.putVariablePosition(featureAnalysisEntry.name, sampleName,
          position, content,
          variantQuality, variantCoverage, variantFrequency, isMFA, additionalAnnotation);
      variantPositionsTable.countAlternateCall(featureAnalysisEntry.name, sampleName, position, content);
    } else if (variantContentArray.length > referenceContentArray.length && referenceContentArray.length == 1) {
      // CASE: Insertion.
      for (int i = 1; i < variantContentArray.length; i++) {
        position = variantPosition + "+" + i;
        content = variantAllele.getBaseString().toCharArray()[i];
        variantPositionsTable.putVariablePosition(featureAnalysisEntry.name, sampleName,
            position, content, variantQuality, variantCoverage, variantFrequency, isMFA, additionalAnnotation);
        variantPositionsTable.countAlternateCall(featureAnalysisEntry.name, sampleName, position, content);
      }
    } else if (variantContentArray.length < referenceContentArray.length && variantContentArray.length == 1) {
      // CASE: Deletion
      for (int i = 1; i < referenceContentArray.length; i++) {
        position = String.valueOf(variantPosition + i);
        content = VariantContent.DELETION;
        variantPositionsTable.putVariablePosition(featureAnalysisEntry.name, sampleName,
            position, content, variantQuality, variantCoverage, variantFrequency, isMFA, additionalAnnotation);
        variantPositionsTable.countAlternateCall(featureAnalysisEntry.name, sampleName, position, content);
      }
    } else {
      // FIXME: Instead of using system printing a logging library should be used.
      Logging.logWarning("Skipping ambiguous variant context of sample " + sampleName + " at position " +
          variantPosition + " within the reference feature " + featureAnalysisEntry.name + IO.LINE_SEPARATOR +
          "> Reference content:\t" + referenceAllele.getBaseString() + IO.LINE_SEPARATOR + "> Sample content:\t" +
          variantAllele.getBaseString());
    }
  }

}
