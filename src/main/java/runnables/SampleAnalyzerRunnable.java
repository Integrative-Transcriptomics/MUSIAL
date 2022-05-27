package runnables;

import static org.forester.util.ForesterUtil.round;

import datastructure.FeatureEntry;
import datastructure.SampleEntry;
import datastructure.VariantsDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import me.tongfei.progressbar.ProgressBar;
import utility.Bio;
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
public final class SampleAnalyzerRunnable implements Runnable {

  private final SampleEntry sampleEntry;
  private final FeatureEntry featureEntry;
  private final VariantsDictionary variantsDictionary;
  private final ProgressBar pb;
  private int ambiguousVariantsCount = 0;

  public SampleAnalyzerRunnable(SampleEntry sampleEntry, FeatureEntry featureEntry,
                                VariantsDictionary variantsDictionary, ProgressBar pb) {
    this.sampleEntry = sampleEntry;
    this.featureEntry = featureEntry;
    this.variantsDictionary = variantsDictionary;
    this.pb = pb;
  }

  /**
   * Runs the threaded analysis of the data specified with the instances properties.
   */
  @Override
  public void run() {
    try {
      Iterator<VariantContext> variantContextIterator = sampleEntry.vcfFileReader.query(
          featureEntry.chromosome, featureEntry.start, featureEntry.end
      );
      // For each variant information contained within the `.vcf` file of the respective sample.
      while (variantContextIterator.hasNext()) {
        VariantContext variantContext = variantContextIterator.next();
        String variantContig = variantContext.getContig();
        if (!variantContig.equals(variantsDictionary.chromosome)) {
          // If the variant is on a chromosome other than the one specified in the variants dictionary it is skipped.
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
            // The allelic frequencies are extracted to order the alleles with respect to decreasing frequency.
            try {
              List<Double> allelesFrequencies = variantContext.getAttributeAsDoubleList("AF", 0.0);
              TreeMap<Double, Allele> sortedAlternateAlleles = new TreeMap<>(Collections.reverseOrder());
              for (int i = 0; i < alternateAlleles.size(); i++) {
                sortedAlternateAlleles.put(round(allelesFrequencies.get(i), 2), alternateAlleles.get(i));
              }
              boolean isPrimary = true;
              for (Map.Entry<Double, Allele> entry : sortedAlternateAlleles.entrySet()) {
                variantFrequency = entry.getKey();
                variantAllele = entry.getValue();
                // TODO: Check for SnpEff annotations; stored as attribute with key 'ANN'.
                processVariantCall(variantPosition, variantQuality, variantCoverage, variantFrequency, variantAllele,
                    referenceAllele, alternateAlleles.size() > 1, isPrimary);
                isPrimary = false;
              }
            } catch (Exception e) {
              throw e;
              /*
              Logging.logWarning(
                  "An error occurred during the analysis of sample " + sampleEntry.name + " at position " +
                      variantPosition + IO.LINE_SEPARATOR + "\t\t> " +
                      e.getMessage());
               */
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
      throw e;
    } finally {
      pb.step();
      if (ambiguousVariantsCount > 0) {
        Logging.logWarning("Sample " + sampleEntry.name + " comprises " + ambiguousVariantsCount + " ambiguous " +
            "(accepted primary) variants for feature " + featureEntry.name + ".");
      }
    }
  }

  /**
   * Internal method to process a variant-call.
   *
   * @param variantPosition  {@link String} the position of the variant.
   * @param variantQuality   {@link Double} the quality of the variant.
   * @param variantCoverage  {@link Double} the coverage of the variant.
   * @param variantFrequency {@link Double} the frequency of the variant.
   * @param variantAllele    {@link Allele} the alternative allele information.
   * @param referenceAllele  {@link Allele} the reference allele information.
   * @param isHetCall        {@link Boolean} to indicate if the variant originates from a het. call; necessary to
   *                         apply the correct filter method.
   * @param isPrimary        {@link Boolean} to indicate if the variant originates from the most frequent allele.
   */
  private void processVariantCall(int variantPosition, double variantQuality, double variantCoverage,
                                  double variantFrequency, Allele variantAllele, Allele referenceAllele,
                                  boolean isHetCall, boolean isPrimary) {
    String variantAlleleContent = variantAllele.getBaseString();
    String referenceAlleleContent = referenceAllele.getBaseString();
    boolean isRejected = variantQuality < variantsDictionary.parameters.minQuality ||
        variantCoverage < variantsDictionary.parameters.minCoverage ||
        (!isHetCall && variantFrequency < variantsDictionary.parameters.minFrequency) ||
        (isHetCall && (variantFrequency < variantsDictionary.parameters.minHetFrequency ||
            variantFrequency > variantsDictionary.parameters.maxHetFrequency));
    if (variantAlleleContent.length() == referenceAlleleContent.length() &&
        referenceAlleleContent.length() == 1) {
      // CASE: Nucleotide substitution.
      variantsDictionary
          .addVariant(featureEntry.name, variantPosition, variantAlleleContent, referenceAlleleContent,
              sampleEntry.name, isPrimary, isRejected, variantQuality, variantCoverage, variantFrequency);
    } else if (variantAlleleContent.length() > referenceAlleleContent.length() &&
        referenceAlleleContent.length() == 1) {
      // CASE: Insertion.
      variantsDictionary
          .addVariant(featureEntry.name, variantPosition, variantAlleleContent, referenceAlleleContent,
              sampleEntry.name, isPrimary, isRejected, variantQuality, variantCoverage, variantFrequency);
    } else if (variantAlleleContent.length() < referenceAlleleContent.length() &&
        variantAlleleContent.length() == 1) {
      // CASE: Deletion
      variantAlleleContent = variantAlleleContent + "-".repeat(referenceAlleleContent.length() - 1);
      variantsDictionary
          .addVariant(featureEntry.name, variantPosition, variantAlleleContent, referenceAlleleContent,
              sampleEntry.name, isPrimary, isRejected, variantQuality, variantCoverage, variantFrequency);
    } else {
      ArrayList<String> resolvedVariants =
          Bio.resolveAmbiguousVariant(
              referenceAllele.getBaseString(),
              variantAllele.getBaseString()
          );
      for (String resolvedVariant : resolvedVariants) {
        int resolvedPosition = variantPosition + Integer.parseInt(resolvedVariant.split("@")[1]) - 1;
        String resolvedVariantAlleleContent = resolvedVariant.split("@")[2];
        String resolvedReferenceAlleleContent = resolvedVariant.split("@")[3];
        variantsDictionary
            .addVariant(featureEntry.name, resolvedPosition, resolvedVariantAlleleContent,
                resolvedReferenceAlleleContent, sampleEntry.name, isPrimary, isRejected, variantQuality,
                variantCoverage, variantFrequency);
        if (isPrimary && !isRejected) {
          this.ambiguousVariantsCount += 1;
        }
      }
    }
  }

}
