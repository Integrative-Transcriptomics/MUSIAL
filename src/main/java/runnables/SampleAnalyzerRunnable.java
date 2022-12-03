package runnables;

import components.Bio;
import components.Logging;
import datastructure.FeatureEntry;
import datastructure.SampleEntry;
import datastructure.VariantsDictionary;
import exceptions.MusialException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.javatuples.Triplet;

import java.util.*;

import static org.forester.util.ForesterUtil.round;

/**
 * Implementation of the {@link Runnable} interface to analyze a sample.
 * <p>
 * Runs a sample analysis in a single thread.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class SampleAnalyzerRunnable implements Runnable {

    /**
     * The {@link SampleEntry} to be analyzed.
     */
    private final SampleEntry sampleEntry;
    /**
     * The {@link FeatureEntry} to be analyzed.
     */
    private final FeatureEntry featureEntry;
    /**
     * The {@link VariantsDictionary} instance the sample and feature originate from.
     */
    private final VariantsDictionary variantsDictionary;

    /**
     * Constructor of the {@link SampleAnalyzerRunnable} class.
     *
     * @param sampleEntry        The {@link SampleEntry} to be analyzed.
     * @param featureEntry       The {@link FeatureEntry} to be analyzed.
     * @param variantsDictionary The {@link VariantsDictionary} instance the sample and feature originate from.
     */
    public SampleAnalyzerRunnable(SampleEntry sampleEntry, FeatureEntry featureEntry,
                                  VariantsDictionary variantsDictionary) {
        this.sampleEntry = sampleEntry;
        this.featureEntry = featureEntry;
        this.variantsDictionary = variantsDictionary;
    }

    /**
     * Runs the threaded analysis of the data specified with the instances properties.
     */
    @Override
    public void run() {
        try {
            // Query variants located in the analyzed features' region.
            Iterator<VariantContext> variantContextIterator = sampleEntry.vcfFileReader.query(
                    featureEntry.chromosome, featureEntry.start, featureEntry.end
            );
            while (variantContextIterator.hasNext()) {
                VariantContext variantContext = variantContextIterator.next();
                String variantContig = variantContext.getContig();
                if (!variantContig.equals(featureEntry.chromosome)) {
                    // If the variant is on a chromosome other than the one specified in the variants' dictionary it is skipped.
                    continue;
                }
                int variantPosition = variantContext.getStart();
                if ( variantsDictionary.excludedPositions.get(variantContig).contains(variantPosition) ) {
                    continue;
                }
                double variantQuality = round(variantContext.getPhredScaledQual(), 2);
                double variantCoverage = variantContext.getAttributeAsDouble("DP", 0.0);
                double variantFrequency;
                // Check if the position is a no call position; true if a quality value of -10.0 is returned.
                boolean isNoCall = variantQuality == -10.0;
                if (isNoCall) {
                    // CASE: The variant context is a no call.
                    // TODO: Possibly extend implementation if no-call information is provided.
                    //noinspection UnnecessaryContinue
                    continue;
                } else {
                    Allele referenceAllele = variantContext.getReference();
                    List<Allele> alternateAlleles = variantContext.getAlternateAlleles();
                    Allele variantAllele;
                    if (alternateAlleles.size() != 0) {
                        // The allelic frequencies are extracted to order the alleles with respect to decreasing frequency.
                        List<Double> allelesFrequencies = variantContext.getAttributeAsDoubleList("AF", 0.0);
                        TreeMap<Double, Allele> sortedAlternateAlleles = new TreeMap<>(Collections.reverseOrder());
                        for (int i = 0; i < alternateAlleles.size(); i++) {
                            sortedAlternateAlleles.put(round(allelesFrequencies.get(i), 2), alternateAlleles.get(i));
                        }
                        boolean isPrimary = true;
                        for (Map.Entry<Double, Allele> entry : sortedAlternateAlleles.entrySet()) {
                            variantFrequency = entry.getKey();
                            variantAllele = entry.getValue();
                            processVariantCall(variantPosition, variantQuality, variantCoverage, variantFrequency, variantAllele,
                                    referenceAllele, alternateAlleles.size() > 1, isPrimary);
                            isPrimary = false;
                        }
                    } else {
                        // CASE: The variant context is a reference call.
                        // TODO: Possibly extend implementation if reference-call information is provided.
                        //noinspection UnnecessaryContinue
                        continue;
                    }
                }
            }
        } catch (Exception e) {
            Logging.logError("An error occurred during the analysis of sample " + sampleEntry.name + " (file: " + sampleEntry.vcfFile.getAbsolutePath() + "); " + e.getMessage());
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
                                    boolean isHetCall, boolean isPrimary) throws MusialException {
        String variantAlleleContent = variantAllele.getBaseString();
        String referenceAlleleContent = referenceAllele.getBaseString();
        //noinspection unused
        boolean isRejected = variantQuality < variantsDictionary.parameters.minQuality ||
                variantCoverage < variantsDictionary.parameters.minCoverage ||
                (!isHetCall && variantFrequency < variantsDictionary.parameters.minFrequency) ||
                (isHetCall && (variantFrequency < variantsDictionary.parameters.minHetFrequency ||
                        variantFrequency > variantsDictionary.parameters.maxHetFrequency));
        if (variantAlleleContent.length() == referenceAlleleContent.length() &&
                referenceAlleleContent.length() == 1) {
            // CASE: Nucleotide substitution.
            variantsDictionary
                    .addNucleotideVariant(featureEntry.name, variantPosition, variantAlleleContent, referenceAlleleContent,
                            sampleEntry.name, isPrimary);
        } else if (variantAlleleContent.length() > referenceAlleleContent.length() &&
                referenceAlleleContent.length() == 1) {
            // CASE: Insertion.
            variantsDictionary
                    .addNucleotideVariant(featureEntry.name, variantPosition, variantAlleleContent, referenceAlleleContent,
                            sampleEntry.name, isPrimary);
        } else if (variantAlleleContent.length() < referenceAlleleContent.length() &&
                variantAlleleContent.length() == 1) {
            // CASE: Deletion
            variantAlleleContent = variantAlleleContent + "-".repeat(referenceAlleleContent.length() - 1);
            variantsDictionary
                    .addNucleotideVariant(featureEntry.name, variantPosition, variantAlleleContent, referenceAlleleContent,
                            sampleEntry.name, isPrimary);
        } else {
            // CASE: Ambiguous call.
            Triplet<Integer, String, String> alignedAlleles = Bio.globalNucleotideSequenceAlignment(
                    referenceAllele.getBaseString(),
                    variantAllele.getBaseString(),
                    Bio.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FORBID,
                    Bio.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE,
                    null
            );
            ArrayList<String> resolvedVariants = Bio.getVariantsOfAlignedSequences(alignedAlleles.getValue1(), alignedAlleles.getValue2());
            for (String resolvedVariant : resolvedVariants) {
                int resolvedPosition = variantPosition + Integer.parseInt(resolvedVariant.split(VariantsDictionary.FIELD_SEPARATOR_1)[0]) - 1;
                String resolvedVariantAlleleContent = resolvedVariant.split(VariantsDictionary.FIELD_SEPARATOR_1)[1];
                String resolvedReferenceAlleleContent = resolvedVariant.split(VariantsDictionary.FIELD_SEPARATOR_1)[2];
                variantsDictionary
                        .addNucleotideVariant(featureEntry.name, resolvedPosition, resolvedVariantAlleleContent,
                                resolvedReferenceAlleleContent, sampleEntry.name, isPrimary);
                /* TODO: Handle resolved ambiguous variants. At least log to file.
                if (isPrimary && !isRejected) { ... }
                 */
            }
        }
    }
}
