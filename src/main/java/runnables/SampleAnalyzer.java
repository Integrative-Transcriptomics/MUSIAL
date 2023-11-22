package runnables;

import datastructure.*;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import main.MusialConstants;
import org.javatuples.Triplet;
import utility.Bio;
import utility.Compression;
import utility.Logger;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentSkipListMap;

import static org.forester.util.ForesterUtil.round;

/**
 * Implementation of the {@link Runnable} interface to execute threaded sample analysis.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.0
 */
public final class SampleAnalyzer implements Runnable {

    /**
     * The {@link Sample} to be analyzed.
     */
    private final Sample sample;
    /**
     * The names of the {@link Feature}s to be analyzed.
     */
    private final Iterator<String> featureNames;
    /**
     * The current {@link Feature} to be analyzed.
     */
    private Feature feature;
    /**
     * The {@link MusialStorage} instance the sample and feature originate from.
     */
    private final MusialStorage musialStorage;
    /**
     * Concatenation of all accepted variants of rank one to infer alleles (see {@link Form}).
     */
    private StringBuilder variantFingerprint = new StringBuilder();
    /**
     * A {@link ConcurrentSkipListMap} to temporarily store novel variant calls.
     * <p>
     * Entries are stored as _:_ separated {@link String} values in the following order:
     * <ul>
     *  <li>position of the variant wrt. the reference genome</li>
     *  <li>rank of the variant in case of multiple calls of the same sample at the same position</li>
     *  <li>if the variant is rejected wrt. {@link MusialStorage#parameters}</li>
     *  <li>alt. base content of the call</li>
     *  <li>ref. base content of the call</li>
     *  <li>frequency of the call wrt. read support</li>
     *  <li>quality of the call</li>
     *  <li>coverage of the call wrt. read support</li>
     * </ul>
     */
    public HashSet<String> variantCalls = new HashSet<>();

    /**
     * Constructor of the {@link SampleAnalyzer} class.
     *
     * @param sample        The {@link Sample} to be analyzed.
     * @param featureNames  The names of the {@link Feature}s to be analyzed.
     * @param musialStorage The {@link MusialStorage} instance the sample and feature originate from.
     */
    public SampleAnalyzer(Sample sample, Iterator<String> featureNames,
                          MusialStorage musialStorage) {
        this.sample = sample;
        this.featureNames = featureNames;
        this.musialStorage = musialStorage;
    }

    /**
     * Run sample variant call analysis.
     */
    @Override
    public void run() {
        String featureName;
        while (this.featureNames.hasNext()) {
            featureName = featureNames.next();
            this.feature = musialStorage.getFeature(featureName);
            this.variantFingerprint = new StringBuilder();
            this.variantCalls = new HashSet<>();
            try {
                // Query variants located in the analyzed features' region.
                Iterator<VariantContext> variantContextIterator = sample.vcfFileReader.query(
                        feature.chromosome, feature.start, feature.end
                );
                while (variantContextIterator.hasNext()) {
                    VariantContext variantContext = variantContextIterator.next();
                    String variantContig = variantContext.getContig();
                    if (!variantContig.equals(feature.chromosome)) {
                        // If the variant is on a chromosome other than the one specified in the variants' dictionary it is skipped.
                        continue;
                    }
                    int variantPosition = variantContext.getStart();
                    if (musialStorage.excludedPositions.containsKey(variantContig)
                            && musialStorage.excludedPositions.get(variantContig).contains(variantPosition)) {
                        continue;
                    }
                    double variantQuality = round(variantContext.getPhredScaledQual(), 2);
                    double variantCoverage = variantContext.getAttributeAsDouble("DP", 0.0);
                    double variantFrequency;
                    Allele referenceAllele = variantContext.getReference();
                    List<Allele> alternateAlleles = variantContext.getAlternateAlleles();
                    Allele variantAllele;
                    if (alternateAlleles.size() != 0) {
                        // The allelic frequencies are extracted to order the alleles with respect to decreasing frequency.
                        // Frequencies are computed as AD / DP, not AF!
                        // Note: If multiple samples are present, the allelic depth will be summed. I.e., one biol. sample per file.
                        ArrayList<int[]> ADsPerGenotype = new ArrayList<>();
                        for (Genotype genotype : variantContext.getGenotypes()) {
                            int[] genotypeADs = genotype.getAD();
                            if (genotypeADs != null)
                                ADsPerGenotype.add(genotypeADs);
                        }
                        TreeMap<Double, Allele> sortedAlternateAlleles = new TreeMap<>(Collections.reverseOrder());
                        int ADSum = 0;
                        for (int i = 0; i < alternateAlleles.size(); i++) {
                            for (int[] ADs : ADsPerGenotype) {
                                ADSum += ADs[i + 1];
                            }
                            sortedAlternateAlleles.put(round(ADSum / variantCoverage, 2), alternateAlleles.get(i));
                        }
                        int rank = 1;
                        for (Map.Entry<Double, Allele> entry : sortedAlternateAlleles.entrySet()) {
                            variantFrequency = entry.getKey();
                            variantAllele = entry.getValue();
                            processVariantCall(variantPosition, variantQuality, variantCoverage, variantFrequency, variantAllele,
                                    referenceAllele, rank, alternateAlleles.size());
                            rank += 1;
                        }
                    }
                }

                transferVariantCalls();
            } catch (Exception e) {
                Logger.logError("An error occurred during the analysis of sample " + sample.name + " (file: " + sample.vcfFile.getAbsolutePath() + ") and feature " + featureName + "; " + e.getMessage());
                e.printStackTrace();
            }
        }
        sample.killVcfFileReader();
    }

    /**
     * Process a single variant call.
     */
    private void processVariantCall(int variantPosition, double variantQuality, double variantCoverage,
                                    double variantFrequency, Allele variantAllele, Allele referenceAllele,
                                    int rank, int maximalRank) {
        String variantAlleleContent = variantAllele.getBaseString();
        String referenceAlleleContent = referenceAllele.getBaseString();
        boolean heterozygous = maximalRank != 1;
        boolean rejected = variantQuality < musialStorage.parameters.minimalQuality ||
                variantCoverage < musialStorage.parameters.minimalCoverage ||
                (!heterozygous && variantFrequency < musialStorage.parameters.minimalHomozygousFrequency) ||
                (heterozygous && (variantFrequency < musialStorage.parameters.minimalHeterozygousFrequency ||
                        variantFrequency > musialStorage.parameters.maximalHeterozygousFrequency));
        if (variantAlleleContent.length() == referenceAlleleContent.length() &&
                referenceAlleleContent.length() == 1) {
            // CASE: Nucleotide substitution.
            addVariantCall(variantPosition,
                    rank,
                    maximalRank,
                    rejected,
                    variantAlleleContent,
                    referenceAlleleContent,
                    variantFrequency,
                    variantQuality,
                    variantCoverage
            );
        } else if (variantAlleleContent.length() > referenceAlleleContent.length() &&
                referenceAlleleContent.length() == 1) {
            // CASE: Insertion.
            addVariantCall(variantPosition,
                    rank,
                    maximalRank,
                    rejected,
                    variantAlleleContent,
                    referenceAlleleContent + "-".repeat(variantAlleleContent.length() - 1),
                    variantFrequency,
                    variantQuality,
                    variantCoverage
            );
        } else if (variantAlleleContent.length() < referenceAlleleContent.length() &&
                variantAlleleContent.length() == 1) {
            // CASE: Deletion
            addVariantCall(variantPosition,
                    rank,
                    maximalRank,
                    rejected,
                    variantAlleleContent + "-".repeat(referenceAlleleContent.length() - 1),
                    referenceAlleleContent,
                    variantFrequency,
                    variantQuality,
                    variantCoverage
            );
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
                int resolvedPosition = variantPosition + Integer.parseInt(resolvedVariant.split(MusialConstants.FIELD_SEPARATOR_1)[0]) - 1;
                String resolvedVariantAlleleContent = resolvedVariant.split(MusialConstants.FIELD_SEPARATOR_1)[1];
                String resolvedReferenceAlleleContent = resolvedVariant.split(MusialConstants.FIELD_SEPARATOR_1)[2];
                addVariantCall(resolvedPosition,
                        rank,
                        maximalRank,
                        rejected,
                        resolvedVariantAlleleContent,
                        resolvedReferenceAlleleContent,
                        variantFrequency,
                        variantQuality,
                        variantCoverage
                );
            }
        }
    }

    /**
     * Add information about a variant call to {@link #variantCalls}.
     */
    private void addVariantCall(int position, int rank, int maximalRank, boolean rejected, String alternativeContent, String referenceContent, double frequency, double quality, double coverage
    ) {
        this.variantCalls.add(position + MusialConstants.FIELD_SEPARATOR_1 +
                rank + "/" + maximalRank + MusialConstants.FIELD_SEPARATOR_1 +
                rejected + MusialConstants.FIELD_SEPARATOR_1 +
                alternativeContent + MusialConstants.FIELD_SEPARATOR_1 +
                referenceContent + MusialConstants.FIELD_SEPARATOR_1 +
                frequency + MusialConstants.FIELD_SEPARATOR_1 +
                quality + MusialConstants.FIELD_SEPARATOR_1 +
                coverage
        );
        if (rank == 1 && !rejected) {
            this.variantFingerprint.append(position).append(MusialConstants.FIELD_SEPARATOR_1).append(alternativeContent).append(MusialConstants.FIELD_SEPARATOR_2);
        }
    }

    /**
     * @return Iterator over all variant calls in {@link #variantCalls}.
     */
    private Iterator<String> getVariantCallsIterator() {
        return this.variantCalls.iterator();
    }

    /**
     * Transfer all collected variant calls to {@link #feature}.
     *
     * @throws IOException If variant encoding fails.
     */
    private void transferVariantCalls() throws IOException {
        String[] variantCall;
        if (this.variantFingerprint.length() != 0) {
            this.variantFingerprint.setLength(this.variantFingerprint.length() - 1);
        }
        // Construct allele name from variant fingerprint.
        String alleleName;
        if (this.variantFingerprint.length() == 0) {
            alleleName = MusialConstants.REFERENCE_ID;
        } else {
            alleleName = constructAlleleName( );
        }
        Form allele = new Form(alleleName);
        // Add default allele annotations.
        if (!alleleName.equals(MusialConstants.REFERENCE_ID)) {
            allele.addAnnotation(MusialConstants.VARIANTS, Compression.brotliEncodeString(this.variantFingerprint.toString()));
        } else {
            allele.addAnnotation(MusialConstants.VARIANTS, "");
        }
        this.feature.addAllele(allele);
        // Add sample information.
        this.feature.getAllele(alleleName).addOccurrence(sample.name);
        this.sample.addAnnotation(MusialConstants.SAMPLE_ANNOTATION_ALLELE_PREFIX + this.feature.name, allele.name);
        // Add nucleotide variant entries.
        for (Iterator<String> variantCallIterator = getVariantCallsIterator(); variantCallIterator.hasNext(); ) {
            variantCall = variantCallIterator.next().split(MusialConstants.FIELD_SEPARATOR_1);
            transferVariantCall(variantCall);
        }
    }

    /**
     * Transfer a single variant call to {@link #feature}.
     * <p>
     * See format of {@link #variantCalls} for field access documentation.
     */
    private void transferVariantCall(String[] variantCall) {
        feature.addNucleotideVariant(Integer.parseInt(variantCall[0]), variantCall[3]);
        VariantAnnotation variantAnnotation = feature.getNucleotideVariantAnnotation(Integer.parseInt(variantCall[0]), variantCall[3]);
        variantAnnotation.addProperty(MusialConstants.REFERENCE_CONTENT, variantCall[4]);
        variantAnnotation.addProperty(MusialConstants.VARIANT_OCCURRENCE_SAMPLE_PREFIX + sample.name, variantCall[1] + MusialConstants.FIELD_SEPARATOR_1 +
                variantCall[2] + MusialConstants.FIELD_SEPARATOR_1 +
                variantCall[5] + MusialConstants.FIELD_SEPARATOR_1 +
                variantCall[6] + MusialConstants.FIELD_SEPARATOR_1 +
                variantCall[7]
        );
    }

    /**
     * Generates a {@link String} representing the name of an allele.
     *
     * @return Allele name.
     */
    @SuppressWarnings("DuplicatedCode")
    private String constructAlleleName( ) {
        String variantFingerprint = this.variantFingerprint.toString();
        int variantFingerprintCode = variantFingerprint.hashCode();
        if ( musialStorage.formNameCodes.containsKey( variantFingerprintCode ) ) {
            return musialStorage.formNameCodes.get( variantFingerprintCode );
        } else {
            int s = 0;
            int i = 0;
            int d = 0;
            for (String variant : variantFingerprint.split(MusialConstants.FIELD_SEPARATOR_2)) {
                String[] variantInformation = variant.split(MusialConstants.FIELD_SEPARATOR_1);
                String alt = variantInformation[1];
                if (alt.contains("-")) {
                    d += 1;
                } else if (alt.length() > 1) {
                    i += 1;
                } else {
                    s += 1;
                }
            }
            String alleleName = ( feature.getAlleleCount() + 1 ) + ".S" + s + ".I" + i + ".D" + d;
            musialStorage.formNameCodes.put( variantFingerprintCode, alleleName );
            return alleleName;
        }
    }
}
