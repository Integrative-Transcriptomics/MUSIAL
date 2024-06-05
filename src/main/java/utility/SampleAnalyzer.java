package utility;

import datastructure.Feature;
import datastructure.Form;
import datastructure.MusialStorage;
import datastructure.Sample;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import main.Constants;
import org.apache.commons.lang3.tuple.Triple;
import org.javatuples.Quartet;
import org.javatuples.Triplet;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Implements static methods to run sample/.vcf analysis.
 *
 * @author Simon Hackl
 */
public final class SampleAnalyzer {

    /**
     * The MusialStorage instance containing genetic variant data.
     */
    private static MusialStorage storage;

    /**
     * The current Feature being analyzed.
     */
    private static Feature feature;

    /**
     * The current Sample being analyzed.
     */
    private static Sample _sample;

    /**
     * Container of all accepted variant sites to infer alleles.
     */
    private static ArrayList<Triple<Integer, String, String>> acceptedSites;

    /**
     * The VariantContext representing a genetic variant.
     */
    private static VariantContext variantContext;

    /**
     * The position of the genetic variant.
     */
    private static int variantPosition;

    /**
     * Runs the sample analysis for a given Sample, iterator of feature names, and MusialStorage.
     *
     * @param sample        The Sample instance to be analyzed.
     * @param featureNames  An iterator of feature names to be analyzed.
     * @param musialStorage The MusialStorage instance containing genetic variant data.
     */
    public static void run(Sample sample, Iterator<String> featureNames, MusialStorage musialStorage) {
        storage = musialStorage;
        _sample = sample;
        while (featureNames.hasNext()) {
            feature = storage.getFeature(featureNames.next());
            try {
                // Query variants located in the analyzed features' region.
                Iterator<VariantContext> variantContextIterator = _sample.vcfFileReader.query(
                        feature.contig, feature.start, feature.end
                );
                acceptedSites = new ArrayList<>();
                String contig;
                while (variantContextIterator.hasNext()) {
                    variantContext = variantContextIterator.next();
                    contig = variantContext.getContig();
                    if (!contig.equals(feature.contig))
                        // If the variant is on a contig other than the one specified in the variants' dictionary it is skipped.
                        continue;
                    variantPosition = variantContext.getStart();
                    if (storage.isPositionExcluded(contig, variantPosition))
                        // Skip excluded positions.
                        continue;
                    // Process the variant, i.e., filter and decide on actual content.
                    processVariantContext();
                }
                // Transfers all variant calls of the analyzed sample to the resp. feature.
                allocateVariantsInformation();
            } catch (Exception e) {
                Logger.logError("An error occurred during the analysis of sample " + _sample.name + " (file: " + _sample.vcfFile.getAbsolutePath() + ") and feature " + feature.name + " at position " + variantPosition + "; " + e.getMessage());
                e.printStackTrace();
            }
        }
        _sample.killVcfFileReader();
    }

    /**
     * Processes a single variant call within the sample analysis.
     */
    @SuppressWarnings("DuplicatedCode")
    private static void processVariantContext() throws MusialException {
        // Representation of variant calls as: ( Ref, Alt, AF, AD )
        LinkedHashMap<Integer, ArrayList<Quartet<String, String, Double, Integer>>> callsPerPosition = new LinkedHashMap<>();
        Allele referenceAllele = variantContext.getReference();
        ArrayList<Allele> alleles = new ArrayList<>();
        alleles.add(referenceAllele);
        alleles.addAll(variantContext.getAlternateAlleles());
        Genotype genotypesContext = variantContext.getGenotypes().get(0); // Only one genotype is expected, i.e., single sample.
        // TODO: Change for multi sample .vcf file support.
        int DP; // Total depth of coverage at site. Filtered reads are not counted.
        int AD; // Allelic depth of coverage, i.e., read support for one specific allele.
        double AF; // Allelic frequency, i.e., frequency of one specific allele wrt. all reads.
        int refAD = -1;
        double refAF = -1;
        // Determine depth of coverage: TODO: Extend for other input formats.
        if (genotypesContext.hasAD()) // (1) GATK HaplotypeCaller and (2) freeBayes
            DP = Arrays.stream(genotypesContext.getAD()).sum();
        else if (variantContext.hasAttribute("DP4")) // (3) bcftools
            DP = variantContext.getAttributeAsIntList("DP4", 0).stream().reduce(0, Integer::sum);
        else // (3) Other
            DP = genotypesContext.getDP();
        Allele allele;
        String ref;
        String alt;
        for (int i = 0; i < alleles.size(); i++) {
            allele = alleles.get(i);
            // Determine allele frequency: TODO: Extend for other input formats.
            if (genotypesContext.hasAD()) // (1) GATK HaplotypeCaller
                AD = genotypesContext.getAD()[i];
            else if (variantContext.hasAttribute("DP4")) { // (2) bcftools
                List<Integer> dp4 = variantContext.getAttributeAsIntList("DP4", 0);
                if (i == 0)
                    AD = dp4.get(0) + dp4.get(1); // Ref. call reads, i.e., reference allele coverage.
                else if (i == 1)
                    AD = dp4.get(2) + dp4.get(3); // Alt. call reads, i.e., alternative allele coverage.
                else {
                    // TODO: Get GT?
                    Logger.logWarning("Inference of allele frequency by `DP4` does not support more than one ALT call. Skipping variant at position " + variantPosition + " of sample " + _sample.name + ".");
                    return;
                }
            } else {
                throw new MusialException("Failed to infer allele frequency; Variant has none of the attributes [AD, DP4]!");
            }
            AF = (double) AD / DP;
            if (i == 0) { // Store reference call values for use in ambiguous call resolution.
                refAD = AD;
                refAF = AF;
            }
            ref = referenceAllele.getBaseString();
            alt = allele.getBaseString();
            // Adjust base strings.
            if (alt.length() == 1 && ref.length() > alt.length()) {
                // 1. Allele represents deletion.
                alt = alt + Constants.DELETION_OR_GAP_STRING.repeat(ref.length() - 1);
                callsPerPosition.putIfAbsent(variantPosition, new ArrayList<>());
                callsPerPosition.get(variantPosition).add(Quartet.with(ref, alt, AF, AD));
                if (alt.charAt(0) != ref.charAt(0)) // Temp. warning for unexpected variant call format.
                    Logger.logWarning("Unexpected variant content for sample " + _sample.name + " at position " + variantPosition + ": " + ref + " " + alt);
            } else if (ref.length() == 1 && alt.length() > ref.length()) {
                // 2. Allele represents insertion.
                ref = ref + Constants.DELETION_OR_GAP_STRING.repeat(alt.length() - 1);
                callsPerPosition.putIfAbsent(variantPosition, new ArrayList<>());
                callsPerPosition.get(variantPosition).add(Quartet.with(ref, alt, AF, AD));
                if (alt.charAt(0) != ref.charAt(0)) // Temp. warning for unexpected variant call format.
                    Logger.logWarning("Unexpected variant content for sample " + _sample.name + " at position " + variantPosition + ": " + ref + " " + alt);
            } else if (ref.length() == 1 && alt.length() == 1) {
                // 3. Substitution.
                callsPerPosition.putIfAbsent(variantPosition, new ArrayList<>());
                callsPerPosition.get(variantPosition).add(Quartet.with(ref, alt, AF, AD));
            } else if (ref.length() > 1 && alt.length() > 1) {
                // 4. Ambiguous call: Mutual variants are derived by aligning alternative allele with reference sequence.
                if (i != 0) { // Skip reference allele to avoid variants that contain itself as ref. and alt.
                    Triplet<Integer, String, String> alignment = SequenceOperations.globalNucleotideSequenceAlignment(ref, alt, SequenceOperations.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FORBID, SequenceOperations.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE, null);
                    ArrayList<Triple<Integer, String, String>> resolvedVariants = SequenceOperations.getVariantsOfAlignedSequences(alignment.getValue1(), alignment.getValue2());
                    for (Triple<Integer, String, String> resolvedVariant : resolvedVariants) {
                        callsPerPosition.putIfAbsent(variantPosition + (resolvedVariant.getLeft()), new ArrayList<>());
                        // TODO: Validate correctness for highly complex variants, e.g. (ref) CACCC (alt) GTCCG, ATCCT.
                        if (callsPerPosition.get(variantPosition + (resolvedVariant.getLeft())).size() == 0) // Add reference only as first call, if not present.
                            callsPerPosition.get(variantPosition + (resolvedVariant.getLeft())).add(Quartet.with(resolvedVariant.getMiddle(), resolvedVariant.getMiddle(), refAF, refAD));
                        callsPerPosition.get(variantPosition + (resolvedVariant.getLeft())).add(Quartet.with(resolvedVariant.getMiddle(), resolvedVariant.getRight(), AF, AD));
                    }
                }
            } else {
                // 5. Unexpected.
                Logger.logWarning("Ignore not processable variant content for sample " + _sample.name + " at position " + variantPosition + ": " + ref + " " + alt);
            }
        }
        // Decide on most frequent call and transfer to sample and form.
        ArrayList<Quartet<String, String, Double, Integer>> calls;
        for (Integer position : callsPerPosition.keySet()) {
            calls = callsPerPosition.get(position);
            //noinspection OptionalGetWithoutIsPresent
            Quartet<String, String, Double, Integer> mostFrequentCall = calls.stream().max(Comparator.comparingDouble(Quartet::getValue2)).get();
            String mostFrequentCallIndex;
            if (DP == 0) {
                mostFrequentCallIndex = Constants.CALL_INFO_NO_INFO;
            } else if (DP < storage.parameters.minimalCoverage || mostFrequentCall.getValue2() < storage.parameters.minimalFrequency) {
                // Conservative wrt. reference.
                mostFrequentCallIndex = Constants.CALL_INFO_REJECTED;
            } else {
                mostFrequentCallIndex = String.valueOf(calls.indexOf(mostFrequentCall));
            }
            if (!Objects.equals(mostFrequentCallIndex, "0")) {
                ref = mostFrequentCall.getValue0();
                alt = mostFrequentCall.getValue1();
                if (mostFrequentCallIndex.equals(Constants.CALL_INFO_NO_VARIANT) || mostFrequentCallIndex.equals(Constants.CALL_INFO_REJECTED)) {
                    int variantLength = Math.max(ref.length(), alt.length());
                    alt = Constants.ANY_NUCLEOTIDE_STRING.repeat(variantLength);
                }
                acceptedSites.add(Triple.of(position, ref, alt));
            }
            _sample.addVariantCall(feature.name, position, new Tuple<>(mostFrequentCallIndex, DP), calls.stream().map(c -> new Tuple<>(c.getValue1(), c.getValue3())).collect(Collectors.toList()));
        }
    }

    /**
     * Allocates variant information for the analyzed sample.
     */
    @SuppressWarnings("DuplicatedCode")
    private static void allocateVariantsInformation() {
        StringBuilder acceptedSitesString = new StringBuilder(acceptedSites.size() * 4);
        String formName;
        if (acceptedSites.size() == 0) {
            formName = Constants.REFERENCE_FORM_NAME;
        } else {
            int pos;
            String ref;
            String alt;
            int count_substitution = 0;
            int count_insertion = 0;
            int count_deletion = 0;
            int count_ambiguous = 0;
            for (Triple<Integer, String, String> s : acceptedSites) {
                pos = s.getLeft();
                ref = s.getMiddle();
                alt = s.getRight();
                if (alt.contains(Constants.ANY_NUCLEOTIDE_STRING)) {
                    count_ambiguous += alt.codePoints().filter(c -> c == Constants.ANY_NUCLEOTiDE_CHAR).count();
                } else if (ref.length() == 1 && alt.length() == 1) {
                    count_substitution += 1;
                } else {
                    count_insertion += ref.codePoints().filter(c -> c == '-').count();
                    count_deletion += alt.codePoints().filter(c -> c == '-').count();
                }
                acceptedSitesString.append(pos).append(Constants.FIELD_SEPARATOR_1).append(alt).append(Constants.FIELD_SEPARATOR_2);
                feature.addNucleotideVariant(pos, alt, ref);
                feature.getNucleotideVariant(pos, alt).addOccurrence(_sample.name);
            }
            acceptedSitesString.setLength(acceptedSitesString.length() - 1); // Delete last ';' separator symbol.
            formName = storage.getFormName(
                    acceptedSitesString.toString().hashCode(),
                    "A" + (feature.getAlleleCount() + 1) + ".s" + count_substitution + ".i" + count_insertion + ".d" + count_deletion + ".a" + count_ambiguous
            );
        }
        Form allele = new Form(formName, acceptedSitesString.toString());
        feature.addAllele(allele);
        feature.getAllele(formName).addOccurrence(_sample.name);
        _sample.setAllele(feature.name, allele.name);
    }

}
