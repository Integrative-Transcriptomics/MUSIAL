package utility;

import datastructure.Feature;
import datastructure.Form;
import datastructure.MusialStorage;
import datastructure.Sample;
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
    private static ArrayList<Triple<Integer, String, String>> variantSites;

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
                variantSites = new ArrayList<>();
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
     * Processes a single variant context.
     */
    @SuppressWarnings("DuplicatedCode")
    private static void processVariantContext() {
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
        if (genotypesContext.hasAD()) // (i) GATK HaplotypeCaller, freeBayes
            DP = Arrays.stream(genotypesContext.getAD()).sum();
        else if (variantContext.hasAttribute("DP4")) // (ii) bcftools
            DP = variantContext.getAttributeAsIntList("DP4", 0).stream().reduce(0, Integer::sum);
        else // (iii) Other
            DP = genotypesContext.getDP();
        Allele allele;
        String ref;
        String alt;
        for (int i = 0; i < alleles.size(); i++) { // Iterate over alleles, first one represents reference.
            allele = alleles.get(i);

            // Determine allele frequency: TODO: Extend for other input formats.
            if (genotypesContext.hasAD()) // (i) GATK HaplotypeCaller, freebayes
                AD = genotypesContext.getAD()[i];
            else if (variantContext.hasAttribute("DP4")) { // (ii) bcftools
                List<Integer> dp4 = variantContext.getAttributeAsIntList("DP4", 0);
                if (i == 0)
                    AD = dp4.get(0) + dp4.get(1); // Ref. call reads, i.e., reference allele coverage.
                else if (i == 1)
                    AD = dp4.get(2) + dp4.get(3); // Alt. call reads, i.e., alternative allele coverage.
                else {
                    Logger.logWarning("Inference of allele frequency by `DP4` does not support more than one ALT call. Skipping variant at position " + variantPosition + " of sample " + _sample.name + ".");
                    return;
                }
            } else {
                Logger.logWarning("Inference of allele frequency without attributes [AD, DP4] is not supported. Skipping variant at position " + variantPosition + " of sample " + _sample.name + ".");
                return;
            }
            AF = (double) AD / DP;

            // Store reference call and AD, AF values - only used to resolve ambiguous calls.
            if (i == 0) {
                refAF = AF;
                refAD = AD;
            }

            // Access reference and alternative allele content.
            ref = referenceAllele.getBaseString();
            alt = allele.getBaseString();

            // Decide on variant type and add content to all position calls.
            if (i == 0) {
                // 1. Allele is reference.
                callsPerPosition.putIfAbsent(variantPosition, new ArrayList<>());
                callsPerPosition.get(variantPosition).add(Quartet.with(ref, ref, refAF, refAD));
            } else if (ref.length() == 1 && alt.length() == 1) {
                // 2. Allele/variant is Substitution.
                callsPerPosition.putIfAbsent(variantPosition, new ArrayList<>());
                callsPerPosition.get(variantPosition).add(Quartet.with(ref, alt, AF, AD));
            } else if (alt.length() == 1) {
                // 3. Allele/variant is (simple) deletion.
                alt = alt + Constants.DELETION_OR_GAP_STRING.repeat(ref.length() - 1);
                callsPerPosition.putIfAbsent(variantPosition, new ArrayList<>());
                callsPerPosition.get(variantPosition).add(Quartet.with(ref, alt, AF, AD));
            } else if (ref.length() == 1) {
                // 4. Allele/variant is (simple) insertion.
                ref = ref + Constants.DELETION_OR_GAP_STRING.repeat(alt.length() - 1);
                callsPerPosition.putIfAbsent(variantPosition, new ArrayList<>());
                callsPerPosition.get(variantPosition).add(Quartet.with(ref, alt, AF, AD));
            } else {
                // 5. Ambiguous call/Complex InDel: Mutual variants are derived by aligning alternative allele sequence with reference allele sequence.
                Triplet<Integer, String, String> alignment = SequenceOperations.globalNucleotideSequenceAlignment(ref, alt, SequenceOperations.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FORBID, SequenceOperations.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE, null);
                ArrayList<Triple<Integer, String, String>> resolvedVariants = SequenceOperations.getVariantsOfAlignedSequences(alignment.getValue1(), alignment.getValue2());

                Logger.logWarning("Resolved ambiguous variant for sample " + _sample.name + " at position " + variantPosition + ": " + ref + ", " + alt + "\n" + resolvedVariants.stream().map(t -> " ".repeat(30) + t.getLeft() + "\t" + t.getMiddle() + "\t" + t.getRight()).collect(Collectors.joining("\n")));

                for (Triple<Integer, String, String> resolvedVariant : resolvedVariants) {
                    ref = resolvedVariant.getMiddle();
                    alt = resolvedVariant.getRight();
                    callsPerPosition.putIfAbsent(variantPosition + (resolvedVariant.getLeft()), new ArrayList<>());
                    // TODO: This may overwrite previous allocations.
                    if (callsPerPosition.get(variantPosition + (resolvedVariant.getLeft())).size() == 0) // Add reference only as first call, if not present.
                        callsPerPosition.get(variantPosition + (resolvedVariant.getLeft())).add(Quartet.with(ref, ref, refAF, refAD));
                    callsPerPosition.get(variantPosition + (resolvedVariant.getLeft())).add(Quartet.with(ref, alt, AF, AD));
                }
            }
        }

        // Decide on most frequent call and transfer to sample and form.
        ArrayList<Quartet<String, String, Double, Integer>> calls;
        for (Integer position : callsPerPosition.keySet()) {
            calls = callsPerPosition.get(position);

            // Skip position, if only one call (reference) is stored.
            if (calls.size() == 1)
                continue;

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
            if (!Objects.equals(mostFrequentCallIndex, "0")) { // If most frequent call is not reference allele.
                ref = mostFrequentCall.getValue0();
                alt = mostFrequentCall.getValue1();
                if (mostFrequentCallIndex.equals(Constants.CALL_INFO_REJECTED)) { // If call does not pass filters, add information as any nucleotide.
                    int variantLength = Math.max(ref.length(), alt.length());
                    alt = Constants.ANY_NUCLEOTIDE_STRING.repeat(variantLength);
                }
                variantSites.add(Triple.of(position, ref, alt));
            }
            _sample.addVariantCall(feature.name, position, new Tuple<>(mostFrequentCallIndex, DP), calls.stream().map(c -> new Tuple<>(c.getValue1(), c.getValue3())).collect(Collectors.toList()));
        }
    }

    /**
     * Allocates variant information for the analyzed sample.
     */
    @SuppressWarnings("DuplicatedCode")
    private static void allocateVariantsInformation() {
        StringBuilder acceptedSitesString = new StringBuilder(variantSites.size() * 4);
        String formName;
        if (variantSites.size() == 0) {
            formName = Constants.REFERENCE_FORM_NAME;
        } else {
            int pos;
            String ref;
            String alt;
            int count_substitution = 0;
            int count_insertion = 0;
            int count_deletion = 0;
            int count_ambiguous = 0;
            for (Triple<Integer, String, String> s : variantSites) {
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
