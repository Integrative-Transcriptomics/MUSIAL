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
import org.javatuples.Quintet;
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
     * Container of variant sites to infer alleles. Stored as (Position, Reference, Alternative, Primary, Pass).
     */
    private static ArrayList<Quintet<Integer, String, String, Boolean, Boolean>> variantSites;

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
        // Check if all variants are excluded/at least one variant alleles remains.
        alleles = (ArrayList<Allele>) alleles.stream().filter(
                (a) -> !storage.isVariantExcluded(variantContext.getContig(), variantPosition, a.getBaseString())
        ).collect(Collectors.toList());
        if (alleles.size() < 2)
            return;
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
            // Access reference and alternative allele content.
            ref = referenceAllele.getBaseString();
            alt = allele.getBaseString();
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
                    if (callsPerPosition.get(variantPosition + (resolvedVariant.getLeft())).size() == 0) // Add reference only as first call, if not present.
                        callsPerPosition.get(variantPosition + (resolvedVariant.getLeft())).add(Quartet.with(ref, ref, refAF, refAD));
                    callsPerPosition.get(variantPosition + (resolvedVariant.getLeft())).add(Quartet.with(ref, alt, AF, AD));
                }
            }
        }
        // Decide on most frequent call and transfer to sample and form.
        int position;
        ArrayList<Quartet<String, String, Double, Integer>> calls;
        for (Map.Entry<Integer, ArrayList<Quartet<String, String, Double, Integer>>> entry : callsPerPosition.entrySet()) {
            position = entry.getKey();
            calls = entry.getValue();
            // Skip position, if only one call (reference) is stored.
            if (calls.size() == 1)
                continue;
            //noinspection OptionalGetWithoutIsPresent
            Quartet<String, String, Double, Integer> primaryCall = calls.stream().max(Comparator.comparingDouble(Quartet::getValue2)).get();
            int primaryCallIndex = calls.indexOf(primaryCall);
            boolean primary;
            boolean pass;
            boolean ambiguous = false;
            int callIndex;
            for (Quartet<String, String, Double, Integer> call : calls) {
                callIndex = calls.indexOf(call);
                primary = callIndex == primaryCallIndex; // Set if call is the primary call.
                pass = call.getValue3() >= storage.parameters.minimalCoverage || call.getValue2() >= storage.parameters.minimalFrequency; // Set if call passes acceptance criteria.
                ambiguous = primary && !pass;
                ref = call.getValue0();
                alt = call.getValue1();
                variantSites.add(Quintet.with(position, ref, alt, primary, pass));
            }
            // Add call information to sample.
            _sample.addVariantCall(
                    feature.name,
                    position,
                    new Tuple<>(
                            ambiguous ? Constants.CALL_INFO_REJECTED : String.valueOf(primaryCallIndex),
                            DP
                    ),
                    calls.stream().map(c -> new Tuple<>(c.getValue1(), c.getValue3())).collect(Collectors.toList())
            );
        }
    }

    /**
     * Allocates variant information for the analyzed sample.
     */
    @SuppressWarnings("DuplicatedCode")
    private static void allocateVariantsInformation() {
        ArrayList<String> primaryVariantsList = new ArrayList<>(variantSites.size());
        String formName;
        int pos;
        String ref;
        String alt;
        boolean isReference;
        boolean primary;
        boolean pass;
        for (Quintet<Integer, String, String, Boolean, Boolean> s : variantSites) {
            pos = s.getValue0();
            ref = s.getValue1();
            alt = s.getValue2();
            isReference = ref.equals(alt);
            primary = s.getValue3();
            pass = s.getValue4();
            if (primary) {
                if (!pass) {
                    // If the primary variant does not pass filter criteria, replace alt content with N to indicate ambiguity.
                    // This might also be the case, if the reference allele is the primary variant.
                    primaryVariantsList.add(pos + Constants.FIELD_SEPARATOR + Constants.ANY_NUCLEOTIDE_STRING);
                    feature.addNucleotideVariant(pos, Constants.ANY_NUCLEOTIDE_STRING, ref);
                    feature.getNucleotideVariant(pos, Constants.ANY_NUCLEOTIDE_STRING).addOccurrence(_sample.name);
                    //noinspection ConstantConditions
                    feature.getNucleotideVariant(pos, Constants.ANY_NUCLEOTIDE_STRING).addInfo(Constants.VARIANT_INFO_PRIMARY, String.valueOf(primary));
                    feature.getNucleotideVariant(pos, Constants.ANY_NUCLEOTIDE_STRING).addInfo(Constants.VARIANT_INFO_ACTUAL_ALT, alt);
                } else if (!isReference) {
                    // If the primary variant passes filter criteria AND is not the reference allele, add information as is.
                    primaryVariantsList.add(pos + Constants.FIELD_SEPARATOR + alt);
                    feature.addNucleotideVariant(pos, alt, ref);
                    feature.getNucleotideVariant(pos, alt).addOccurrence(_sample.name);
                    //noinspection ConstantConditions
                    feature.getNucleotideVariant(pos, alt).addInfo(Constants.VARIANT_INFO_PRIMARY, String.valueOf(primary));
                }
                // Otherwise the information is ignored, i.e., the reference allele is the primary variant.
            }
            if (!primary && pass && !isReference) {
                // If the variant is not primary and not reference, but passes the filter criteria, add information to feature anyway.
                feature.addNucleotideVariant(pos, alt, ref);
                feature.getNucleotideVariant(pos, alt).addOccurrence(_sample.name);
                //noinspection ConstantConditions
                feature.getNucleotideVariant(pos, alt).addInfo(Constants.VARIANT_INFO_PRIMARY, String.valueOf(primary));
            }
        }
        String primaryVariantsString = String.join(Constants.ENTRY_SEPARATOR, primaryVariantsList);
        if (primaryVariantsList.size() == 0)
            formName = Constants.REFERENCE_FORM_NAME;
        else
            formName = storage.getFormName(
                    primaryVariantsString.hashCode(),
                    "A" + (feature.getAlleleCount() + 1)
            );
        Form allele = new Form(formName, primaryVariantsString);
        feature.addAllele(allele);
        feature.getAllele(formName).addOccurrence(_sample.name);
        _sample.setAllele(feature.name, allele.name);
    }

}
