package op;

import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import model.*;
import util.Bio;
import util.Constants;
import util.IO;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Encapsulates complex update logic for genomic data storage from the {@link Storage} class itself.
 * <p>
 * This class provides methods to update sample attributes, process variant calls from {@link Sample} instances into {@link Variant}s in
 * stored {@link Contig}s, calculate sequence types, and generate statistics.
 * <p>
 * The single update methods should be called only once the relevant data is loaded into the {@link Storage} instance, e.g.
 * {@link #updateSequenceTypes()} without prior processing of VCF files will not have any effect.
 */
public class StorageUpdater {

    /**
     * The storage object for managing genomic data.
     */
    private final Storage storage;

    /**
     * Constructs a new StorageUpdater instance.
     * <p>
     * This constructor initializes the StorageUpdater with the provided storage object, which is used to manage genomic data and perform
     * updates on samples, features, and contigs.
     *
     * @param storage The storage object that contains genomic data to be updated.
     */
    public StorageUpdater(Storage storage) {
        this.storage = storage;
    }

    /**
     * Updates the attributes of samples in the storage.
     * <p>
     * This method iterates through the provided map of attributes, where each entry consists of a sample identifier and a map of
     * attributes. If the sample exists in the storage, it adds the attributes to the sample only if they are not already present.
     * <p>
     * <i>This method should only ba called after {@link VCFProcessor#processFiles()} was called in the context of the
     * {@code build} or {@code expand} tasks.</i>
     *
     * @param attributes A map where the key is the sample identifier, and the value is another map containing attribute key-value pairs to
     *                   be added to the sample.
     */
    public void updateSampleAttributes(Map<String, Map<String, String>> attributes) {
        for (var entry : attributes.entrySet()) {
            String sampleIdentifier = entry.getKey();
            if (storage.hasSample(sampleIdentifier)) {
                storage.getSample(sampleIdentifier).setAttributesIfAbsent(entry.getValue());
            }
        }
    }

    /**
     * Updates the sequence types for all features and samples in the storage.
     * <p>
     * This method performs the following tasks:
     * <ul>
     *   <li>Iterates through all features and samples in the storage.</li>
     *   <li>Updates alleles for each feature and sample based on variants.</li>
     *   <li>Calculates sequence effects such as frameshifts and updates allele attributes.</li>
     *   <li>If the feature is coding and the contig has a sequence, updates proteoforms.</li>
     *   <li>Handles proteoform sequence alignment, variant extraction, and effect annotation.</li>
     * </ul>
     * <p>
     * <i>This method should only ba called after {@link VCFProcessor#processFiles()} and a respective annotation method
     * implemented in {@link VariantAnnotator} was called in the context of the {@code build} or {@code expand} tasks.</i>
     *
     * @throws IOException     If an I/O error occurs during processing.
     * @throws MusialException If a specific error related to the Musial library occurs.
     */
    public void updateSequenceTypes() throws IOException, MusialException {
        // Iterate through all features in the storage.
        for (Feature feature : storage.getFeatures()) {
            // Retrieve the contig associated with the feature.
            Contig contig = storage.getContig(feature.contig);

            // Iterate through all active samples to infer allele sequence types.
            for (Sample sample : storage.getActiveSamples()) {
                // Update allele.
                Allele allele;
                String alleleIdentifier = sample.getRelatedAllele(feature._id);
                List<Variant> variants;

                // Note: In this special case, a returned reference allele may also indicate the presence of a new allele.
                if (alleleIdentifier.equals(Constants.REFERENCE)) {
                    // Retrieve variants for the sample within the feature's range.
                    variants = contig.getVariantsOfSamplesWithin(feature.start, feature.end, Collections.singleton(sample._id));

                    // Construct stubs and collect effects.
                    List<Variant.Stub> stubs = new ArrayList<>(variants.size());
                    Set<String> effects = new HashSet<>();
                    for (Variant variant : variants) {
                        if (!variant.isFiltered(sample._id)) {
                            stubs.add(variant.asStub());
                            effects.addAll(variant.getAttributeSet(Constants.SNP_EFF_PREFIX + Constants.SNP_EFF_KEYS.get(1)));
                        } else if (storage.parameters.maskFiltered()) {
                            stubs.add(variant.asMaskedStub());
                        }
                    }

                    // If no variants/stubs are present, continue to the next sample.
                    if (stubs.isEmpty()) continue;

                    // Determine the allele based on the stubs.
                    allele = new Allele(stubs);

                    // Check if the allele already exists for the feature, i.e., was created by another sample.
                    if (feature.hasAllele(allele._id)) {
                        allele = feature.getAllele(allele._id);
                    } else { // If the allele is new, annotate its effects and add relations.
                        // Calculate sequence length deviation and frameshift effects.
                        int lengthDelta = Integer.parseInt(allele.getAttribute(Constants.AttributesKeys.SEQUENCE_LENGTH_DEVIATION));
                        int netFrameshift = Math.abs(lengthDelta % 3);

                        // Add frameshift effects if applicable.
                        if (netFrameshift != 0) {
                            effects.add(lengthDelta > 0 ? "plus_%d_frameshift".formatted(netFrameshift)
                                    : "minus_%d_frameshift".formatted(netFrameshift));
                        }
                        allele.setAttribute(Constants.AttributesKeys.SO_EFFECTS, String.join(Constants.COMMA, effects));

                        // Add relations between variants and the allele.
                        for (Variant.Stub stub : stubs) {
                            if (stub.alternative().equals(Constants.ANY_NUCLEOTIDE)) continue;
                            contig.getVariant(stub.position(), stub.alternative()).addAlleleRelation(feature._id, allele._id);
                        }

                        // Add the new allele to the feature.
                        feature.addAllele(allele);
                    }

                    // Add relations between the allele, sample, and feature.
                    allele.addSampleRelation(sample._id);
                    sample.addRelation(feature._id, allele._id);
                }
            }

            // If the feature is coding and the contig has a sequence, update proteoform sequence types.
            if (feature.isCoding() && contig.hasSequence()) {

                // Iterate through all alleles of the feature to infer proteoform sequence types.
                for (Allele allele : feature.getAlleles()) {
                    Proteoform proteoform;

                    // Skip if a proteoform is already assigned.
                    String proteoformIdentifier = allele.getRelatedProteoform();
                    if (proteoformIdentifier != null) continue;

                    // Construct a position-sorted map of variants/stubs of the allele.
                    NavigableMap<Integer, String> variants = new TreeMap<>();
                    allele.getStubs().forEach(stub -> variants.put(stub.position() - feature.start, stub.alternative()));

                    // Translate the reference and allele sequences.
                    String referenceNucleotideSequence = contig.getSequence(feature.start, feature.end);
                    String referenceAminoacidSequence = Bio.translateSequence(referenceNucleotideSequence, feature.isReverse());
                    String proteoformSequence = Bio.translateSequence(Bio.integrateVariants(referenceNucleotideSequence, variants, true),
                            feature.isReverse());
                    assert proteoformSequence != null;

                    // Check if the proteoform sequence is synonymous with the reference sequence.
                    if (referenceAminoacidSequence.equals(proteoformSequence)) {
                        proteoformIdentifier = Constants.SYNONYMOUS;
                    } else {
                        // Align sequences and extract amino acid variants.
                        Tuple<String, String> alignment = Bio.globalProteinSequenceAlignment(
                                referenceAminoacidSequence, proteoformSequence, Math.max(referenceAminoacidSequence.length(),
                                        proteoformSequence.length()),
                                6, true, false, Math.abs(referenceAminoacidSequence.length() - proteoformSequence.length())
                        );

                        List<Variant.Stub> aaVariants = Bio.getCanonicalVariants(alignment.a, alignment.b).stream()
                                .map(aav -> new Variant.Stub(aav.getLeft() + 1, aav.getRight()))
                                .collect(Collectors.toList());
                        if (aaVariants.isEmpty()) throw new IllegalArgumentException("Proteoform has no variants.");

                        // Create a new proteoform or retrieve an existing one.
                        proteoform = feature.hasProteoform(new Proteoform(aaVariants)._id)
                                ? feature.getProteoform(new Proteoform(aaVariants)._id)
                                : new Proteoform(aaVariants);

                        // Annotate proteoform effects if it is newly created.
                        if (!feature.hasProteoform(proteoform._id)) {
                            Set<String> effects = new HashSet<>();
                            String alleleEffects = allele.getAttribute(Constants.AttributesKeys.SO_EFFECTS);

                            // Add specific effects based on the proteoform's variants.
                            if (alleleEffects.contains("frameshift")) {
                                effects.add("frameshift_sequence_variation");
                            }
                            if (aaVariants.stream().anyMatch(aav -> Bio.isInsertion(aav.alternative()))) {
                                effects.add("amino_acid_insertion");
                            }
                            if (aaVariants.stream().anyMatch(aav -> Bio.isDeletion(aav.alternative()))) {
                                effects.add("amino_acid_deletion");
                            }
                            if (aaVariants.stream().anyMatch(aav -> Bio.isSubstitution(aav.alternative()))) {
                                effects.add("amino_acid_substitution");
                            }
                            if (proteoform.hasVariantAt(1) && proteoform.getVariant(1).charAt(0) != referenceAminoacidSequence.charAt(0)) {
                                effects.add("start_lost");
                            }
                            aaVariants.stream()
                                    .filter(aav -> aav.alternative().contains(Constants.TERMINAL_AA))
                                    .findFirst()
                                    .ifPresent(stopCodonVariant -> {
                                        int stopCodonPosition = stopCodonVariant.position() +
                                                stopCodonVariant.alternative().indexOf(Constants.TERMINAL_AA);
                                        String effect = stopCodonPosition <= referenceAminoacidSequence.length() ? "stop_gained" :
                                                "redundant_inserted_stop_gained";
                                        effects.add(effect);
                                    });
                            proteoform.setAttribute(Constants.AttributesKeys.SO_EFFECTS, String.join(Constants.COMMA, effects));
                            feature.addProteoform(proteoform);
                        }

                        // Add the new proteoform to the feature and establish relations.
                        proteoformIdentifier = proteoform._id;
                        proteoform.addRelation(allele._id);
                    }

                    // Set the proteoform identifier for the allele.
                    allele.setProteoformRelation(proteoformIdentifier);
                }
            }
        }
    }

    /**
     * Updates statistical attributes for samples, contigs, and features in the storage.
     * <p>
     * This method calculates and updates various statistics, including:
     * <ul>
     *   <li>Number of calls, filtered calls, mean coverage, and mean entropy for each sample.</li>
     *   <li>Frequency of reference alleles and disrupted coding features for each sample.</li>
     *   <li>Variant frequencies for each contig and sample-specific variant counts.</li>
     *   <li>Allelic frequencies, proteoform frequencies, and disrupted proteoform frequencies for each feature.</li>
     * </ul>
     * <p>
     * <i>This method should only ba called as the last step before serializing a {@link Storage} instance.</i>
     */
    public void updateStatistics() {
        // Calculate the total number of features, coding features, and samples.
        int noFeatures = storage.getFeatures().size();
        int noCodingFeatures = (int) storage.getFeatures().stream().filter(Feature::isCoding).count();
        int noSamples = storage.getSamples().size();

        // Update variant statistics.
        statisticsOfVariants(noSamples);

        // Update sample statistics.
        statisticsOfSamples(noFeatures, noCodingFeatures);

        // Update feature statistics.
        statisticsOfFeatures(noSamples);
    }

    /**
     * Updates statistical attributes in the {@link Storage} from iterating through all {@link Contig#variants} in all
     * {@link Storage#contigs}.
     * <p>
     * This method calculates and updates various statistics for variants and samples, including:
     * <ul>
     *   <li>Variant frequency based on the number of related samples.</li>
     *   <li>Per-sample counts of SNVs, insertions/deletions, and filtered calls.</li>
     *   <li>Mean coverage and entropy for each variant.</li>
     *   <li>Mean coverage for each sample.</li>
     * </ul>
     *
     * @param noSamples The total number of samples used for statistical calculations.
     */
    private void statisticsOfVariants(int noSamples) {
        // Initialize per-sample count maps for various statistics.
        Map<String, Integer> sampleSubstitutions = new HashMap<>();
        Map<String, Integer> sampleInDels = new HashMap<>();
        Map<String, Integer> sampleCalls = new HashMap<>();
        Map<String, Integer> sampleFiltered = new HashMap<>();
        Map<String, Integer> sampleCoverage = new HashMap<>();

        // Iterate through all contigs in the storage.
        for (Contig contig : storage.getContigs()) {
            // Iterate through all variants in the current contig.
            for (Variant variant : contig.getAllVariants()) {
                // Retrieve the set of samples related to the current variant.
                Set<Tuple<String, String>> relatedSamples = variant.getRelatedSamples();

                // Calculate and set the variant frequency based on the number of related samples.
                variant.setAttribute(Constants.AttributesKeys.VARIANT_FREQUENCY,
                        IO.formatFrequency((float) relatedSamples.size() / noSamples));

                // Determine if the variant is an SNV or an insertion/deletion.
                int cSNV = variant.type.equals(Variant.Type.SNV) ? 1 : 0;
                int cInDel = (variant.type.equals(Variant.Type.INSERTION) || variant.type.equals(Variant.Type.DELETION)) ? 1 : 0;

                // Initialize counters for atomic calls, coverage sum, and entropy sum.
                int atomicCalls = 0;
                short coverageSum = 0;
                float entropySum = 0;

                // Iterate through all related samples for the current variant.
                for (Tuple<String, String> entry : relatedSamples) {
                    String sampleIdentifier = entry.a; // Sample identifier.
                    String[] variantCallStrings = entry.b.split(Constants.PIPE); // Variant call strings.

                    // Update per-sample counts for SNVs, insertions/deletions, and total calls.
                    sampleSubstitutions.merge(sampleIdentifier, cSNV, Integer::sum);
                    sampleInDels.merge(sampleIdentifier, cInDel, Integer::sum);
                    sampleCalls.merge(sampleIdentifier, variantCallStrings.length, Integer::sum);

                    // Update the count of filtered calls for the sample if the variant is filtered.
                    if (variant.isFiltered(sampleIdentifier)) {
                        sampleFiltered.merge(sampleIdentifier, 1, Integer::sum);
                    }

                    // Process each variant call string to calculate coverage and entropy.
                    for (String variantCallString : variantCallStrings) {
                        atomicCalls++;
                        String[] fields = variantCallString.split(Constants.SEMICOLON);
                        short coverage = Short.parseShort(fields[1]); // Coverage value.
                        coverageSum += coverage;
                        entropySum += Float.parseFloat(fields[2]); // Entropy value.
                        sampleCoverage.merge(sampleIdentifier, (int) coverage, Integer::sum);
                    }
                }

                // Calculate and set the mean coverage and entropy for the variant.
                variant.setAttribute(Constants.AttributesKeys.MEAN_COVERAGE, IO.formatNumber((float) coverageSum / atomicCalls));
                variant.setAttribute(Constants.AttributesKeys.MEAN_ENTROPY, IO.formatNumber(entropySum / atomicCalls));
            }
        }

        // Update sample attributes with the calculated statistics
        for (Sample sample : storage.getSamples()) {
            sample.setAttribute(Constants.AttributesKeys.NUMBER_OF_SNVS,
                    String.valueOf(sampleSubstitutions.getOrDefault(sample._id, 0)));
            sample.setAttribute(Constants.AttributesKeys.NUMBER_OF_INDELS,
                    String.valueOf(sampleInDels.getOrDefault(sample._id, 0)));
            sample.setAttribute(Constants.AttributesKeys.MEAN_COVERAGE,
                    IO.formatNumber((float) sampleCoverage.getOrDefault(sample._id, 0) /
                            sampleCalls.getOrDefault(sample._id, 1)));
            sample.setAttribute(Constants.AttributesKeys.FREQUENCY_FILTERED_CALLS,
                    IO.formatFrequency((float) sampleFiltered.getOrDefault(sample._id, 0) /
                            sampleCalls.getOrDefault(sample._id, 1)));
        }
    }

    /**
     * Updates statistical attributes in the {@link Storage} from iterating through all {@link Storage#samples}.
     * <p>
     * This method calculates and updates the following statistics for each sample:
     * <ul>
     *   <li>Frequency of reference alleles: The proportion of features where the sample has the reference allele.</li>
     *   <li>Frequency of disrupted coding features: The proportion of coding features where the sample's proteoform is disrupted.</li>
     * </ul>
     *
     * @param noFeatures       The total number of features in the storage, used to calculate reference allele frequency.
     * @param noCodingFeatures The total number of coding features in the storage, used to calculate disrupted feature frequency.
     */
    private void statisticsOfSamples(int noFeatures, int noCodingFeatures) {
        // Iterate through all samples in the storage.
        for (Sample sample : storage.getSamples()) {
            // Retrieve the set of alleles related to the current sample.
            Set<Tuple<String, String>> alleles = sample.getRelatedAlleles();

            // Calculate and set the frequency of reference alleles for the sample.
            sample.setAttribute(Constants.AttributesKeys.FREQUENCY_REFERENCE,
                    IO.formatFrequency(1 - (float) alleles.size() / noFeatures));

            int disrupted = 0; // Counter for disrupted coding features.

            // Iterate through all alleles related to the sample.
            for (var entry : alleles) {
                // Retrieve the feature associated with the allele.
                Feature feature = storage.getFeature(entry.a);

                // Check if the feature is coding and its proteoform is disrupted.
                if (feature.isCoding()) {
                    String proteoformIdentifier = feature.getAllele(entry.b).getRelatedProteoform();
                    if (!Objects.equals(proteoformIdentifier, Constants.SYNONYMOUS)) {
                        if (feature.getProteoform(proteoformIdentifier).isDisrupted()) {
                            disrupted++; // Increment the disrupted feature counter.
                        }
                    }
                }
            }

            // Calculate and set the frequency of disrupted coding features for the sample.
            sample.setAttribute(Constants.AttributesKeys.FREQUENCY_DISRUPTED,
                    IO.formatFrequency((float) disrupted / noCodingFeatures));
        }
    }

    /**
     * Updates statistical attributes in the {@link Storage} from iterating through all {@link Feature}.
     * <p>
     * This method calculates and updates various statistics for each feature, including:
     * <ul>
     *   <li>Allelic frequency and diversity.</li>
     *   <li>Proteoform frequency, diversity, and disrupted proteoform counts (if the feature is coding).</li>
     * </ul>
     * <p>
     * The {@code diversity} is based on the Simpson-Index and yields the probability that two randomly drawn samples of the population
     * do not have the same allele/proteoform, i.e., a value of {@code 0} indicates no diversity (all samples have the same
     * allele/proteoform), whereas a value of {@code 1} indicates maximum diversity (all samples have different alleles/proteoforms).
     *
     * @param noSamples The total number of samples used for statistical calculations.
     */
    private void statisticsOfFeatures(int noSamples) {
        // Calculate the denominator for diversity calculations.
        float complexityDenominator = noSamples * (noSamples - 1);

        // Iterate through all features in the storage.
        for (Feature feature : storage.getFeatures()) {
            // Determine if the feature is coding and typing is not skipped.
            boolean typedCoding = !storage.parameters.skipTyping() && feature.isCoding();
            int alternativeAlleleCount = 0; // Count of alternative alleles across all samples.
            float alleleDiversity = 0; // Diversity measure for alleles.

            int disruptedProteoformCount = 0; // Count of disrupted proteoforms.
            float proteoformDiversity = 0; // Diversity measure for proteoforms.
            Map<String, Integer> proteoformCounts = new HashMap<>(); // Map to store proteoform counts.

            // Iterate through all alleles of the feature.
            for (Allele allele : feature.getAlleles()) {
                int count = allele.getRelatedSamplesCount(); // Number of samples related to the allele.
                allele.setAttribute(Constants.AttributesKeys.ALLELIC_FREQUENCY, IO.formatFrequency(count / (float) noSamples));
                alternativeAlleleCount += count;
                alleleDiversity += (count * (count - 1)) / complexityDenominator;

                // If the feature is coding, process proteoform-related statistics.
                if (typedCoding) {
                    String proteoformIdentifier = allele.getRelatedProteoform();
                    if (!Objects.equals(proteoformIdentifier, Constants.SYNONYMOUS)) {
                        Collection<String> effects = feature.getProteoform(proteoformIdentifier)
                                .getAttributeSet(Constants.AttributesKeys.SO_EFFECTS);
                        if (effects.contains("start_lost") || effects.contains("stop_gained")) {
                            disruptedProteoformCount++;
                        }
                        proteoformCounts.merge(proteoformIdentifier, count, Integer::sum);
                    }
                }
            }

            // Calculate reference allele count and allele diversity.
            int referenceAlleleCount = noSamples - alternativeAlleleCount;
            alleleDiversity = 1 - (alleleDiversity + (referenceAlleleCount * (referenceAlleleCount - 1)) / complexityDenominator);

            // Update feature attributes with allele statistics.
            feature.setAttribute(Constants.AttributesKeys.NUMBER_OF_ALLELES, String.valueOf(feature.getAlleleCount()));
            feature.setAttribute(Constants.AttributesKeys.FREQUENCY_REFERENCE,
                    IO.formatFrequency(1 - (alternativeAlleleCount / (float) noSamples)));
            feature.setAttribute(Constants.AttributesKeys.DIVERSITY_ALLELE, IO.formatNumber(alleleDiversity));

            // If the feature is coding, calculate and update proteoform statistics.
            if (typedCoding) {
                int alternativeProteoformCount = 0;
                for (var entry : proteoformCounts.entrySet()) {
                    int count = entry.getValue();
                    feature.getProteoform(entry.getKey()).setAttribute(Constants.AttributesKeys.ALLELIC_FREQUENCY,
                            IO.formatFrequency(count / (float) noSamples));
                    alternativeProteoformCount += count;
                    proteoformDiversity += (count * (count - 1)) / complexityDenominator;
                }

                // Calculate reference proteoform count and proteoform diversity.
                int referenceProteoformCount = noSamples - alternativeProteoformCount;
                proteoformDiversity =
                        1 - (proteoformDiversity + (referenceProteoformCount * (referenceProteoformCount - 1)) / complexityDenominator);

                // Update feature attributes with proteoform statistics.
                feature.setAttribute(Constants.AttributesKeys.NUMBER_OF_PROTEOFORMS, String.valueOf(feature.getProteoformCount()));
                feature.setAttribute(Constants.AttributesKeys.FREQUENCY_DISRUPTED,
                        IO.formatFrequency(disruptedProteoformCount / (float) feature.getProteoformCount()));
                feature.setAttribute(Constants.AttributesKeys.DIVERSITY_PROTEOFORM, IO.formatNumber(proteoformDiversity));
            }
        }
    }

}
