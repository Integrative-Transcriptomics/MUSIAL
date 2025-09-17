package op;

import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import model.*;
import org.apache.commons.lang3.tuple.ImmutableTriple;
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
 * stored {@link Contig}s, calculate sequence types, and generate statistics. This is mainly to keep the {@link Storage} class cleaner and
 * more close to a POJO.
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
                storage.getSample(sampleIdentifier).addAttributesIfAbsent(entry.getValue());
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

            // Iterate through all samples to infer allele sequence types.
            for (Sample sample : storage.getSamples()) {
                // Update allele.
                Allele allele;
                String alleleIdentifier = sample.getRelatedAllele(feature._id);
                List<Variant.Stub> variants;

                // Note: In this special case, a returned reference allele may also indicate the presence of a new allele.
                if (alleleIdentifier.equals(Constants.REFERENCE)) {
                    // Retrieve variants for the sample within the feature's range.
                    variants = Bio.variantsAsStub(contig.getVariants(feature.start, feature.end, sample._id));
                    if (variants.isEmpty()) continue; // The sample has the reference allele; skip further processing.

                    // Determine the allele from the sample's variants.
                    allele = new Allele(variants);

                    // Check if the allele already exists for the feature, i.e., was created by another sample.
                    if (feature.hasAllele(allele._id)) {
                        allele = feature.getAllele(allele._id);
                    } else { // If the allele is new, annotate its effects and add relations.
                        // Calculate sequence length deviation and frameshift effects.
                        int lengthDelta = Integer.parseInt(allele.getAttribute(Constants.AttributesKeys.SEQUENCE_LENGTH_DEVIATION));
                        int netFrameshift = Math.abs(lengthDelta % 3);
                        Set<String> effects = contig.getVariantsEffects(variants);

                        // Add frameshift effects if applicable.
                        if (netFrameshift != 0) {
                            effects.add(lengthDelta > 0 ? "plus_%d_frameshift".formatted(netFrameshift)
                                    : "minus_%d_frameshift".formatted(netFrameshift));
                        }
                        allele.setAttribute(Constants.AttributesKeys.SO_EFFECTS, String.join(Constants.COMMA, effects));

                        // Add relations between variants and the allele.
                        for (Variant.Stub stub : variants) {
                            contig.getVariant(stub.position(), stub.alternative()).addRelation(feature._id, allele._id);
                        }

                        // Add the new allele to the feature.
                        feature.addAllele(allele);
                    }

                    // Add relations between the allele, sample, and feature.
                    allele.addRelation(sample._id);
                    sample.addRelation(feature._id, allele._id);
                }
            }

            // If the feature is coding and the contig has a sequence, update proteoform sequence types.
            if (feature.isCoding() && contig.hasSequence()) {

                // Iterate through all alleles of the feature to infer proteoform sequence types.
                for (Allele allele : feature.getAlleles()) {
                    Proteoform proteoform;

                    // Skip if a proteoform is already assigned.
                    String proteoformIdentifier = allele.getProteoform();
                    if (proteoformIdentifier != null) continue;

                    // Construct a position-sorted map of variants for the allele.
                    NavigableMap<Integer, String> variants = new TreeMap<>();
                    allele.getVariants().forEach(stub -> variants.put(stub.position(), stub.alternative()));

                    // Translate the reference and allele sequences.
                    String referenceSequence = Bio.translateSequence(contig.getSequence(feature.start, feature.end), feature.isReverse());
                    String proteoformSequence = Bio.translateSequence(Bio.integrateVariants(contig, feature, variants, true),
                            feature.isReverse());
                    assert referenceSequence != null && proteoformSequence != null;

                    // Check if the proteoform sequence is synonymous with the reference sequence.
                    if (referenceSequence.equals(proteoformSequence)) {
                        proteoformIdentifier = Constants.SYNONYMOUS;
                    } else {
                        // Align sequences and extract amino acid variants.
                        Tuple<String, String> alignment = Bio.globalProteinSequenceAlignment(
                                referenceSequence, proteoformSequence, Math.max(referenceSequence.length(), proteoformSequence.length()),
                                6, true, false, Math.abs(referenceSequence.length() - proteoformSequence.length())
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
                            if (alleleEffects.contains("frameshift")) effects.add("frameshift_sequence_variation");
                            if (aaVariants.stream().anyMatch(aav -> Bio.isInsertion(aav.alternative())))
                                effects.add("amino_acid_insertion");
                            if (aaVariants.stream().anyMatch(aav -> Bio.isDeletion(aav.alternative())))
                                effects.add("amino_acid_deletion");
                            if (proteoform.hasVariantAt(1) && proteoform.getVariant(1).charAt(0) != referenceSequence.charAt(0))
                                effects.add("start_lost");
                            aaVariants.stream()
                                    .filter(aav -> aav.alternative().contains(Constants.TERMINAL_AA))
                                    .findFirst()
                                    .ifPresent(stopCodonVariant -> {
                                        int stopCodonPosition = stopCodonVariant.position() +
                                                stopCodonVariant.alternative().indexOf(Constants.TERMINAL_AA);
                                        effects.add(stopCodonPosition <= referenceSequence.length() ? "stop_gained" :
                                                "redundant_inserted_stop_gained");
                                    });
                            proteoform.setAttribute(Constants.AttributesKeys.SO_EFFECTS, String.join(Constants.COMMA, effects));
                            feature.addProteoform(proteoform);
                        }

                        // Add the new proteoform to the feature and establish relations.
                        proteoformIdentifier = proteoform._id;
                        proteoform.addRelation(allele._id);
                    }

                    // Set the proteoform identifier for the allele.
                    allele.setProteoform(proteoformIdentifier);
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
        long noFeatures = storage.getFeatures().size();
        long noCodingFeatures = storage.getFeatures().stream().filter(Feature::isCoding).count();
        long noSamples = storage.getSamples().size();

        // Initialize maps to store per-sample substitution and InDel counts.
        Map<String, Integer> perSampleSubstitutions = new HashMap<>();
        Map<String, Integer> perSampleInDels = new HashMap<>();

        // Iterate through all samples to calculate sample-specific statistics.
        for (Sample sample : storage.getSamples()) {
            int totalCalls = 0, filteredCalls = 0;
            List<Short> coverages = new ArrayList<>();
            List<Float> entropy = new ArrayList<>();
            perSampleSubstitutions.put(sample._id, 0);
            perSampleInDels.put(sample._id, 0);

            // Process variant calls for the sample.
            for (ImmutableTriple<String, Integer, VariantCall> item : sample.getVariantCalls(false)) {
                VariantCall variantCall = item.right;
                totalCalls++;
                coverages.add(variantCall.totalDepth());
                if (!variantCall.isFiltered()) filteredCalls++;
                else entropy.add(variantCall.callEntropy());
            }

            // Update sample attributes with calculated statistics.
            sample.setAttribute(Constants.AttributesKeys.NUMBER_OF_CALLS, String.valueOf(totalCalls));
            sample.setAttribute(Constants.AttributesKeys.NUMBER_OF_FILTERED_CALLS, String.valueOf(filteredCalls));
            sample.setAttribute(Constants.AttributesKeys.MEAN_COVERAGE,
                    IO.formatNumber((short) coverages.stream().mapToInt(Short::intValue).average().orElse(0)));
            sample.setAttribute(Constants.AttributesKeys.MEAN_ENTROPY,
                    IO.formatNumber((float) entropy.stream().mapToDouble(Float::doubleValue).average().orElse(0)));
            sample.setAttribute(Constants.AttributesKeys.FREQUENCY_REFERENCE,
                    IO.formatFrequency(1 - (sample.getRelatedAllelesCount() / (float) noFeatures)));

            // Calculate disrupted coding feature frequency if typing is not skipped.
            if (!storage.parameters.skipTyping()) {
                int disrupted = 0;
                for (Tuple<String, String> item : sample.getRelatedAlleles()) {
                    Feature feature = storage.getFeature(item.a);
                    if (feature.isCoding()) {
                        String proteoformIdentifier = feature.getAllele(item.b).getProteoform();
                        if (!Constants.SYNONYMOUS.equals(proteoformIdentifier)) {
                            var effects = feature.getProteoform(proteoformIdentifier)
                                    .getAttributeSet(Constants.AttributesKeys.SO_EFFECTS);
                            if (effects.contains("start_lost") || effects.contains("stop_gained")) disrupted++;
                        }
                    }
                }
                sample.setAttribute(Constants.AttributesKeys.FREQUENCY_DISRUPTED,
                        IO.formatFrequency(disrupted / (float) noCodingFeatures));
            }
        }

        // Iterate through all contigs to calculate variant frequencies and sample-specific variant counts.
        for (Contig contig : storage.getContigs()) {
            for (Variant variant : contig.getVariants()) {
                int sampleCount = variant.getRelatedSamples().size();
                variant.setAttribute(Constants.AttributesKeys.VARIANT_FREQUENCY,
                        IO.formatFrequency(sampleCount / (float) noSamples));
                for (String sampleIdentifier : variant.getRelatedSamples()) {
                    Map<String, Integer> targetMap = switch (variant.type) {
                        case SNV -> perSampleSubstitutions;
                        case INSERTION, DELETION -> perSampleInDels;
                    };
                    targetMap.put(sampleIdentifier, targetMap.get(sampleIdentifier) + 1);
                }
            }
        }

        // Update sample attributes with substitution and InDel counts.
        perSampleSubstitutions.forEach((sampleIdentifier, count) -> storage.getSample(sampleIdentifier)
                .setAttribute(Constants.AttributesKeys.NUMBER_OF_SNVS, String.valueOf(count)));
        perSampleInDels.forEach((sampleIdentifier, count) -> storage.getSample(sampleIdentifier)
                .setAttribute(Constants.AttributesKeys.NUMBER_OF_INDELS, String.valueOf(count)));

        // Initialize a map to store proteoform occurrences for features.
        Map<String, Integer> perProteoformOccurrence = new HashMap<>();

        // Iterate through all features to calculate feature-specific statistics.
        for (Feature feature : storage.getFeatures()) {
            float nonReferenceOccurrence = 0;
            int disrupted = 0;
            perProteoformOccurrence.clear();

            // Process alleles for the feature.
            for (Allele allele : feature.getAlleles()) {
                int alleleOccurrence = allele.getRelatedSamplesCount();
                allele.setAttribute(Constants.AttributesKeys.ALLELIC_FREQUENCY,
                        IO.formatFrequency(alleleOccurrence / (float) noSamples));
                nonReferenceOccurrence += alleleOccurrence;

                // Calculate proteoform frequencies and disrupted proteoform counts if typing is not skipped.
                if (!storage.parameters.skipTyping() && feature.isCoding()) {
                    String proteoformUid = allele.getProteoform();
                    if (!Objects.equals(proteoformUid, Constants.SYNONYMOUS)) {
                        Collection<String> effects = feature.getProteoform(proteoformUid)
                                .getAttributeSet(Constants.AttributesKeys.SO_EFFECTS);
                        if (effects.contains("start_lost") || effects.contains("stop_gained")) disrupted++;
                        perProteoformOccurrence.merge(proteoformUid, alleleOccurrence, Integer::sum);
                    }
                }
            }

            // Update feature attributes with calculated statistics.
            feature.setAttribute(Constants.AttributesKeys.FREQUENCY_REFERENCE,
                    IO.formatFrequency(1 - (nonReferenceOccurrence / noSamples)));
            feature.setAttribute(Constants.AttributesKeys.NUMBER_OF_ALLELES, String.valueOf(feature.getAlleleCount()));

            if (!storage.parameters.skipTyping() && feature.isCoding()) {
                int proteoformCount = feature.getProteoformCount();
                float disruptedFrequency = proteoformCount == 0 ? 0 : disrupted / (float) proteoformCount;
                feature.setAttribute(Constants.AttributesKeys.FREQUENCY_DISRUPTED,
                        IO.formatFrequency(disruptedFrequency));
                feature.setAttribute(Constants.AttributesKeys.NUMBER_OF_PROTEOFORMS, String.valueOf(proteoformCount));
                perProteoformOccurrence.forEach((proteoformUid, count) ->
                        feature.getProteoform(proteoformUid).setAttribute(Constants.AttributesKeys.ALLELIC_FREQUENCY,
                                IO.formatFrequency(count / (float) noSamples)));
            }
        }
    }

}
