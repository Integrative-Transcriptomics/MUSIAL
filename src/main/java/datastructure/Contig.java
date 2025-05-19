package datastructure;

import htsjdk.samtools.util.Tuple;
import utility.Constants;
import utility.IO;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents a reference sequence segment.
 * <p>
 * This class models a segment of a reference sequence, which can represent a complete genome,
 * a plasmid, a single contig, or a scaffold. It extends the {@link Attributable} class to
 * inherit functionality for managing attributes associated with the contig.
 * </p>
 * <p>
 * Each instance of this class is uniquely identified by its {@code name} and contains
 * information about its nucleotide sequence, variants, and other relevant properties.
 * </p>
 */
public class Contig extends Attributable {

    /**
     * The name or internal identifier of this contig.
     * <p>
     * This field uniquely identifies the contig within the context of the application.
     * It is a final field, meaning its value is immutable once assigned during the
     * construction of the {@link Contig} instance.
     * </p>
     */
    public final String name;

    /**
     * The sequence of this contig.
     * <p>
     * This field stores the nucleotide sequence of the contig. The sequence is expected to be
     * stored as a GZIP-compressed string to optimize storage. It may be empty or null if no
     * sequence is available for the contig.
     * </p>
     * <p>
     * <b>Note:</b> The sequence is not validated against the variants stored in the {@code variants} map.
     * </p>
     */
    protected final String sequence;

    /**
     * Hierarchical map structure to store variants located on this contig.
     * <p>
     * This map organizes variants in a hierarchical structure:
     * <ul>
     *     <li>The first level (key: {@link Integer}) represents the position of the variant on the contig.</li>
     *     <li>The second level (key: {@link String}) represents the alternative variant content (e.g., alternate alleles).</li>
     *     <li>The third level (value: {@link VariantInformation}) contains additional information about the variant,
     *         such as associations with {@link SequenceType}s or {@link Sample}s.</li>
     * </ul>
     * </p>
     * <p>
     * This structure allows storage and retrieval of variant data, enabling queries by position,
     * alternative content, and associated metadata.
     * </p>
     */
    protected final TreeMap<Integer, Map<String, VariantInformation>> variants;

    /**
     * Cache to store the sequence of a contig for a specific start and end position.
     * <p>
     * This field is a transient {@link HashMap} used to cache subsequences of the contig's sequence.
     * The keys in the map are {@link Tuple} objects representing the start and end positions of the subsequence,
     * and the values are the corresponding subsequences as {@link String}.
     * </p>
     * <p>
     * The cache is transient because it is not intended to be serialized, as it is dynamically populated
     * during runtime to optimize performance by avoiding redundant sequence decompression or retrieval.
     * </p>
     */
    protected transient HashMap<Tuple<Integer, Integer>, String> sequenceCache;

    /**
     * Constructs a new {@link Contig} instance with the specified name and sequence.
     * <p>
     * This constructor initializes a contig with its name and nucleotide sequence. It also
     * initializes the {@link #variants} map to store variant information and the {@link #sequenceCache}
     * map to cache subsequences for optimized retrieval. The sequence is expected to be stored
     * as a GZIP-compressed string to reduce storage requirements.
     * </p>
     *
     * @param name     The name or identifier of the contig.
     * @param sequence The nucleotide sequence of the contig, stored as a GZIP-compressed string.
     */
    protected Contig(String name, String sequence) {
        super();
        this.name = name;
        this.sequence = sequence;
        this.variants = new TreeMap<>(Integer::compare);
        this.sequenceCache = new HashMap<>();
    }

    /**
     * Checks if this contig has an associated nucleotide sequence.
     * <p>
     * This method determines whether the contig has a stored sequence by checking
     * if the {@code sequence} field is not empty. A non-empty sequence indicates
     * that the contig has an associated nucleotide sequence.
     * </p>
     *
     * @return {@code true} if the contig has a sequence (i.e., the sequence length is not zero),
     * {@code false} otherwise.
     */
    public boolean hasSequence() {
        return !sequence.isEmpty();
    }

    /**
     * Retrieves the nucleotide sequence of this contig or an empty string if no sequence is stored.
     * <p>
     * This method decompresses the GZIP-compressed sequence stored in the {@code sequence} field
     * and returns it as a string. If no sequence is stored, it returns an empty string.
     * </p>
     *
     * @return The decompressed nucleotide sequence of this contig, or an empty string if no sequence is stored.
     * @throws IOException If an error occurs during the decompression of the sequence.
     */
    public String getSequence() throws IOException {
        if (hasSequence())
            return IO.gzipDecompress(this.sequence);
        else
            return Constants.EMPTY;
    }

    /**
     * Retrieves a subsequence of this contig, caching the result to optimize performance.
     * <p>
     * This method extracts a subsequence from the nucleotide sequence of the contig based on the
     * specified start and end positions. The subsequence is cached to avoid redundant decompression
     * and substring operations for the same range. If the subsequence is already cached, it is
     * retrieved directly from the cache. Otherwise, it is computed, stored in the cache, and returned.
     * </p>
     *
     * <p>
     * The start and end positions are 1-based indices, meaning the first nucleotide in the sequence
     * is at position 1. If no sequence is stored for the contig, the method returns an empty string.
     * </p>
     *
     * @param start The 1-based indexed start position of the subsequence (inclusive).
     * @param end   The 1-based indexed end position of the subsequence (exclusive).
     * @return The subsequence of this contig, or an empty string if no sequence is stored.
     * @throws IOException If an error occurs during the decompression of the sequence.
     */
    public String getSubsequence(int start, int end) throws IOException {
        if (hasSequence()) {
            Tuple<Integer, Integer> cacheKey = new Tuple<>(start, end);
            if (sequenceCache.containsKey(cacheKey)) {
                return sequenceCache.get(cacheKey);
            } else {
                String subsequence = getSequence().substring(start - 1, end);
                sequenceCache.put(cacheKey, subsequence);
                return subsequence;
            }
        } else {
            return Constants.EMPTY;
        }
    }

    /**
     * Calculates the total number of variants located on this contig.
     * <p>
     * This method iterates through the hierarchical map of variants stored in the {@code variants} field.
     * It computes the total count by summing up the sizes of all inner maps, where each inner map represents
     * the alternative sequences for a specific position on the contig.
     * </p>
     *
     * @return The total number of variants located on this contig.
     */
    public int getVariantsCount() {
        return this.variants.values().stream().mapToInt(Map::size).sum();
    }

    /**
     * Retrieves the {@link VariantInformation} associated with a specific variant on this contig.
     * <p>
     * This method accesses the hierarchical map of variants to retrieve the {@link VariantInformation}
     * for a variant located at the specified position with the given alternative bases. The returned
     * {@link VariantInformation} contains details about the variant, including its occurrences in
     * samples and features, as well as any associated attributes.
     * </p>
     *
     * @param position         The 1-based position of the variant on the contig.
     * @param alternativeBases The alternative base sequence of the variant.
     * @return The {@link VariantInformation} associated with the specified variant, or {@code null}
     * if no such variant exists at the given position with the specified alternative bases.
     */
    public VariantInformation getVariantInformation(int position, String alternativeBases) {
        return this.variants.get(position).get(alternativeBases);
    }

    /**
     * Extracts and aggregates the SnpEff effects from a list of variants.
     * <p>
     * This method processes a list of {@link Tuple} objects, where each tuple contains:
     * <ul>
     *   <li>The position of the variant (field {@code a} of the tuple).</li>
     *   <li>The alternate allele of the variant (field {@code b} of the tuple).</li>
     * </ul>
     * </p>
     *
     * <p>
     * For each variant, the method retrieves the associated {@link VariantInformation} using
     * the position and alternate allele. It then checks if the {@code Constants.EFFECTS} attribute
     * is present. If the attribute is found, its value (a comma-separated string of effects) is split
     * into individual effects, which are trimmed and aggregated into a {@link Set} to ensure uniqueness.
     * </p>
     *
     * @param variants A list of {@link Tuple} objects representing the variants. Each tuple contains:
     *                 <ul>
     *                   <li>The position of the variant.</li>
     *                   <li>The alternate allele sequence.</li>
     *                 </ul>
     * @return A {@link Set} of unique SnpEff effects extracted from the variants.
     */
    public Set<String> getVariantsEffects(List<Tuple<Integer, String>> variants) {
        return variants.stream()
                .map(variant -> getVariantInformation(variant.a, variant.b))
                .flatMap(variantInfo -> variantInfo.getAttributeAsCollection(Constants.snpEffAttributeKeyPrefix + Constants.snpEffKeys.get(1)).stream())
                .collect(Collectors.toSet());
    }

    /**
     * Retrieves all variants located on this contig.
     * <p>
     * This method iterates through the hierarchical {@code variants} map, which organizes
     * variants by their positions and alternative sequences. For each variant, it creates
     * a {@link Tuple} containing:
     * <ul>
     *   <li>The position of the variant on the contig.</li>
     *   <li>The alternative base sequence of the variant.</li>
     * </ul>
     * </p>
     *
     * @return An {@link ArrayList} of {@link Tuple} objects, where each tuple contains
     * the position and alternative base sequence of a variant.
     */
    public ArrayList<Tuple<Integer, String>> getVariants() {
        ArrayList<Tuple<Integer, String>> variants = new ArrayList<>();
        this.variants.forEach((position, innerMap) ->
                innerMap.forEach((alternativeBases, variantInformation) ->
                        variants.add(new Tuple<>(position, alternativeBases))
                )
        );
        return variants;
    }

    /**
     * Retrieves variants located within a specified range on this contig.
     * <p>
     * This method filters the {@code variants} map to identify variants that fall within the
     * specified start and end positions. For each variant in the range, it creates a {@link Tuple}
     * containing:
     * <ul>
     *   <li>The position of the variant on the contig.</li>
     *   <li>The alternative base sequence of the variant.</li>
     * </ul>
     * </p>
     *
     * @param start The 1-based indexed inclusive start position of the range.
     * @param end   The 1-based indexed inclusive end position of the range.
     * @return An {@link ArrayList} of {@link Tuple} objects, where each tuple contains
     * the position and alternative base sequence of a variant within the specified range.
     */
    public ArrayList<Tuple<Integer, String>> getVariantsByLocation(int start, int end) {
        ArrayList<Tuple<Integer, String>> variants = new ArrayList<>();
        this.variants.subMap(start, end + 1).forEach((position, innerMap) ->
                innerMap.forEach((alternativeBases, variantInformation) ->
                        variants.add(new Tuple<>(position, alternativeBases))
                )
        );
        return variants;
    }

    /**
     * Retrieves variants associated with specific alleles of a feature.
     * <p>
     * This method filters the {@code variants} map to identify variants that are associated with the specified
     * feature and alleles. It first retrieves a sub-map of variants within the location range of the feature
     * and then filters the inner maps to find variants matching the alleles. Each matching variant is represented
     * as a {@link Tuple} containing:
     * <ul>
     *   <li>The position of the variant on the contig.</li>
     *   <li>The alternative base sequence of the variant.</li>
     * </ul>
     * </p>
     *
     * @param feature    The feature to filter variants by.
     * @param alleleUids A set of allele unique identifiers to filter variants by.
     * @return An {@link ArrayList} of {@link Tuple} objects, where each tuple contains
     * the position and alternative base sequence of a variant associated with the specified feature
     * and alleles within the specified location range.
     */
    public ArrayList<Tuple<Integer, String>> getVariantsByAlleles(Feature feature, Set<String> alleleUids) {
        ArrayList<Tuple<Integer, String>> variants = new ArrayList<>();
        this.variants.subMap(feature.start, feature.end + 1).forEach((position, innerMap) ->
                innerMap.forEach((alternativeBases, variantInformation) -> {
                            if (alleleUids.stream().anyMatch(uid -> variantInformation.hasOccurrence(feature.name, uid)))
                                variants.add(new Tuple<>(position, alternativeBases));
                        }
                )
        );
        return variants;
    }

    /**
     * Retrieves variants associated with a specific sample.
     * <p>
     * This method filters the {@code variants} map to identify variants that are associated
     * with the specified sample. For each position in the map, it checks the inner map of
     * alternative sequences and retrieves the first variant that matches the sample. Each
     * matching variant is represented as a {@link Tuple} containing:
     * <ul>
     *   <li>The position of the variant on the contig.</li>
     *   <li>The alternative base sequence of the variant.</li>
     * </ul>
     * </p>
     *
     * @param sampleName The name of the sample to filter variants by.
     * @return An {@link ArrayList} of {@link Tuple} objects, where each tuple contains
     * the position and alternative base sequence of a variant associated with the specified sample.
     */
    public ArrayList<Tuple<Integer, String>> getVariantsBySample(String sampleName) {
        ArrayList<Tuple<Integer, String>> sampleVariants = new ArrayList<>();
        this.variants.forEach((position, innerMap) ->
                innerMap.entrySet().stream()
                        .filter(e -> e.getValue().hasOccurrence(Constants.$Attributable_samplesOccurrence, sampleName))
                        .findFirst()
                        .ifPresent(e -> sampleVariants.add(new Tuple<>(position, e.getKey())))
        );
        return sampleVariants;
    }

    /**
     * Retrieves variants associated with a specific sample within a specified location range.
     * <p>
     * This method filters the {@code variants} map to identify variants that are associated
     * with the specified sample and fall within the given start and end positions. It first
     * retrieves a sub-map of variants within the location range and then filters the inner
     * maps to find variants matching the sample. Each matching variant is represented as a
     * {@link Tuple} containing:
     * <ul>
     *   <li>The position of the variant on the contig.</li>
     *   <li>The alternative base sequence of the variant.</li>
     * </ul>
     * </p>
     *
     * @param sampleName The name of the sample to filter variants by.
     * @param start      The 1-based indexed inclusive start position of the location range.
     * @param end        The 1-based indexed inclusive end position of the location range.
     * @return An {@link ArrayList} of {@link Tuple} objects, where each tuple contains
     * the position and alternative base sequence of a variant associated with the specified sample
     * and within the specified location range.
     */
    public ArrayList<Tuple<Integer, String>> getVariantsBySampleAndLocation(String sampleName, int start, int end) {
        ArrayList<Tuple<Integer, String>> sampleVariants = new ArrayList<>();
        this.variants.subMap(start, end + 1).forEach((position, innerMap) ->
                innerMap.entrySet().stream()
                        .filter(e -> e.getValue().hasOccurrence(Constants.$Attributable_samplesOccurrence, sampleName))
                        .findFirst()
                        .ifPresent(e -> sampleVariants.add(new Tuple<>(position, e.getKey())))
        );
        return sampleVariants;
    }

}
