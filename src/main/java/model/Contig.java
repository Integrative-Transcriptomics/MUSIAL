package model;

import htsjdk.samtools.util.Tuple;
import uk.co.omegaprime.btreemap.BTreeMap;
import util.Constants;
import util.IO;
import util.Logging;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents a reference genomic location.
 * <p>
 * Models a segment of a genomic sequence, i.e., a complete genome, plasmid, single contig or scaffold. It extends the {@link Attributes}
 * class to inherit functionality for managing attributes associated with the contig. In addition, an inner map is used to store
 * {@link Variant}s associated with the contig.
 * <p>
 * Contigs are stored in the {@link Storage#contigs} property in the model.
 */
public class Contig extends Attributes {

    /**
     * Unique identifier of this contig.
     * <p>
     * This field serves as the unique identifier for the contig and is used to reference it in the model. This should be at best a database
     * identifier, such as a NCBI accession number.
     */
    public final String _id;

    /**
     * Length of the contig's sequence.
     */
    private final int length;

    /**
     * Sequence of this contig.
     * <p>
     * This field stores the nucleotide sequence of the contig. The sequence is expected to be stored as a GZIP-compressed string to
     * optimize storage. It may be empty or null if no sequence is available for the contig.
     * <p>
     * <b>Note:</b> The sequence is not validated against the variants stored in the {@code variants} map.
     */
    private final String sequence;

    /**
     * Map to store variants associated with this contig.
     * <p>
     * The map organizes {@link Variant} instances by their positions to allow more efficient queries.
     */
    private final NavigableMap<Integer, Map<String, Variant>> variants;

    /**
     * Cache to store the (sub-)sequence of this contig given a start and end position.
     * <p>
     * This field is a transient {@link HashMap} used to cache subsequences of the contig's sequence. The keys in the map are {@link Tuple}
     * objects representing the start and end positions of the subsequence, and the values are the corresponding subsequences as
     * {@link String}.
     * <p>
     * This is not intended to be serialized, as it is dynamically populated during runtime to optimize performance by avoiding redundant
     * sequence decompression or retrieval.
     */
    protected transient Map<Tuple<Integer, Integer>, String> sequenceCache;

    /**
     * Constructs a new {@link Contig} instance with the specified sequence (optionally empty).
     * <p>
     * This constructor initializes a contig with its identifier and nucleotide sequence. The {@link #variants} map to store variant
     * information and the {@link #sequenceCache} map to cache subsequences for optimized retrieval are initialized as empty instances. If a
     * non-empty sequence is provided, it is compressed using GZIP to reduce storage requirements.
     *
     * @param identifier The unique identifier for this contig, which is used to reference it in the context of a {@link Storage} instance.
     * @param sequence   The nucleotide sequence of the contig, stored as a GZIP-compressed string.
     * @throws IOException If an error occurs during the compression of the sequence.
     */
    Contig(String identifier, String sequence) throws IOException {
        super(); // Call the constructor of the parent class (Attributes).
        this._id = identifier; // Assign the unique identifier to the contig.

        String compressedSequence; // Variable to store the compressed sequence.

        // Check if the sequence is non-null and not empty.
        if (Objects.nonNull(sequence) && !sequence.isEmpty()) {
            compressedSequence = IO.gzipCompress(sequence); // Compress the sequence using GZIP.
            this.length = sequence.length(); // Calculate the length of the sequence.
        } else {
            compressedSequence = Constants.EMPTY; // Assign an empty string if the sequence is null or empty.
            this.length = 0; // Set the length to 0 for an empty sequence.
        }

        this.sequence = compressedSequence; // Store the compressed sequence.
        this.variants = BTreeMap.create(); // Initialize the map to store variants.
        this.sequenceCache = new HashMap<>(); // Initialize the cache for subsequences.
    }

    /**
     * Checks if this contig has an associated nucleotide sequence.
     * <p>
     * This method determines whether the contig has a stored sequence by checking if the {@code sequence} field is not empty. A non-empty
     * sequence indicates that the contig has an associated nucleotide sequence.
     *
     * @return {@code true} if the contig has a sequence or {@code false} otherwise.
     */
    public boolean hasSequence() {
        return !sequence.isEmpty();
    }

    /**
     * Retrieves the full nucleotide sequence of this contig or an empty string if no sequence is stored.
     * <p>
     * This method decompresses the GZIP-compressed sequence stored in the {@code sequence} field and returns it as a string. If no sequence
     * is stored, it returns an empty string.
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
     * This method extracts a subsequence from the nucleotide sequence of the contig based on the specified start and end positions. The
     * subsequence is cached to avoid redundant decompression and substring operations for the same range. If the subsequence is already
     * cached, it is retrieved directly from the cache. Otherwise, it is computed, stored in the cache, and returned.
     * <p>
     * The start and end positions are 1-based indices, meaning the first nucleotide in the sequence is at position 1. If no sequence is
     * stored for the contig, the method returns an empty string.
     *
     * @param start The 1-based indexed start position of the subsequence (inclusive).
     * @param end   The 1-based indexed end position of the subsequence (exclusive).
     * @return The subsequence of this contig, or an empty string if no sequence is stored.
     * @throws IOException If an error occurs during the decompression of the sequence.
     */
    public String getSequence(int start, int end) throws IOException {
        if (Objects.isNull(sequenceCache)) sequenceCache = new HashMap<>();
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
     * Retrieves the length of the contig's sequence.
     * <p>
     * This method returns the length of the nucleotide sequence associated with this contig. The length is determined during the
     * initialization of the contig and reflects the number of bases in the sequence. If the contig does not have an associated sequence,
     * the length is 0.
     *
     * @return The length of the contig's sequence as an integer.
     */
    public int getSequenceLength() {
        return this.length;
    }

    /**
     * Retrieves a {@link Variant} from the contig based on the specified position and alternative base sequence.
     * <p>
     * This method checks if the {@code variants} map contains the specified position as a key. If the position exists, it retrieves the
     * list of variants at that position, filters the list to find the first variant that matches the specified alternative base sequence,
     * and returns it. If no matching variant is found, the method returns {@code null}.
     *
     * @param position    The 1-based position of the variant on the contig.
     * @param alternative The alternative base sequence of the variant.
     * @return The {@link Variant} that matches the specified position and alternative base sequence, or {@code null} if no match is found.
     */
    public Variant getVariant(int position, String alternative) {
        if (!this.variants.containsKey(position)) return null;
        return this.variants.get(position).get(alternative);
    }

    /**
     * Retrieves all variants associated with this contig.
     * <p>
     * This method flattens the {@code variants} map, which organizes variants by their positions, into a single list of {@link Variant}
     * objects. The returned list is unmodifiable.
     *
     * @return A {@link List} containing all {@link Variant} objects associated with this contig.
     */
    public List<Variant> getVariants() {
        return this.variants.values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .toList();
    }

    /**
     * Retrieves all novel variants associated with this contig.
     * <p>
     * This method flattens the {@code variants} map, which organizes variants by their positions, into a single list of {@link Variant}
     * objects. It then filters the list to include only those variants that are marked as novel (i.e., have the {@code novel} property set
     * to {@code true}). The returned list is unmodifiable.
     *
     * @return A {@link List} containing all novel {@link Variant} objects associated with this contig.
     */
    public List<Variant> getNovelVariants() {
        return this.variants.values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .filter(v -> v.novel)
                .toList();
    }

    /**
     * Retrieves variants within the specified range of positions.
     * <p>
     * This method retrieves variants from the {@code variants} map that fall within the specified start and end positions (inclusive of
     * start, exclusive of end). The resulting variants are flattened into a single list. The returned list is unmodifiable.
     *
     * @param start The 1-based start position of the range (inclusive).
     * @param end   The 1-based end position of the range (exclusive).
     * @return A {@link List} of {@link Variant} objects within the specified range.
     */
    public List<Variant> getVariants(int start, int end) {
        return this.variants.subMap(start, end + 1).values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .toList();
    }

    /**
     * Retrieves variants based on the provided set of sample names.
     * <p>
     * This method filters the {@code variants} map to include only those variants that are associated with at least one of the specified
     * sample names. The returned list is unmodifiable.
     *
     * @param relations A variable-length array of identifiers to filter the variants.
     * @return A {@link List} of {@link Variant} objects that match the specified sample names.
     */
    public List<Variant> getVariants(String... relations) {
        return this.variants.values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .filter(variant -> Arrays.stream(relations).anyMatch(variant::hasRelation))
                .toList();
    }

    /**
     * Retrieves variants within the specified range of positions and filters them based on the provided set of sample names.
     * <p>
     * This method retrieves variants from the {@code variants} map that fall within the specified start and end positions (inclusive of
     * start, exclusive of end), and filters them to include only those variants that are associated with at least one of the specified
     * sample names. The returned list is unmodifiable.
     *
     * @param start     The 1-based start position of the range (inclusive).
     * @param end       The 1-based end position of the range (exclusive).
     * @param relations A variable-length array of identifiers to filter the variants.
     * @return A {@link List} of {@link Variant} objects within the specified range and matching the specified sample names.
     */
    public List<Variant> getVariants(int start, int end, String... relations) {
        return this.variants.subMap(start, end + 1).values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .filter(variant -> Arrays.stream(relations).anyMatch(variant::hasRelation))
                .toList();
    }

    /**
     * Calculates the total number of variants associated with this contig.
     * <p>
     * This method iterates through the {@code variants} map and sums up the sizes of all the lists of variants. The result represents the
     * total count of {@link Variant} objects stored in this contig.
     *
     * @return The total number of {@link Variant} objects associated with this contig.
     */
    public int getVariantsCount() {
        return this.variants.values().stream().mapToInt(Map::size).sum();
    }

    /**
     * Calculates the total number of novel variants associated with this contig.
     * <p>
     * This method flattens the {@code variants} map into a stream of {@link Variant} objects. It then filters the stream to include only
     * those variants that are marked as novel (i.e., have the {@code novel} property set to {@code true}). The method counts the filtered
     * variants and returns the total count as an integer.
     *
     * @return The total number of novel {@link Variant} objects associated with this contig.
     */
    public int getNovelVariantsCount() {
        return (int) this.variants.values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .filter(v -> v.novel)
                .count();
    }

    /**
     * Retrieves the set of variant effects for a given list of variant stubs.
     * <p>
     * This method processes a list of {@link Variant.Stub} objects, retrieves the corresponding {@link Variant} objects from the contig,
     * and extracts their associated effects. The effects are determined by accessing the attribute set of each variant using a predefined
     * key.
     * </p>
     *
     * @param variants A {@link List} of {@link Variant.Stub} objects representing the variants to process.
     * @return A {@link Set} of {@link String} containing the effects associated with the given variants.
     */
    public Set<String> getVariantsEffects(List<Variant.Stub> variants) {
        return variants.stream()
                .map(v -> getVariant(v.position(), v.alternative())) // Retrieve the Variant object for each stub.
                .filter(Objects::nonNull) // Ensure only non-null Variant objects are processed.
                .flatMap(V -> V.getAttributeSet(Constants.SNP_EFF_PREFIX + Constants.SNP_EFF_KEYS.get(1)).stream())
                .collect(Collectors.toSet()); // Collect the effects into a set to ensure uniqueness.
    }

    /**
     * Adds a variant to the contig's variant map.
     * <p>
     * This method ensures that the {@code variants} map contains an entry for the specified position. If no entry exists, a new
     * {@link HashMap} is created for that position. The method then adds the provided {@link Variant} to the map at the specified position,
     * using the variant's alternative base sequence as the key. If a variant with the same alternative base sequence already exists at the
     * position, it checks whether the reference base of the existing variant matches the new variant. If the reference bases differ, a
     * warning is logged.
     *
     * @param variant The {@link Variant} object to be added to the contig's variant map.
     */
    void addVariant(Variant variant) {
        this.variants.computeIfAbsent(variant.position, p -> new HashMap<>());
        Variant previous = this.variants.get(variant.position).putIfAbsent(variant.alternative, variant);
        if (Objects.nonNull(previous) && !variant.reference.equals(previous.reference)) {
            Logging.logWarning((("Variant at position %d on contig %s with alternative '%s' already exists with a different reference " +
                    "base" +
                    "(%s and %s).").formatted(variant.position, this._id, variant.alternative, previous.reference, variant.reference)));
        }
    }

    /**
     * Returns the string representation of this contig.
     * <p>
     * This method overrides the {@code toString} method to return the unique identifier of the contig.
     *
     * @return The unique identifier of this contig as a {@link String}.
     */
    public String toString() {
        return "%s%s%s".formatted(this._id, Constants.COLON, IO.md5Hash(this.sequence));
    }

    /**
     * Computes the hash code for this contig.
     * <p>
     * This method overrides the {@code hashCode} method to compute the hash code based on the unique identifier of the contig.
     *
     * @return The hash code of this contig.
     */
    public int hashCode() {
        return this._id.hashCode();
    }

    /**
     * Compares this contig to another object for equality.
     * <p>
     * This method overrides the {@code equals} method to compare the unique identifier of this contig with another object. Two contigs are
     * considered equal if they are of the same class and have the same unique identifier.
     *
     * @param obj The object to compare with this contig.
     * @return {@code true} if the objects are equal; {@code false} otherwise.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        Contig that = (Contig) obj;
        return this._id.equals(that._id);
    }

}
