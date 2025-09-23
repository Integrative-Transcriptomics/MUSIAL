package model;

import htsjdk.samtools.util.Tuple;
import uk.co.omegaprime.btreemap.BTreeMap;
import util.Constants;
import util.IO;
import util.Logging;

import java.io.IOException;
import java.util.*;

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
     * Retrieves a variant associated with the specified position and alternative base sequence.
     * <p>
     * This method checks if the {@code variants} map contains an entry for the given position. If no entry exists, it returns {@code null}.
     * Otherwise, it retrieves the {@link Variant} object associated with the specified alternative base sequence at the given position.
     *
     * @param position    The 1-based position of the variant to retrieve.
     * @param alternative The alternative base sequence of the variant to retrieve.
     * @return The {@link Variant} object associated with the specified position and alternative base sequence, or {@code null} if no such
     * variant exists.
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
    public List<Variant> getAllVariants() {
        return this.variants.values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .toList();
    }

    /**
     * Retrieves all active variants associated with this contig.
     * <p>
     * This method flattens the {@code variants} map, which organizes variants by their positions, into a single list of {@link Variant}
     * objects. It then filters the list to include only those variants that are marked as active (i.e., have the {@code active} property
     * set to {@code true}). The returned list is unmodifiable.
     *
     * @return A {@link List} containing all active {@link Variant} objects associated with this contig.
     */
    public List<Variant> getActiveVariants() {
        return this.variants.values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .filter(v -> v.active)
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
    public List<Variant> getVariantsWithin(int start, int end) {
        return this.variants.subMap(start, end + 1).values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .toList();
    }

    /**
     * Retrieves variants at the specified positions.
     * <p>
     * This method retrieves variants from the {@code variants} map that are located at the specified positions. The resulting variants are
     * flattened into a single list. The returned list is unmodifiable.
     *
     * @param positions An array of 1-based positions to retrieve variants from.
     * @return A {@link List} of {@link Variant} objects located at the specified positions.
     */
    public List<Variant> getVariantsAt(int... positions) {
        return Arrays.stream(positions)
                .boxed()
                .filter(this.variants::containsKey)
                .flatMap(pos -> this.variants.get(pos).values().stream())
                .toList();
    }

    /**
     * Retrieves variants associated with the specified sample identifiers.
     * <p>
     * This method filters the variants stored in the contig to include only those that are associated with at least one of the specified
     * sample identifiers. The resulting list is unmodifiable.
     *
     * @param sampleIdentifiers An array of sample identifiers to filter the variants.
     * @return A {@link List} of {@link Variant} objects associated with the specified sample identifiers.
     */
    public List<Variant> getVariantsOfSamples(String... sampleIdentifiers) {
        return this.variants.values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .filter(variant -> variant.ofSamples(sampleIdentifiers))
                .toList();
    }

    /**
     * Retrieves variants associated with the specified sample identifiers within a given range of positions.
     * <p>
     * This method filters the variants stored in the contig to include only those that are associated with at least one of the specified
     * sample identifiers and fall within the specified start and end positions (inclusive of start, exclusive of end). The resulting list
     * is unmodifiable.
     *
     * @param start             The 1-based start position of the range (inclusive).
     * @param end               The 1-based end position of the range (exclusive).
     * @param sampleIdentifiers An array of sample identifiers to filter the variants.
     * @return A {@link List} of {@link Variant} objects associated with the specified sample identifiers within the given range.
     */
    public List<Variant> getVariantsOfSamplesWithin(int start, int end, String... sampleIdentifiers) {
        return this.variants.subMap(start, end + 1).values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .filter(variant -> variant.ofSamples(sampleIdentifiers))
                .toList();
    }

    /**
     * Retrieves variants associated with the specified allele identifiers.
     * <p>
     * This method filters the variants stored in the contig to include only those that are associated with at least one of the specified
     * allele identifiers. The resulting list is unmodifiable.
     *
     * @param alleleIdentifiers An array of allele identifiers to filter the variants.
     * @return A {@link List} of {@link Variant} objects associated with the specified allele identifiers.
     */
    public List<Variant> getVariantsOfAlleles(String... alleleIdentifiers) {
        return this.variants.values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .filter(variant -> variant.ofAlleles(alleleIdentifiers))
                .toList();
    }

    /**
     * Retrieves variants associated with the specified allele identifiers within a given range of positions.
     * <p>
     * This method filters the variants stored in the contig to include only those that are associated with at least one of the specified
     * allele identifiers and fall within the specified start and end positions (inclusive of start, exclusive of end). The resulting list
     * is unmodifiable.
     *
     * @param start             The 1-based start position of the range (inclusive).
     * @param end               The 1-based end position of the range (exclusive).
     * @param alleleIdentifiers An array of allele identifiers to filter the variants.
     * @return A {@link List} of {@link Variant} objects associated with the specified allele identifiers within the given range.
     */
    public List<Variant> getVariantsOfAllelesWithin(int start, int end, String... alleleIdentifiers) {
        return this.variants.subMap(start, end + 1).values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .filter(variant -> variant.ofAlleles(alleleIdentifiers))
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
     * Calculates the total number of active variants associated with this contig.
     * <p>
     * This method flattens the {@code variants} map into a stream of {@link Variant} objects. It then filters the stream to include only
     * those variants that are marked as active (i.e., have the {@code active} property set to {@code true}). The method counts the filtered
     * variants and returns the total count as an integer.
     *
     * @return The total number of active {@link Variant} objects associated with this contig.
     */
    public int getActiveVariantsCount() {
        return (int) this.variants.values().stream()
                .flatMap(variantMap -> variantMap.values().stream())
                .filter(v -> v.active)
                .count();
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
     * Removes a variant from the contig's variant map at the specified position and alternative base sequence.
     * <p>
     * This method checks if the {@code variants} map contains an entry for the given position. If no entry exists, the method returns
     * without performing any operation. If an entry exists, it removes the variant associated with the specified alternative base
     * sequence.
     * <p>
     * After removing the variant, the method checks if the map at the given position is empty. If it is, the position entry is also removed
     * from the {@code variants} map to maintain a clean structure.
     *
     * @param position    The 1-based position of the variant to remove.
     * @param alternative The alternative base sequence of the variant to remove.
     */
    void removeVariant(int position, String alternative) {
        if (!this.variants.containsKey(position)) return;
        this.variants.get(position).remove(alternative);
        if (this.variants.get(position).isEmpty()) this.variants.remove(position);
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
