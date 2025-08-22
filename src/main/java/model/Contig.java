package model;

import htsjdk.samtools.util.Tuple;
import utility.Constants;
import utility.IO;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Representation of a genomic location.
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
    private final TreeMap<Integer, ArrayList<Variant>> variants;

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
    protected transient HashMap<Tuple<Integer, Integer>, String> cache;

    /**
     * Constructs a new {@link Contig} instance with the specified sequence (optionally empty).
     * <p>
     * This constructor initializes a contig with its identifier and nucleotide sequence. The {@link #variants} map to store variant
     * information and the {@link #cache} map to cache subsequences for optimized retrieval are initialized as empty instances. If a
     * non-empty sequence is provided, it is compressed using GZIP to reduce storage requirements.
     *
     * @param identifier The unique identifier for this contig, which is used to reference it in the context of a {@link Storage} instance.
     * @param sequence   The nucleotide sequence of the contig, stored as a GZIP-compressed string.
     * @throws IOException If an error occurs during the compression of the sequence.
     */
    protected Contig(String identifier, String sequence) throws IOException {
        super(); // Call the constructor of the parent class (Attributes).
        this._id = identifier; // Assign the unique identifier to the contig.

        String compressedSequence; // Variable to store the compressed sequence.
        int length; // Variable to store the length of the sequence.

        // Check if the sequence is non-null and not empty.
        if (Objects.nonNull(sequence) && !sequence.isEmpty()) {
            compressedSequence = IO.gzipCompress(sequence); // Compress the sequence using GZIP.
            length = sequence.length(); // Calculate the length of the sequence.
        } else {
            compressedSequence = Constants.empty; // Assign an empty string if the sequence is null or empty.
            length = 0; // Set the length to 0 for an empty sequence.
        }

        // Add the length of the sequence as an attribute to the contig.
        addAttribute(Constants.$Contig_length, String.valueOf(length));

        this.sequence = compressedSequence; // Store the compressed sequence.
        this.variants = new TreeMap<>(Integer::compare); // Initialize the map to store variants.
        this.cache = new HashMap<>(); // Initialize the cache for subsequences.
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
            return Constants.empty;
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
    public String getSubsequence(int start, int end) throws IOException {
        if (hasSequence()) {
            Tuple<Integer, Integer> cacheKey = new Tuple<>(start, end);
            if (cache.containsKey(cacheKey)) {
                return cache.get(cacheKey);
            } else {
                String subsequence = getSequence().substring(start - 1, end);
                cache.put(cacheKey, subsequence);
                return subsequence;
            }
        } else {
            return Constants.empty;
        }
    }

    /**
     * Adds a new {@link Variant} to this contig at the specified position with the given alternative bases and reference bases.
     * <p>
     * This method creates a new {@link Variant} instance with the specified alternative bases and reference bases, and adds it to the
     * {@link #variants} map.
     *
     * @param position    The 1-based position of the variant on the contig.
     * @param reference   The reference base sequence of the variant.
     * @param alternative The alternative base sequence of the variant.
     * @return The newly created {@link Variant} instance.
     */
    protected Variant addVariant(int position, String reference, String alternative) {
        Variant variant = new Variant(position, reference, alternative);
        this.variants.computeIfAbsent(position, k -> new ArrayList<>()).add(variant);
        return variant;
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
        return this.variants.get(position).stream().findFirst().filter(v -> v.alternative.equals(alternative)).orElse(null);
    }

    /**
     * Retrieves all variants associated with this contig.
     * <p>
     * This method flattens the {@code variants} map, which organizes variants by their positions, into a single list of {@link Variant}
     * objects.
     *
     * @return A {@link List} containing all {@link Variant} objects associated with this contig.
     */
    public List<Variant> getVariants() {
        // Flatten the variants map into a single list.
        return this.variants.values().stream()
                .flatMap(Collection::stream)
                .collect(Collectors.toList());
    }

    /**
     * Retrieves variants within the specified range of positions.
     * <p>
     * This method retrieves variants from the {@code variants} map that fall within the specified start and end positions (inclusive of
     * start, exclusive of end). The resulting variants are flattened into a single list.
     *
     * @param start The 1-based start position of the range (inclusive).
     * @param end   The 1-based end position of the range (exclusive).
     * @return A {@link List} of {@link Variant} objects within the specified range.
     */
    public List<Variant> getVariants(int start, int end) {
        // Retrieve variants within the specified range and flatten them into a single list.
        return this.variants.subMap(start, end + 1).values().stream()
                .flatMap(Collection::stream)
                .collect(Collectors.toList());
    }

    /**
     * Retrieves variants based on the provided set of sample names.
     * <p>
     * This method filters the {@code variants} map to include only those variants that are associated with at least one of the specified
     * sample names.
     *
     * @param relations A variable-length array of identifiers to filter the variants.
     * @return A {@link List} of {@link Variant} objects that match the specified sample names.
     */
    public List<Variant> getVariants(String... relations) {
        // Filter variants based on the provided set of sample names.
        return this.variants.values().stream()
                .flatMap(Collection::stream)
                .filter(variant -> Arrays.stream(relations).anyMatch(variant::hasRelation))
                .collect(Collectors.toList());
    }

    /**
     * Retrieves variants within the specified range of positions and filters them based on the provided set of sample names.
     * <p>
     * This method retrieves variants from the {@code variants} map that fall within the specified start and end positions (inclusive of
     * start, exclusive of end), and filters them to include only those variants that are associated with at least one of the specified
     * sample names.
     *
     * @param start     The 1-based start position of the range (inclusive).
     * @param end       The 1-based end position of the range (exclusive).
     * @param relations A variable-length array of identifiers to filter the variants.
     * @return A {@link List} of {@link Variant} objects within the specified range and matching the specified sample names.
     */
    public List<Variant> getVariants(int start, int end, String... relations) {
        // Retrieve variants within the specified range and filter them based on the provided set of sample names.
        return this.variants.subMap(start, end + 1).values().stream()
                .flatMap(Collection::stream)
                .filter(variant -> Arrays.stream(relations).anyMatch(variant::hasRelation))
                .collect(Collectors.toList());
    }

    /**
     * Calculates the total number of variants associated with this contig.
     * <p>
     * This method iterates through the {@code variants} map, which organizes variants by their positions, and sums up the sizes of all the
     * lists of variants. The result represents the total count of {@link Variant} objects stored in this contig.
     *
     * @return The total number of {@link Variant} objects associated with this contig.
     */
    public int getVariantsCount() {
        return this.variants.values().stream().mapToInt(ArrayList::size).sum();
    }

}
