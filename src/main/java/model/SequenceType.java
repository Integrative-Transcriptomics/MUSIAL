package model;

import htsjdk.samtools.util.Tuple;
import utility.Constants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents a sequence type with associated variants, samples, and attributes.
 * <p>
 * It extends the {@link Attributes} class to inherit functionality for managing attributes associated with the sample. This class is
 * extended by the {@link Feature.Allele} and {@link Feature.Proteoform} classes.
 * <p>
 * Sequence types are stored in the {@link Feature#alleles} and {@link Feature#proteoforms} properties of the model.
 */
public class SequenceType extends Attributes {

    /**
     * Unique identifier of this sequence type.
     * <p>
     * This field serves as the unique identifier for the feature and is used to reference it in the model.
     */
    public final String _id;

    /**
     * Variants defining this sequence type.
     * <p>
     * Navigable map of positions to canonical variants that define this sequence type.
     */
    protected final NavigableMap<Integer, String> variants = new TreeMap<>();

    /**
     * A set of sample names associated with this sequence type.
     * <p>
     * This set is used to track which samples correspond to this sequence type.
     */
    protected final HashSet<String> samples = new HashSet<>();

    /**
     * Constructs a new {@link SequenceType} instance with the specified identifier and variants.
     * <p>
     * This constructor initializes the sequence type with a unique identifier and a list of variants. The variants are sorted by their
     * positions in ascending order and then added to the {@link #variants} map.
     *
     * @param identifier The unique identifier for this sequence type.
     * @param variants   A list of {@link Tuple} objects representing the variants, where each tuple contains:
     *                   <ul>
     *                     <li>{@code a}: The position of the variant.</li>
     *                     <li>{@code b}: The variant's canonical base string.</li>
     *                   </ul>
     */
    public SequenceType(String identifier, List<Tuple<Integer, String>> variants) {
        super();
        this._id = identifier;
        variants.sort(Comparator.comparingInt(variant -> variant.a));
        variants.forEach(variant ->
                this.variants.put(variant.a, variant.b)
        );
    }

    /**
     * Associates a sample with this sequence type.
     *
     * @param sampleIdentifier The identifier of the sample to associate with this sequence type.
     */
    public void addRelation(String sampleIdentifier) {
        this.samples.add(sampleIdentifier);
    }

    /**
     * Retrieves a collection of sample identifiers that are related to this sequence type.
     *
     * @return A collection of sample identifiers that are related to this sequence type.
     */
    public Collection<String> getRelatedSamples() {
        return this.samples;
    }

    /**
     * Retrieves the count of occurrences of this sequence type.
     *
     * @return The number of unique identifiers associated with this sequence type.
     */
    public int getRelatedSamplesCount() {
        return this.samples.size();
    }

    /**
     * Checks if this sequence type is associated with an entity by its {@code identifier}.
     *
     * @param identifier Unique identifier to check for.
     * @return {@code true} if the entity is associated with this sequence type, {@code false} otherwise.
     */
    public boolean hasRelation(String identifier) {
        return this.samples.contains(identifier);
    }

    /**
     * Retrieves the variant at the specified position associated with this sequence type.
     * <p>
     * This method looks up the variant at the given position in the {@link #variants} map. If a variant exists at the specified position,
     * it returns the variant's canonical base string and returns {@code null} else.
     *
     * @param position The position to retrieve the variant for.
     * @return The variant's canonical base string at the specified position, or {@code null} if no variant is present.
     */
    public String getVariant(int position) {
        return this.variants.getOrDefault(position, null);
    }

    /**
     * Retrieves a list of variants associated with this sequence type.
     * <p>
     * This method converts the {@link #variants} map, which stores positions as keys and alternate alleles as values, into a list of
     * {@link Tuple} objects. Each tuple contains a position and its corresponding variant's canonical base string. If the map is empty, an
     * empty list is returned.
     *
     * @return A {@link List} of {@link Tuple} objects, where each tuple represents a variant with its position and alternate allele.
     */
    public List<Tuple<Integer, String>> getVariants() {
        // Convert the navigable map of variants to a list of tuples for easier access.
        if (this.variants.isEmpty()) {
            return Collections.emptyList();
        }
        List<Tuple<Integer, String>> variants = new ArrayList<>(this.variants.size());
        for (Map.Entry<Integer, String> entry : this.variants.entrySet()) {
            variants.add(new Tuple<>(entry.getKey(), entry.getValue()));
        }
        return variants;
    }

    /**
     * Checks if this sequence type has a variant at the specified position.
     *
     * @param position The position to check for a variant.
     * @return {@code true} if a variant exists at the specified position, {@code false} otherwise.
     */
    public boolean hasVariant(int position) {
        return this.variants.containsKey(position);
    }

    /**
     * Converts a list of variants to a string representation.
     * <p>
     * This method takes a list of {@link Tuple} objects, where each tuple contains a position and an alternate allele. It converts the list
     * into a string representation in the format {@code (POS0)(ALT0).(POS1)(ALT1)...}.
     *
     * @param variants A list of {@link Tuple} objects representing the variants.
     * @return A {@link String} representation of the variants in the format {@code (POS0)(ALT0).(POS1)(ALT1)...}.
     */
    public static String variantsToString(List<Tuple<Integer, String>> variants) {
        return variants.stream()
                .map(v -> v.a + v.b)
                .collect(Collectors.joining(Constants.dot));
    }

    /**
     * Computes the net shift in sequence length caused by variants.
     * <p>
     * This method calculates the cumulative effect of insertions and deletions on the sequence length. Each variant is analyzed to
     * determine whether it represents an insertion or a deletion:
     * <ul>
     *   <li>If the variant is an insertion, its length (number of bases minus one) is added to the net shift.</li>
     *   <li>If the variant is a deletion, its length (number of bases minus one) is subtracted from the net shift.</li>
     *   <li>Other types of variants do not affect the net shift.</li>
     * </ul>
     *
     * @param variants A list of {@link Tuple} objects, where each tuple contains:
     *                 <ul>
     *                   <li>{@code a}: The position of the variant (not used in this method).</li>
     *                   <li>{@code b}: The alternate allele of the variant.</li>
     *                 </ul>
     * @return The net shift in sequence length as an {@code int}.
     */
    protected static int computeLengthVariation(List<Tuple<Integer, String>> variants) {
        return variants.stream().mapToInt(variant -> {
            int length = variant.b.length() - 1;
            return Variant.isInsertion(variant.b) ? length :
                    Variant.isDeletion(variant.b) ? -length : 0;
        }).sum();
    }
}