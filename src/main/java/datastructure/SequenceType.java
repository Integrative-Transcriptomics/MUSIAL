package datastructure;

import htsjdk.samtools.util.Tuple;
import utility.Constants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents a sequence type with associated variants, occurrences, and attributes.
 * <p>
 * This class extends {@link Attributable} to manage metadata.
 * <p>
 * This class is extended by the {@link Feature.Allele} and {@link Feature.Proteoform} classes.
 *
 * @noinspection unused
 */
public class SequenceType extends Attributable {

    /**
     * Optional name to describe this sequence type.
     * <p>
     * This field stores a human-readable name for the sequence type. It is optional and can be
     * set to provide additional context or description for the sequence type. If not set, the
     * sequence type is identified solely by its unique identifier (_uid).
     * </p>
     */
    protected String name;

    /**
     * The unique identifier of this entity.
     * <p>
     * This field serves as a final and immutable unique identifier for the sequence type.
     * It is assigned during the construction of the {@link SequenceType} instance and cannot
     * be modified afterward. The identifier is used to uniquely distinguish this entity
     * from other sequence types.
     * </p>
     */
    public final String _uid;

    /**
     * Variants defining this entity.
     * <p>
     * This field stores a map of variants associated with this sequence type. The map is ordered
     * and navigable, with the keys representing the positions of the variants and the values
     * representing the alternate alleles. All variants must be represented in their canonical form.
     * </p>
     */
    protected final NavigableMap<Integer, String> variants;

    /**
     * A set representing the occurrences of this sequence type.
     * <p>
     * This field stores unique identifiers of entities where this sequence type occurs.
     * It is used to track and manage the presence of this sequence type across different contexts.
     * </p>
     */
    protected final HashSet<String> occurrence = new HashSet<>();

    /**
     * Constructs a new {@link SequenceType} instance with the specified unique identifier and variants.
     * <p>
     * This constructor initializes a {@link SequenceType} object with a unique identifier and a list
     * of variants. The variants are provided as a list of {@link Tuple} objects, where each tuple contains:
     * <ul>
     *   <li>The position of the variant (field {@code a} of the tuple).</li>
     *   <li>The alternate allele of the variant (field {@code b} of the tuple).</li>
     * </ul>
     * </p>
     *
     * <p>
     * The constructor populates the {@code variants} field, which is a {@link TreeMap}, by iterating
     * through the provided list of tuples. The positions and alternate alleles are extracted from
     * each tuple and added to the map, ensuring that the variants are stored in a sorted order
     * based on their positions.
     * </p>
     *
     * @param uid      The unique identifier for this sequence type.
     * @param variants A list of {@link Tuple} objects representing the variants associated
     *                 with this sequence type.
     */
    public SequenceType(String uid, List<Tuple<Integer, String>> variants) {
        super();
        this._uid = uid;
        this.variants = new TreeMap<>();
        for (Tuple<Integer, String> variant : variants) {
            this.variants.put(variant.a, variant.b);
        }
    }

    /**
     * Sets the name of this sequence type.
     *
     * @param name The name to set for this sequence type.
     */
    protected void setName(String name) {
        this.name = name;
    }

    /**
     * Retrieves the name of this sequence type.
     *
     * @return The name of this sequence type, or {@code null} if it has not been set.
     */
    public String getName() {
        return this.name;
    }

    /**
     * Adds an occurrence to this sequence type.
     *
     * @param identifier The unique identifier of the entity to add as an occurrence.
     */
    public void addOccurrence(String identifier) {
        this.occurrence.add(identifier);
    }

    /**
     * Retrieves the entities (by their unique identifiers) associated with this sequence type.
     *
     * @return A set of unique identifiers associated with this sequence type.
     */
    public HashSet<String> getOccurrence() {
        return this.occurrence;
    }

    /**
     * Checks if this sequence type is associated with an entity by its {@code identifier}.
     *
     * @param identifier Unique identifier to check for.
     * @return {@code true} if the entity is associated with this sequence type, {@code false} otherwise.
     */
    public boolean hasOccurrence(String identifier) {
        return this.occurrence.contains(identifier);
    }

    /**
     * Converts the occurrences of this sequence type to a comma-separated string.
     * <p>
     * This method joins all unique identifiers stored in the {@code occurrence} set
     * into a single string, separated by commas. It uses the delimiter defined in
     * {@link Constants#COMMA}.
     * </p>
     *
     * @return A {@link String} representation of the occurrences, separated by commas.
     */
    public String occurrenceAsString() {
        return String.join(Constants.COMMA, this.occurrence);
    }

    /**
     * Retrieves the variant at the specified position associated with this sequence type.
     * <p>
     * This method looks up the variant at the given position in the {@code variants} map.
     * If a variant exists at the specified position, it returns the corresponding alternate allele.
     * If no variant is found, it returns {@code null}.
     * </p>
     *
     * @param position The position to retrieve the variant for.
     * @return The alternate allele at the specified position, or {@code null} if no variant is present.
     */
    public String getVariant(int position) {
        return this.variants.getOrDefault(position, null);
    }

    /**
     * Retrieves the variants associated with this sequence type.
     * <p>
     * This method returns the map of variants that define this sequence type. The map is navigable,
     * with the keys representing the positions of the variants and the values representing the
     * alternate base sequences. The returned map is immutable and reflects the canonical form
     * of the variants.
     * </p>
     *
     * @return A {@link NavigableMap} of variants, where the keys are positions and the values are
     * the alternate base sequences.
     */
    public NavigableMap<Integer, String> getVariants() {
        return this.variants;
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
     * Checks if this sequence type has a specific variant at the specified position.
     *
     * @param position The position to check for a variant.
     * @param content  The content of the variant to check for.
     * @return {@code true} if the specified variant exists at the position, {@code false} otherwise.
     */
    public boolean hasVariant(int position, String content) {
        return Objects.equals(getVariant(position), content);
    }

    /**
     * Converts the variants of this sequence type to a string representation.
     * <p>
     * This method uses {@link #variantsAsString(Map)} to convert the variants map
     * into a string representation in the format {@code (POS0)(ALT0).(POS1)(ALT1)...}.
     * </p>
     *
     * @return A {@link String} representation of the variants in the format
     * {@code (POS0)(ALT0).(POS1)(ALT1)...}.
     */
    public String variantsAsString() {
        return variantsAsString(this.variants);
    }

    /**
     * Returns a string representation of this sequence type in the format {@code identifier    attributes  variants}.
     * <p>
     * The string representation includes:
     * <ul>
     *   <li>The identifier, which is either the {@code name} (if set) or the unique identifier {@code _uid}.</li>
     *   <li>The attributes of this sequence type, formatted using {@link Attributable#attributesAsString()}.</li>
     *   <li>The variants associated with this sequence type, formatted using {@link #variantsAsString()}.</li>
     * </ul>
     * </p>
     *
     * @return A {@link String} representing this sequence type in the format {@code identifier    attributes  variants}.
     */
    public String toString() {
        return (Objects.isNull(name) ? "%s\t".formatted(_uid) : "%s\t".formatted(name)) +
                this.attributesAsString() + "\t" +
                this.variantsAsString();
    }

    /**
     * Converts a map of variants to a string representation.
     * <p>
     * This method takes a map of variants, where the keys are positions and the values are
     * alternate alleles. It converts the map into a string representation in the format
     * {@code (POS0)(ALT0).(POS1)(ALT1)...}.
     * </p>
     *
     * @param variants A map of variants, where the keys are positions and the values are alternate alleles.
     * @return A {@link String} representation of the variants in the format {@code (POS0)(ALT0).(POS1)(ALT1)...}..
     */
    public static String variantsAsString(Map<Integer, String> variants) {
        return variants.entrySet().stream()
                .map(e -> e.getKey() + e.getValue())
                .collect(Collectors.joining(Constants.DOT));
    }

    /**
     * Converts a list of variants to a string representation.
     * <p>
     * This method takes a list of {@link Tuple} objects, where each tuple contains a position
     * and an alternate allele. It converts the list into a string representation in the format
     * {@code (POS0)(ALT0).(POS1)(ALT1)...}.
     * </p>
     *
     * @param variants A list of {@link Tuple} objects representing the variants.
     * @return A {@link String} representation of the variants in the format {@code (POS0)(ALT0).(POS1)(ALT1)...}.
     */
    public static String variantsAsString(List<Tuple<Integer, String>> variants) {
        return variants.stream()
                .map(e -> e.a + e.b)
                .collect(Collectors.joining(Constants.DOT));
    }

    /**
     * Computes the net shift in sequence length caused by variants.
     * <p>
     * This method calculates the cumulative effect of insertions and deletions
     * on the sequence length. Each variant is analyzed to determine whether it
     * represents an insertion or a deletion:
     * <ul>
     *   <li>If the variant is an insertion, its length (number of bases minus one) is added to the net shift.</li>
     *   <li>If the variant is a deletion, its length (number of bases minus one) is subtracted from the net shift.</li>
     *   <li>Other types of variants do not affect the net shift.</li>
     * </ul>
     * </p>
     *
     * @param variants A list of {@link Tuple} objects, where each tuple contains:
     *                 <ul>
     *                   <li>{@code a}: The position of the variant (not used in this method).</li>
     *                   <li>{@code b}: The alternate allele of the variant.</li>
     *                 </ul>
     * @return The net shift in sequence length as an {@code int}.
     */
    public static int computeLengthVariation(List<Tuple<Integer, String>> variants) {
        return variants.stream().mapToInt(variant -> {
            int length = variant.b.length() - 1;
            return VariantInformation.isInsertion(variant.b) ? length :
                    VariantInformation.isDeletion(variant.b) ? -length : 0;
        }).sum();
    }

    /**
     * Generates a FASTA header for the sequence type.
     * <p>
     * This method constructs a FASTA header string for the sequence type using the provided feature name
     * and sequence identifier. The header includes an identifier in the format {@code lcl|<featureName>_<sequenceIdentifier>}.
     * Additionally, it appends optional properties to the header if they are present as attributes:
     * <ul>
     *   <li>{@code allelic_frequency}: The allelic frequency of the sequence type.</li>
     *   <li>{@code so_effects}: Sequence ontology effects associated with the sequence type.</li>
     * </ul>
     * </p>
     *
     * @param featureName        The name of the feature to include in the FASTA header.
     * @param sequenceIdentifier The identifier of the sequence to include in the FASTA header.
     * @return A {@link String} representing the FASTA header, including the identifier and optional properties.
     */
    public String getFastaHeader(String featureName, String sequenceIdentifier) {
        // Construct the identifier for the FASTA header in the format "lcl|<featureName>_<sequenceIdentifier>"
        String id = "lcl|%s_%s".formatted(featureName, sequenceIdentifier);

        // Initialize a collection to store optional properties for the FASTA header
        Collection<String> properties = new HashSet<>();

        // Add the "allelic_frequency" property if it exists and is not empty
        if (hasAttribute("allelic_frequency") && !getAttribute("allelic_frequency").isEmpty())
            properties.add("allelic_frequency=%s".formatted(getAttribute("allelic_frequency")));

        // Add the "so_effects" property if it exists and is not empty
        if (hasAttribute(Constants.$SequenceType_effects) && !getAttribute(Constants.$SequenceType_effects).isEmpty())
            properties.add("so_effects=%s".formatted(getAttribute(Constants.$SequenceType_effects)));

        // Construct and return the FASTA header, appending properties if they exist
        return ">%s%s\n".formatted(id, properties.isEmpty() ? Constants.EMPTY : " " + String.join(" ", properties));
    }
}