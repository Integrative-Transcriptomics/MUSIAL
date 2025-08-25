package model;

import htsjdk.samtools.util.Tuple;
import org.apache.commons.lang3.tuple.MutableTriple;
import utility.Constants;

import java.util.*;
import java.util.regex.Pattern;

/**
 * Representation of a collection of variant calls from a single biological sample.
 * <p>
 * This class provides a structure to store and manage variant calls associated with a biological sample. Furthermore, it allows the
 * association with an allele per feature. It extends the {@link Attributes} class to inherit functionality for managing attributes
 * associated with the sample.
 * <p>
 * Samples are stored in the {@link Storage#samples} property of the model.
 */
public class Sample extends Attributes {

    /**
     * Unique identifier of this sample.
     * <p>
     * This field serves as the unique identifier for the sample and is used to reference it in the model. This should be at best a database
     * identifier, such as a NCBI accession number, or strain identifier.
     */
    public final String _id;

    /**
     * A map that assigns features to their corresponding alleles.
     * <p>
     * This {@link Map} stores the relationship between feature names and their associated allele identifiers. The keys represent the names
     * of the features, and the values represent the unique identifiers of the alleles. This structure is used to track which allele is
     * associated with each feature in the sample.
     */
    private final Map<String, String> alleles;

    /**
     * Regular expression pattern to match variant call strings.
     * <p>
     * This pattern is designed to parse variant call strings that conform the expected format {@code CI;DP;CE;REF_0:ALT_0:AD_0,...}.
     * <ul>
     *     <li>{@code CI}: Indicates whether the reference (0) or an alternative (1) is called. An optional prefix character of
     *     either {@code f} (low frequency), {@code x} (low coverage), or {@code u} (missing upstream deletion) can be added.</li>
     *     <li>{@code DP}: The read depth at the variant site.</li>
     *     <li>{@code CE}: The calls normalized information entropy.</li>
     *     <li>{@code REF_i:ALT_i:AD_i}: The reference content, alternative content - a dot (.) in case of the reference call -
     *     and the allelic depth; separated by commas and sorted by read support/frequency.</li>
     * </ul>
     * <pre>
     * Example: {@code 1;49;0.592;G:C:42,G:.:7}
     * </pre>
     */
    public static final Pattern callPattern =
            Pattern.compile("([fxu]?[01]);[0-9]+;[.0-9]+;([ACGTN-]+:(.|[ACGTN*-]+):[0-9]+(,[ACGTN-]+:(.|[ACGTN*-]+):[0-9]+)*)");

    /**
     * Hierarchical map structure to store variant calls.
     * <p>
     * This map organizes variant calls in a two-layered map structure:
     * <ul>
     *     <li>Outer key: The key is the id of the contig ({@link Contig#_id}).</li>
     *     <li>Inner key: The key is the position of the variant on the contig.</li>
     *     <li>Inner value: The value is a string representing the variant call, formatted as:
     *         {@code CI;DP;CE;REF_0:ALT_0:AD_0,...}.
     *         <ul>
     *             <li>{@code CI}: Indicates whether the reference (0) or an alternative (1) is called. An optional prefix character of
     *             either {@code f} (low frequency), {@code x} (low coverage), or {@code u} (missing upstream deletion) can be added.</li>
     *             <li>{@code DP}: The read depth at the variant site.</li>
     *             <li>{@code CE}: The calls normalized information entropy.</li>
     *             <li>{@code REF_i:ALT_i:AD_i}: The reference content, alternative content - a dot (.) in case of the reference call -
     *             and the allelic depth; separated by commas and sorted by read support/frequency.</li>
     *         </ul>
     *     </li>
     * </ul>
     */
    private final HashMap<String, TreeMap<Integer, String>> calls;

    /**
     * Constructs a new {@link Sample} instance with the specified id and initial capacity for the alleles map.
     * <p>
     * This constructor initializes a {@link Sample} object with the given id and allocates a {@link HashMap} for the {@link #alleles} field
     * with the specified initial capacity. The {@link #_id} field is set to the provided id, and the superclass constructor is invoked to
     * initialize inherited properties.
     *
     * @param identifier The id of the sample, used as its unique identifier.
     * @param capacity   The expected initial capacity of the {@link #alleles} map.
     */
    protected Sample(String identifier, int capacity) {
        super();
        this._id = identifier;
        this.calls = new HashMap<>(2);
        this.alleles = new HashMap<>(capacity);
    }

    /**
     * Associates a specific allele with a feature in this sample.
     * <p>
     * This method updates the {@link #alleles} map by setting the sequence type (allele) for the specified feature. The feature is
     * identified by its id, and the allele is identified by its unique identifier.
     *
     * @param featureIdentifier The id of the feature ({@link Feature#name}) to associate with the allele.
     * @param alleleIdentifier  The unique identifier of the allele ({@link SequenceType#_id}) to set for the feature.
     */
    protected void addRelation(String featureIdentifier, String alleleIdentifier) {
        this.alleles.put(featureIdentifier, alleleIdentifier);
    }

    /**
     * Retrieves the allele associated with a specific feature in this sample.
     * <p>
     * This method looks up the {@link #alleles} map to find the allele associated with the given feature identifier. If no association
     * exists, it returns the default reference allele defined in {@link Constants#reference}.
     *
     * @param featureIdentifier The unique identifier of the feature ({@link Feature#name}) to retrieve the associated allele for.
     * @return The unique identifier of the allele ({@link SequenceType#_id}) associated with the feature, or the default reference allele
     * if no association exists.
     */
    public String getRelatedAllele(String featureIdentifier) {
        return this.alleles.getOrDefault(featureIdentifier, Constants.reference);
    }

    /**
     * Retrieves all feature-allele associations in this sample.
     * <p>
     * This method converts the {@link #alleles} map, which stores feature names as keys and their associated allele identifiers as values,
     * into a collection of {@link Tuple} objects. Each tuple contains a feature name and its corresponding allele identifier.
     *
     * @return A {@link Collection} of {@link Tuple} objects, where each tuple represents a feature-allele association.
     */
    public Collection<Tuple<String, String>> getRelatedAlleles() {
        return this.alleles.entrySet().stream()
                .map(entry -> new Tuple<>(entry.getKey(), entry.getValue()))
                .toList();
    }

    /**
     * Retrieves the number of alternate alleles in this sample.
     * <p>
     * This method returns the size of the {@link #alleles} map, which represents the number of unique alleles associated with features in
     * this sample. This corresponds to the number of non-reference alleles present in the sample.
     *
     * @return The number of alleles in this sample.
     */
    public int getRelatedAllelesCount() {
        return this.alleles.size();
    }

    /**
     * Adds a variant call to the specified contig at the given position.
     * <p>
     * This method validates the format of the provided variant call string using the {@link #callPattern} regular expression. If the format
     * is invalid, an {@link IllegalArgumentException} is thrown. If the format is valid, the method ensures that the contig exists in the
     * {@link #calls} map, and then adds the variant call at the specified position.
     *
     * @param contigIdentifier The unique identifier of the contig to which the variant call belongs.
     * @param position         The position of the variant on the contig.
     * @param callString       The variant call string, which must conform to the expected format defined by {@link #callPattern}.
     * @throws IllegalArgumentException If the provided call string does not match the expected format.
     */
    protected void addVariantCall(String contigIdentifier, int position, String callString) {
        if (!callPattern.matcher(callString).matches())
            throw new IllegalArgumentException("Invalid call format %s for sample %s at position %d on contig %s. Expected: %s."
                    .formatted(callString, this._id, position, contigIdentifier, Sample.callPattern.pattern()));
        calls.computeIfAbsent(contigIdentifier, k -> new TreeMap<>()); // Ensure the contig exists
        calls.get(contigIdentifier).put(position, callString);
    }

    /**
     * Checks if there are any variant calls associated with this sample.
     * <p>
     * This method verifies whether the {@link #calls} map contains any entries, indicating the presence of variant calls.
     *
     * @return {@code true} if there are variant calls, {@code false} otherwise.
     */
    public boolean hasVariantCall() {
        return !this.calls.isEmpty();
    }

    /**
     * Checks if there are any variant calls for a specific contig.
     * <p>
     * This method checks whether the {@link #calls} map contains the specified contig identifier as a key and if the associated map of
     * positions is not empty.
     *
     * @param contigIdentifier The unique identifier of the contig to check for variant calls.
     * @return {@code true} if there are variant calls for the specified contig, {@code false} otherwise.
     */
    public boolean hasVariantCall(String contigIdentifier) {
        return this.calls.containsKey(contigIdentifier) && !this.calls.get(contigIdentifier).isEmpty();
    }

    /**
     * Retrieves the variant call at a specific position on a given contig.
     * <p>
     * This method checks if there are variant calls for the specified contig. If no calls exist, it returns an empty string defined by
     * {@link Constants#empty}. Otherwise, it retrieves the variant call at the specified position, or returns {@link Constants#empty} if no
     * call exists at that position.
     *
     * @param contigIdentifier The unique identifier of the contig.
     * @param position         The position of the variant on the contig.
     * @return The variant call string at the specified position, or {@link Constants#empty} if no call exists.
     */
    public String getVariantCall(String contigIdentifier, int position) {
        if (!hasVariantCall(contigIdentifier)) return Constants.empty;
        return calls.get(contigIdentifier).getOrDefault(position, Constants.empty);
    }

    /**
     * Retrieves a list of variant calls associated with this sample.
     * <p>
     * This method processes the {@link #calls} map, which organizes variant calls by contig identifiers and positions, and converts it into
     * a list of {@link MutableTriple} objects. Each triple contains:
     * <ul>
     *   <li>The contig identifier as a {@link String}.</li>
     *   <li>The position of the variant as an {@link Integer}.</li>
     *   <li>The variant call string as a {@link String}.</li>
     * </ul>
     * The resulting list provides a flattened representation of all variant calls in the sample.
     *
     * @return A {@link List} of {@link MutableTriple} objects, where each triple represents a variant call with its contig identifier,
     * position, and call string.
     */
    public List<MutableTriple<String, Integer, String>> getVariantCalls() {
        return this.calls.entrySet().stream()
                .flatMap(entry -> entry.getValue().entrySet().stream()
                        .map(callEntry -> new MutableTriple<>(entry.getKey(), callEntry.getKey(), callEntry.getValue())))
                .toList();
    }

    /**
     * Retrieves a list of variant calls for a specific contig.
     * <p>
     * This method processes the {@link #calls} map to extract variant calls associated with the specified contig identifier. It converts
     * the entries into a list of {@link MutableTriple} objects, where each triple contains:
     * <ul>
     *   <li>The contig identifier as a {@link String}.</li>
     *   <li>The position of the variant as an {@link Integer}.</li>
     *   <li>The variant call string as a {@link String}.</li>
     * </ul>
     * If the contig identifier does not exist in the {@link #calls} map, an empty {@link TreeMap} is used as a default.
     *
     * @param contigIdentifier The unique identifier of the contig to retrieve variant calls for.
     * @return A {@link List} of {@link MutableTriple} objects representing the variant calls for the specified contig.
     */
    public List<MutableTriple<String, Integer, String>> getVariantCalls(String contigIdentifier) {
        return this.calls.getOrDefault(contigIdentifier, new TreeMap<>()).entrySet().stream()
                .map(callEntry -> new MutableTriple<>(contigIdentifier, callEntry.getKey(), callEntry.getValue()))
                .toList();
    }

    /**
     * Parses a variant call string into a list of alternative alleles.
     * <p>
     * This method splits the provided variant call string into its components and extracts the alternative alleles. Each alternative is
     * represented as a {@link MutableTriple} containing the reference content, alternative content, and allelic depth.
     *
     * @param callString The variant call string to parse.
     * @return A {@link List} of {@link MutableTriple} objects representing the alternative alleles.
     */
    public static List<MutableTriple<String, String, Integer>> callStringToAlternatives(String callString) {
        List<MutableTriple<String, String, Integer>> alternatives = new ArrayList<>();
        if (callString.isEmpty()) return alternatives;
        for (String alternative : callString.split(Constants.semicolon)[3].split(Constants.comma)) {
            String[] parts = alternative.split(Constants.colon);
            alternatives.add(new MutableTriple<>(parts[0], parts[1], Integer.parseInt(parts[2])));
        }
        return alternatives;
    }

    /**
     * Extracts the reference allele from a variant call string.
     * <p>
     * This method parses the provided variant call string and retrieves the reference allele from the first alternative entry.
     *
     * @param callString The variant call string to parse.
     * @return The reference allele as a {@link String}.
     */
    public static String callStringToReference(String callString) {
        return callString.split(Constants.semicolon)[3].split(Constants.comma)[0].split(Constants.colon)[0].substring(0, 1);
    }

    /**
     * Converts this sample to its string representation.
     * <p>
     * This method generates a string representation of the sample, including its id and attributes. The attributes are formatted as
     * key-value pairs separated by an equals sign (`=`) and delimited by semicolons (`;`). If the last character of the generated string is
     * a semicolon, it is removed to ensure proper formatting.
     *
     * @return A {@link String} representing the sample, including its id and attributes.
     */
    public String toString() {
        StringBuilder sb = new StringBuilder(_id).append("\t");
        this.getAttributes().forEach((key, value) ->
                sb.append(key).append(Constants.equal).append(value).append(Constants.semicolon)
        );
        if (sb.charAt(sb.length() - 1) == Constants.semicolon.charAt(0)) {
            sb.setLength(sb.length() - 1);
        }
        return sb.toString();
    }

}
