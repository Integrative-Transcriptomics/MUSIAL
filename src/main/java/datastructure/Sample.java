package datastructure;

import utility.Constants;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;

/**
 * Represents a sample containing variant calls from a single biological sample.
 * <p>
 * This class extends {@link Attributable} to inherit functionality for managing attributes.
 * It provides fields and methods to store and manipulate variant calls, alleles, and other
 * sample-specific data. Each instance of this class is uniquely identified by its {@code name}.
 */
public class Sample extends Attributable {

    /**
     * The name or internal identifier of this sample.
     * <p>
     * This field uniquely identifies the sample within the context of the application.
     * It is a final field, meaning its value is immutable once assigned during the
     * construction of the {@link Sample} instance.
     */
    public final String name;

    /**
     * Hierarchical map structure to store variant calls.
     * <p>
     * This map organizes variant calls in a hierarchical structure:
     * <ul>
     *     <li>First level: The key is the name of the contig ({@link Contig#name}).</li>
     *     <li>Second level: The key is the position of the variant on the contig.</li>
     *     <li>Third level: The value is a string representing the variant call, formatted as:
     *         {@code CALL_INDEX;DP;GQ;REF_0:ALT_0:AD_0:PL_0,...}.
     *         <ul>
     *             <li>{@code CALL_INDEX}: Indicates the numeric index of the alternative allele with an
     *             optional prefix character of either {@code f} (low frequency) or {@code x} (low coverage).</li>
     *             <li>{@code DP}: The read depth at the variant site.</li>
     *             <li>{@code GQ}: The genotype quality score.</li>
     *             <li>{@code REF_0:.:AD_0:PL_0}: The reference allele, a placeholder (`.`), the allele depth, and the phred-scaled likelihood.</li>
     *             <li>{@code REF_1:ALT_1:AD_1:PL_1,...}: One or more alternate alleles, each with their respective
     *             reference allele, alternate allele, allele depth, and phred-scaled likelihoods.</li>
     *         </ul>
     *         The format follows the <a href="https://samtools.github.io/hts-specs/VCFv4.2.pdf">VCFv4.2</a> specification.
     *     </li>
     * </ul>
     */
    protected final HashMap<String, TreeMap<Integer, String>> variantCalls = new HashMap<>(2);

    /**
     * A map that assigns features to their corresponding alleles.
     * <p>
     * This {@link Map} stores the relationship between feature names and their associated allele identifiers.
     * The keys represent the names of the features, and the values represent the unique identifiers of the alleles.
     * This structure is used to track which allele is associated with each feature in the sample.
     */
    protected final Map<String, String> alleles;

    /**
     * Regular expression pattern to match variant call strings.
     * <p>
     * This pattern is designed to parse variant call strings that conform to the VCF specification.
     * The expected format includes fields separated by semicolons (`;`), with the following structure:
     * <ul>
     *   <li>{@code CALL_INDEX}: An optional prefix indicating the call index, which can be `f` (low frequency),
     *       `x` (low coverage), or a numeric index of the alternative allele.</li>
     *   <li>{@code DP}: The read depth at the variant site.</li>
     *   <li>{@code GQ}: The genotype quality score.</li>
     *   <li>{@code REF_0:.:AD_0:PL_0}: The reference allele, a placeholder (`.`), the allele depth, and the phred-scaled likelihood.</li>
     *   <li>{@code REF_1:ALT_1:AD_1:PL_1,...}: One or more alternate alleles, each with their respective
     *       reference allele, alternate allele, allele depth, and phred-scaled likelihoods.</li>
     * </ul>
     *
     * <pre>
     * Example match: {@code 1;13;99;TTC:.:0:585,TTC:T--:13:0}
     * </pre>
     */
    public static final Pattern variantCallPattern =
            Pattern.compile("([fx]?[0-9]+);[0-9]+;[0-9]+;([ACGTN-]+:(.|[ACGTN*-]+):[0-9]+:[0-9]+(,[ACGTN-]+:(.|[ACGTN*-]+):[0-9]+:[0-9]+)*)");

    /**
     * Constructs a new {@link Sample} instance with the specified name and initial capacity for the alleles map.
     * <p>
     * This constructor initializes a {@link Sample} object with the given name and allocates a {@link HashMap}
     * for the {@link #alleles} field with the specified initial capacity. The {@link #name} field is set to the
     * provided name, and the superclass constructor is invoked to initialize inherited properties.
     *
     * @param name     The name of the sample, used as its unique identifier.
     * @param capacity The expected initial capacity of the {@link #alleles} map.
     */
    protected Sample(String name, int capacity) {
        super();
        this.name = name;
        this.alleles = new HashMap<>(capacity);
    }

    /**
     * Associates a specific allele with a feature in this sample.
     * <p>
     * This method updates the {@link #alleles} map by setting the sequence type (allele)
     * for the specified feature. The feature is identified by its name, and the allele
     * is identified by its unique identifier.
     *
     * @param featureName The name of the feature ({@link Feature#name}) to associate with the allele.
     * @param alleleUid   The unique identifier of the allele ({@link SequenceType#name}) to set for the feature.
     */
    protected void setAllele(String featureName, String alleleUid) {
        this.alleles.put(featureName, alleleUid);
    }

    /**
     * Retrieves the allele identifier for the specified feature in this sample.
     * <p>
     * This method looks up the allele identifier associated with the given feature name
     * in the {@link #alleles} map. If no allele is set for the specified feature, the method
     * returns the default reference value defined in {@link Constants#synonymous}.
     *
     * @param featureName The name of the feature ({@link Feature#name}) to retrieve the associated allele identifier for.
     * @return The allele identifier ({@link Feature.Allele#_uid}) associated with the feature, or the reference value if not set.
     */
    public String getAllele(String featureName) {
        return this.alleles.getOrDefault(featureName, Constants.synonymous);
    }

    /**
     * Retrieves the entries of the alleles map for this sample.
     * <p>
     * This method returns a collection view of the mappings contained in the {@link #alleles} map.
     * Each entry in the collection represents a feature name and its associated allele identifier.
     * Modifications to the returned collection will reflect in the underlying map.
     *
     * @return A {@link Collection} of {@link Map.Entry} objects representing the entries in the {@link #alleles} map.
     */
    public Collection<Map.Entry<String, String>> getAlleles() {
        return this.alleles.entrySet();
    }

    /**
     * Retrieves the variant calls for the specified contig in this sample.
     * <p>
     * This method returns a {@link TreeMap} containing the variant calls for the given contig.
     * The keys in the map represent the positions of the variants on the contig, and the values
     * are the corresponding variant call strings. If no variant calls exist for the specified
     * contig, an empty {@link TreeMap} is returned.
     *
     * @param contig The name of the contig to retrieve the variant calls for.
     * @return A {@link TreeMap} where the keys are variant positions and the values are variant call strings.
     */
    public TreeMap<Integer, String> getVariantCalls(String contig) {
        return this.variantCalls.getOrDefault(contig, new TreeMap<>());
    }

    /**
     * Extracts the reference base character from the starting position of a variant call string.
     * <p>
     * This method processes a variant call string formatted as per the VCF specification and retrieves
     * the reference base character from the starting position. The call string is expected to follow
     * the structure defined in {@link Sample#variantCallPattern}, where fields are separated by semicolons,
     * commas, and colons.
     *
     * <pre>
     * Example call string: {@code 1;13;99;TTC:.:0:585,TTC:T--:13:0}
     * </pre>
     *
     * @param call The variant call string to process.
     * @return The reference base character from the starting position of the specified call.
     * @throws ArrayIndexOutOfBoundsException If the call string does not conform to the expected format.
     */
    public static String getReferenceOfCall(String call) {
        return call.split(Constants.SEMICOLON)[3].split(Constants.COMMA)[0].split(Constants.COLON)[0].substring(0, 1);
    }

    /**
     * Converts this sample to its string representation.
     * <p>
     * This method generates a string representation of the sample, including its name and attributes.
     * The attributes are formatted as key-value pairs separated by an equals sign (`=`) and delimited
     * by semicolons (`;`). If the last character of the generated string is a semicolon, it is removed
     * to ensure proper formatting.
     *
     * @return A {@link String} representing the sample, including its name and attributes.
     */
    public String toString() {
        StringBuilder sb = new StringBuilder(name).append("\t");
        this.getAttributes().forEach((key, value) ->
                sb.append(key).append(Constants.EQUAL).append(value).append(Constants.SEMICOLON)
        );
        if (sb.charAt(sb.length() - 1) == Constants.SEMICOLON.charAt(0)) {
            sb.setLength(sb.length() - 1);
        }
        return sb.toString();
    }

}
