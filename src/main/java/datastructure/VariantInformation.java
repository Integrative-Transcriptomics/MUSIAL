package datastructure;

import utility.Constants;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.Collectors;

/**
 * Stores information associated with a nucleotide variant.
 * <p>
 * The actual alternative content is not stored in this class!
 * <p>
 * This class represents a nucleotide variant, including its reference base content,
 * type (e.g., SNV, insertion, deletion), and occurrences in samples and alleles.
 * It provides methods to determine the type of the variant, check its canonical
 * or padded canonical status, and manage occurrences in samples and features.
 */
public class VariantInformation extends Attributable {

    /**
     * The reference base content of this variant.
     */
    public final String reference;

    /**
     * A mapping of occurrences of this variant in samples and alleles.
     * <p>
     * The `occurrence` map is structured as follows:
     * <ul>
     *     <li>The key is either {@link Attributable#sampleOccurrence} (representing sample occurrences)
     *         or the name of a {@link Feature} (representing feature occurrences).</li>
     *     <li>The value is a {@link HashSet} containing names of {@link Sample} or {@link SequenceType}
     *         associated with the key.</li>
     * </ul>
     * This structure allows efficient tracking of where the variant occurs in terms of samples and features.
     */
    protected final HashMap<String, HashSet<String>> occurrence = new HashMap<>(2);

    /**
     * The type of this variant (e.g., SNV, insertion, deletion).
     */
    public final Type type;

    /**
     * Enum representing the type of variant.
     */
    public enum Type {
        /**
         * Single Nucleotide Variant.
         */
        SNV,
        /**
         * An insertion of one or more nucleotides
         */
        INSERTION,
        /**
         * A deletion of one or more nucleotides
         */
        DELETION
    }

    /**
     * Constructor for {@link VariantInformation}.
     * <p>
     * Initializes the variant with its reference and alternative content, determining
     * its type (SNV, insertion, or deletion) based on the provided content.
     * <p>
     * This constructor checks if the reference and alternative content match any padded
     * canonical content type. If they do not, an {@link IllegalArgumentException} is thrown.
     *
     * @param referenceContent   The reference base content of the variant.
     * @param alternativeContent The alternative base content of the variant.
     * @throws IllegalArgumentException If the reference and alternative content do not match
     *                                  any padded canonical content type.
     */
    protected VariantInformation(String referenceContent, String alternativeContent) {
        super();
        if (isSubstitution(referenceContent, alternativeContent)) {
            this.type = Type.SNV;
        } else if (isInsertion(referenceContent, alternativeContent, true)) {
            this.type = Type.INSERTION;
        } else if (isDeletion(referenceContent, alternativeContent, true)) {
            this.type = Type.DELETION;
        } else {
            throw new IllegalArgumentException(
                    "Failed to construct `VariantInformation` instance. Contents (ref) %s and (alt) %s do not match any padded canonical content type."
                            .formatted(referenceContent, alternativeContent)
            );
        }
        this.reference = referenceContent;
        this.occurrence.put(Attributable.sampleOccurrence, new HashSet<>());
    }

    /**
     * Determines whether a variant is a substitution; i.e., both the reference and alternative base
     * content match a single base of {@link Constants#baseSymbols}.
     *
     * @param ref The reference base content.
     * @param alt The alternative base content.
     * @return {@code true} if the variant is a substitution, {@code false} otherwise.
     */
    public static boolean isSubstitution(String ref, String alt) {
        return ref.matches("^[%s]$".formatted(Constants.baseSymbols))
                && isSubstitution(alt);
    }

    /**
     * Determines whether a given alternative base content represents a substitution.
     * <p>
     * A substitution is defined as a single base from the set of valid nucleotide symbols
     * defined in {@link Constants#baseSymbols}.
     *
     * @param alt The alternative base content to check.
     * @return {@code true} if the alternative content represents a substitution, {@code false} otherwise.
     */
    public static boolean isSubstitution(String alt) {
        return alt.matches("^[%s]$".formatted(Constants.baseSymbols));
    }

    /**
     * Determines whether a variant is an insertion, i.e.,
     * <ul>
     *     <li>either the alternative base content is a string of any length of {@link Constants#baseSymbols}
     *     and the reference base content is a single base of {@link Constants#baseSymbols} followed by {@link Constants#gapString}s
     *     matching the alternative content's length (padded canonical),</li>
     *     <li>or the reference base content is a single base of {@link Constants#baseSymbols} and the alternative
     *     content is a string of any length of {@link Constants#baseSymbols} (un-padded canonical).</li>
     * </ul>
     *
     * @param ref    The reference base content.
     * @param alt    The alternative base content.
     * @param padded Whether the variant is padded by gap symbols.
     * @return {@code true} if the variant is an insertion, {@code false} otherwise.
     */
    public static boolean isInsertion(String ref, String alt, boolean padded) {
        if (padded) {
            return ref.length() == alt.length()
                    && ref.matches("^[%s]%s+$".formatted(Constants.baseSymbols, Constants.gapString))
                    && isInsertion(alt);
        } else {
            return ref.length() == 1
                    && alt.length() > 1
                    && ref.matches("^[%s]$".formatted(Constants.baseSymbols))
                    && alt.matches("^[%s]+$".formatted(Constants.baseSymbols));
        }
    }

    /**
     * Determines whether a variant is an insertion based on its alternative content.
     * <p>
     * This method checks if the alternative base content represents an insertion.
     * An insertion is defined as a string of at least two consecutive bases
     * from the set of valid nucleotide symbols defined in {@link Constants#baseSymbols}.
     *
     * @param alt The alternative base content to check.
     * @return {@code true} if the alternative content represents an insertion, {@code false} otherwise.
     */
    public static boolean isInsertion(String alt) {
        return alt.matches("^[%s]{2,}$".formatted(Constants.baseSymbols));
    }

    /**
     * Determines whether a variant is a deletion, i.e.,
     * <ul>
     *     <li>either the reference base content is a string of any length of {@link Constants#baseSymbols}
     *     and the alternative base content is a single base of {@link Constants#baseSymbols} followed by {@link Constants#gapString}s
     *     matching the reference content's length (padded canonical),</li>
     *     <li>or the reference base content is a string of any length of {@link Constants#baseSymbols} and the
     *     alternative content is a single base of {@link Constants#baseSymbols} (un-padded canonical).</li>
     * </ul>
     *
     * @param ref    The reference base content.
     * @param alt    The alternative base content.
     * @param padded Whether the variant is padded by gap symbols.
     * @return {@code true} if the variant is a deletion, {@code false} otherwise.
     */
    public static boolean isDeletion(String ref, String alt, boolean padded) {
        if (padded) {
            return ref.length() == alt.length()
                    && ref.matches("^[%s]+$".formatted(Constants.baseSymbols))
                    && isDeletion(alt);
        } else {
            return ref.length() > 1
                    && alt.length() == 1
                    && ref.matches("^[%s]+$".formatted(Constants.baseSymbols))
                    && alt.matches("^[%s]$".formatted(Constants.baseSymbols));
        }
    }

    /**
     * Determines whether a variant is a deletion based on its alternative content.
     * <p>
     * This method checks if the alternative base content represents a deletion.
     * A deletion is defined as a string that starts with a valid nucleotide base
     * (from {@link Constants#baseSymbols}) followed by one or more gap symbols
     * (defined in {@link Constants#gapString}).
     *
     * @param alt The alternative base content to check.
     * @return {@code true} if the alternative content represents a deletion, {@code false} otherwise.
     */
    public static boolean isDeletion(String alt) {
        return alt.matches("^[%s]%s+$".formatted(Constants.baseSymbols, Constants.gapString));
    }

    /**
     * Determines whether a variant is canonical.
     * <p>
     * A variant is canonical if it is:
     * <ul>
     *     <li>a single nucleotide variant (SNV) ({@link #isSubstitution}),</li>
     *     <li>an un-padded canonical insertion ({@link #isInsertion}), or</li>
     *     <li>an un-padded canonical deletion ({@link #isDeletion}).</li>
     * </ul>
     *
     * @param referenceContent   The reference base content.
     * @param alternativeContent The alternative base content.
     * @return {@code true} if the variant is canonical, {@code false} otherwise.
     */
    public static boolean isCanonicalVariant(String referenceContent, String alternativeContent) {
        return isSubstitution(referenceContent, alternativeContent)
                || isInsertion(referenceContent, alternativeContent, false)
                || isDeletion(referenceContent, alternativeContent, false);
    }

    /**
     * Determines whether a variant is padded canonical.
     * <p>
     * A variant is padded canonical if it is:
     * <ul>
     *     <li>a single nucleotide variant (SNV) ({@link #isSubstitution}),</li>
     *     <li>a padded canonical insertion ({@link #isInsertion}), or</li>
     *     <li>a padded canonical deletion ({@link #isDeletion}).</li>
     * </ul>
     *
     * @param referenceContent   The reference base content.
     * @param alternativeContent The alternative base content.
     * @return {@code true} if the variant is padded canonical, {@code false} otherwise.
     */
    public static boolean isPaddedCanonicalVariant(String referenceContent, String alternativeContent) {
        return isSubstitution(referenceContent, alternativeContent)
                || isInsertion(referenceContent, alternativeContent, true)
                || isDeletion(referenceContent, alternativeContent, true);
    }

    /**
     * Adds a sample occurrence to this variant.
     *
     * @param name The name of the sample to add.
     */
    protected void addSampleOccurrence(String name) {
        this.occurrence.get(Attributable.sampleOccurrence).add(name);
    }

    /**
     * Adds a feature occurrence to this variant.
     *
     * @param name The name of the feature to add.
     */
    protected void addFeatureOccurrence(String name) {
        this.occurrence.putIfAbsent(name, new HashSet<>(2));
    }

    /**
     * Adds an allele occurrence to this variant for a specific feature.
     *
     * @param featureName The name of the feature.
     * @param alleleUid   The unique identifier of the allele to add.
     */
    protected void addAlleleOccurrence(String featureName, String alleleUid) {
        addFeatureOccurrence(featureName);
        this.occurrence.get(featureName).add(alleleUid);
    }

    /**
     * Checks whether this variant has an occurrence in a sample or allele.
     *
     * @param of   Either {@code samples} or the name of a {@link Feature}.
     * @param name The name of the sample or allele to check for.
     * @return {@code true} if the sample or allele is associated with this variant, {@code false} otherwise.
     */
    public boolean hasOccurrence(String of, String name) {
        return this.occurrence.containsKey(of) && this.occurrence.get(of).contains(name);
    }

    /**
     * Checks whether this variant has an occurrence in a specific feature.
     *
     * @param name The name of the feature to check for.
     * @return {@code true} if the feature is associated with this variant, {@code false} otherwise.
     */
    public boolean hasOccurrence(String name) {
        return this.occurrence.containsKey(name);
    }

    /**
     * Retrieves the occurrences of this variant in samples.
     *
     * @return A {@link Collection} of sample names.
     */
    public Collection<String> getSampleOccurrence() {
        return this.occurrence.get(Attributable.sampleOccurrence);
    }

    /**
     * Retrieves the features associated with this variant.
     *
     * @return A {@link Collection} of feature names.
     */
    public Collection<String> getFeatureOccurrence() {
        return this.occurrence.keySet().stream().filter(key -> !key.equals(Attributable.sampleOccurrence)).collect(Collectors.toSet());
    }

    /**
     * Retrieves the reference base content of this variant.
     *
     * @param strip Whether to strip gap symbols from the reference base content.
     * @return The reference base content of this variant.
     */
    public String getReferenceBaseString(boolean strip) {
        if (strip) return this.reference.replaceAll("-", "");
        else return this.reference;
    }
}