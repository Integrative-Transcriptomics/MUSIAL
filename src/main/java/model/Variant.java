package model;

import htsjdk.samtools.util.Tuple;
import utility.Constants;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Representation of a nucleotide variant.
 * <p>
 * This class represents a nucleotide variant, including its reference base content, type (e.g., SNV, insertion, deletion), and occurrences
 * in samples and alleles. It provides methods to determine the type of the variant, check its canonical or padded canonical status. It
 * extends the {@link Attributes} class to inherit functionality for managing attributes associated with the variant.
 * <p>
 * In contrast to other entities in the model, this class does not implement an identifier, but the combination of {@code position},
 * {@code reference}, and {@code alternative} is used as such. Variants are stored in the {@link Contig#variants} property of the model.
 */
public class Variant extends Attributes {

    /**
     * The 1-based position of this variant on a contig.
     */
    public final int position;

    /**
     * The reference base content of this variant.
     */
    public final String reference;

    /**
     * The alternative base content of this variant.
     */
    public final String alternative;

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
     * The type of this variant (e.g., SNV, insertion, deletion).
     */
    public final Type type;

    /**
     * A set of sample names associated with this variant.
     * <p>
     * This set is used to track which samples have occurrences of this variant.
     */
    private final HashSet<String> samples = new HashSet<>();

    /**
     * A map of feature occurrences associated with this variant.
     * <p>
     * The keys are feature names, and the values are sets of allele identifiers associated with those features. The initial capacity is set
     * to one, as variants typically have a single feature associated with them.
     */
    private final HashMap<String, HashSet<String>> features = new HashMap<>(1);

    /**
     * Constructs a new {@link Variant} instance, based on the provided position, reference, and alternative content.
     * <p>
     * The constructor determines the type of variant based on the reference and alternative content as well as if the reference and
     * alternative content match any padded canonical content type. If they do not, an {@link IllegalArgumentException} is thrown.
     *
     * @param position    The 1-based position of the variant on a contig.
     * @param reference   The reference base content of the variant.
     * @param alternative The alternative base content of the variant.
     * @throws IllegalArgumentException If the reference and alternative content do not match any padded canonical content type.
     */
    protected Variant(int position, String reference, String alternative) {
        super();
        this.position = position;
        this.reference = reference;
        this.alternative = alternative;
        if (isSubstitution(reference, alternative)) {
            this.type = Type.SNV;
        } else if (isInsertion(reference, alternative, true)) {
            this.type = Type.INSERTION;
        } else if (isDeletion(reference, alternative, true)) {
            this.type = Type.DELETION;
        } else {
            throw new IllegalArgumentException(
                    ("Failed to construct `VariantInformation` instance. Contents (ref) %s and (alt) %s do not match any padded canonical" +
                            " content type.")
                            .formatted(reference, alternative)
            );
        }
    }

    /**
     * Associates a sample with this variant.
     *
     * @param sampleIdentifier The identifier of the sample to associate with this variant.
     */
    protected void addRelation(String sampleIdentifier) {
        this.samples.add(sampleIdentifier);
    }

    /**
     * Associates an allele and its parent feature with this variant.
     *
     * @param featureIdentifier The identifier of the feature to associate with this variant.
     * @param alleleIdentifier  The identifier of the allele to associate with the feature.
     */
    protected void addRelation(String featureIdentifier, String alleleIdentifier) {
        this.features.putIfAbsent(featureIdentifier, new HashSet<>(8));
        this.features.get(featureIdentifier).add(alleleIdentifier);
    }

    /**
     * Checks if this variant has a related sample of the given identifier.
     * <p>
     * This method checks if the provided identifier is present in the samples associated with this variant or in the features and their
     * alleles.
     *
     * @param identifier The identifier to check for occurrences in this variant.
     * @return {@code true} if the identifier is found in samples or features, {@code false} otherwise.
     */
    public boolean hasRelation(String identifier) {
        return this.samples.contains(identifier) || this.features.containsKey(identifier)
                || this.features.values().stream().anyMatch(alleles -> alleles.contains(identifier));
    }

    /**
     * Retrieves a collection of sample identifiers that have occurrences of this variant.
     *
     * @return A collection of sample identifiers that have occurrences of this variant.
     */
    public Collection<String> getRelatedSamples() {
        return this.samples;
    }

    /**
     * Retrieves a collection of tuples representing the feature and allele occurrences associated with this variant.
     * <p>
     * Each tuple contains a feature identifier and an allele identifier, representing the association of alleles with their parent
     * features.
     *
     * @return A collection of tuples representing the feature and allele occurrences.
     */
    public Collection<Tuple<String, String>> getRelatedAlleles() {
        HashSet<Tuple<String, String>> alleles = new HashSet<>();
        for (String feature : this.features.keySet()) {
            for (String allele : this.features.get(feature)) {
                alleles.add(new Tuple<>(feature, allele));
            }
        }
        return alleles;
    }

    /**
     * Retrieves the reference base content of this variant with all gap symbols removed.
     *
     * @return The reference base content with gap symbols removed.
     */
    public String getReferenceStripped() {
        return this.reference.replaceAll("-", "");
    }

    /**
     * Retrieves the alternative base content of this variant with all gap symbols removed.
     *
     * @return The alternative base content with gap symbols removed.
     */
    public String getAlternativeStripped() {
        return this.alternative.replaceAll("-", "");
    }

    /**
     * Determines whether a variant is a substitution; i.e., both the reference and alternative base content match a single base of
     * {@link Constants#baseSymbols}.
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
     * A substitution is defined as a single base from the set of valid nucleotide symbols defined in {@link Constants#baseSymbols}.
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
     *     and the reference base content is a single base of {@link Constants#baseSymbols} followed by {@link Constants#gap}s
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
                    && ref.matches("^[%s]%s+$".formatted(Constants.baseSymbols, Constants.gap))
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
     * This method checks if the alternative base content represents an insertion. An insertion is defined as a string of at least two
     * consecutive bases from the set of valid nucleotide symbols defined in {@link Constants#baseSymbols}.
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
     *     and the alternative base content is a single base of {@link Constants#baseSymbols} followed by {@link Constants#gap}s
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
     * This method checks if the alternative base content represents a deletion. A deletion is defined as a string that starts with a valid
     * nucleotide base (from {@link Constants#baseSymbols}) followed by one or more gap symbols (defined in {@link Constants#gap}).
     *
     * @param alt The alternative base content to check.
     * @return {@code true} if the alternative content represents a deletion, {@code false} otherwise.
     */
    public static boolean isDeletion(String alt) {
        return alt.matches("^[%s]%s+$".formatted(Constants.baseSymbols, Constants.gap));
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

}