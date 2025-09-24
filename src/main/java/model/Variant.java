package model;

import htsjdk.samtools.util.Tuple;
import util.Bio;
import util.Constants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents a nucleotide variant.
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
     * Represents a simplified variant stub with position and alternative allele information.
     * <p>
     * This record encapsulates the position and alternative allele of a variant and provides methods for generating its string
     * representation, computing its hash code, and checking equality.
     * <p>
     * Variant stubs are used to represent variants in a simplified form, to be precisely identify, they constitute the index of a variant
     * in a contig, as variants do not have their own identifiers.
     *
     * @param position    The 1-based position of the variant on a contig.
     * @param alternative The alternative allele of the variant.
     */
    public record Stub(int position, String alternative) {

        /**
         * Converts the variant stub to its string representation.
         * <p>
         * The string representation is formatted as "position + alternative".
         *
         * @return A {@link String} representing the variant stub.
         */
        public String toString() {
            return "%d%s%s".formatted(position, Constants.GREATER_THAN, alternative);
        }

        /**
         * Computes the hash code for this variant stub.
         * <p>
         * The hash code is calculated based on the string representation of the variant stub.
         *
         * @return The hash code of the variant stub.
         */
        public int hashCode() {
            return this.toString().hashCode();
        }

        /**
         * Compares this variant stub to another object for equality.
         * <p>
         * This method checks if the provided object is the same instance as this object. If not, it verifies that the object is of the same
         * class and compares their string representations.
         *
         * @param obj The object to compare with this {@link Stub} instance.
         * @return {@code true} if the objects are the same instance or if their string representations are equal; {@code false} otherwise.
         */
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (obj == null || getClass() != obj.getClass()) return false;
            Stub that = (Stub) obj;
            return Objects.equals(this.toString(), that.toString());
        }

    }

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
    private final Map<String, String> samples = new HashMap<>(100);

    /**
     * A map of feature occurrences associated with this variant.
     * <p>
     * The keys are feature names, and the values are sets of allele identifiers associated with those features. The initial capacity is set
     * to one, as variants typically have a single feature associated with them.
     */
    private final Map<String, Set<String>> features = new HashMap<>(1);

    /**
     * Indicates whether this variant is active.
     * <p>
     * This boolean flag is used to mark variants that are newly identified and not present in an existing storage or modified. It is set to
     * {@code true} for newly created or touched entries, but will not be serialized. During deserialization, it is assumed that all
     * variants are idle and should be set to {@code false} (see {@link Storage#typeAdapter()}).
     */
    protected transient boolean active;

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
        if (Bio.isSubstitution(reference, alternative)) {
            this.type = Type.SNV;
        } else if (Bio.isInsertion(reference, alternative, true)) {
            this.type = Type.INSERTION;
        } else if (Bio.isDeletion(reference, alternative, true)) {
            this.type = Type.DELETION;
        } else {
            throw new IllegalArgumentException(
                    ("Failed to construct `VariantInformation` instance. Contents (ref) %s and (alt) %s do not match any padded canonical" +
                            " content type.")
                            .formatted(reference, alternative)
            );
        }
        this.active = true;
    }

    /**
     * Checks if the variant is associated with a specific sample.
     *
     * @param sampleIdentifier The unique identifier of the sample to check.
     * @return {@code true} if the variant is associated with the given sample identifier, {@code false} otherwise.
     */
    public boolean ofSample(String sampleIdentifier) {
        return this.samples.containsKey(sampleIdentifier);
    }

    /**
     * Checks if the variant is associated with any of the specified samples.
     *
     * @param sampleIdentifiers A collection of sample identifiers to check.
     * @return {@code true} if the variant is associated with at least one of the given sample identifiers, {@code false} otherwise.
     */
    public boolean ofSamples(Collection<String> sampleIdentifiers) {
        return sampleIdentifiers.stream().anyMatch(this.samples::containsKey);
    }

    /**
     * Checks if the variant is associated with a specific feature.
     * <p>
     * This method determines whether the given feature identifier exists in the {@code features} map. The {@code features} map contains
     * associations between feature identifiers and their related alleles.
     *
     * @param featureIdentifier The unique identifier of the feature to check.
     * @return {@code true} if the feature identifier exists in the {@code features} map; {@code false} otherwise.
     */
    public boolean ofFeature(String featureIdentifier) {
        return this.features.containsKey(featureIdentifier);
    }

    /**
     * Checks if the variant is associated with any of the specified features.
     *
     * @param featureIdentifiers A collection of feature identifiers to check.
     * @return {@code true} if the variant is associated with at least one of the given feature identifiers, {@code false} otherwise.
     */
    public boolean ofFeatures(Collection<String> featureIdentifiers) {
        return featureIdentifiers.stream().anyMatch(this.features::containsKey);
    }

    /**
     * Checks if the variant is associated with a specific allele.
     *
     * @param alleleIdentifier The unique identifier of the allele to check.
     * @return {@code true} if the variant is associated with the given allele identifier, {@code false} otherwise.
     */
    public boolean ofAllele(String alleleIdentifier) {
        return this.features.values().stream().anyMatch(alleles -> alleles.contains(alleleIdentifier));
    }

    /**
     * Checks if the variant is associated with any of the specified alleles.
     *
     * @param alleleIdentifiers A collection of allele identifiers to check.
     * @return {@code true} if the variant is associated with at least one of the given allele identifiers, {@code false} otherwise.
     */
    public boolean ofAlleles(Collection<String> alleleIdentifiers) {
        return alleleIdentifiers.stream().anyMatch(this::ofAllele);
    }

    /**
     * Checks if the variant is filtered with respect to a specific sample.
     * <p>
     * Retrieves the variant call string associated with the given sample identifier from the `samples` map. It then checks if the variant
     * call is filtered using the {@link VariantCall#isFiltered(String)} method. If the sample identifier does not exist in the map, the
     * method returns {@code false}, indicating that the variant is considered not filtered.
     *
     * @param sampleIdentifier The identifier of the sample to check.
     * @return {@code true} if the variant is filtered for the specified sample; {@code false} otherwise.
     */
    public boolean isFiltered(String sampleIdentifier) {
        String call = this.samples.get(sampleIdentifier);
        return call != null && VariantCall.isFiltered(call);
    }

    /**
     * Checks if the variant is filtered with respect to all specified samples.
     * <p>
     * See {@link #isFiltered(String)}.
     *
     * @param sampleIdentifiers A collection of sample identifiers to check.
     * @return {@code true} if the variant is filtered for all specified samples; {@code false} otherwise.
     */
    public boolean isFiltered(Collection<String> sampleIdentifiers) {
        return sampleIdentifiers.stream().allMatch(this::isFiltered);
    }

    /**
     * Checks if all variant calls for this variant are filtered.
     * <p>
     * This method iterates through all variant calls in the `samples` map and checks if each call is filtered using the
     * {@link VariantCall#isFiltered(String)} method. If all calls are filtered, it returns {@code true}; otherwise, it returns
     * {@code false}.
     *
     * @return {@code true} if all variant calls for this variant are filtered; {@code false} otherwise.
     */
    public boolean isFiltered() {
        return this.samples.values().stream().allMatch(VariantCall::isFiltered);
    }

    /**
     * Retrieves a set of sample identifiers that have occurrences of this variant.
     * <p>
     * This method returns an unmodifiable set of sample identifiers associated with this variant. The collection ensures that external
     * modifications are not allowed, preserving the integrity of the data.
     *
     * @return A set of sample identifiers that have occurrences of this variant.
     */
    public Set<Tuple<String, String>> getRelatedSamples() {
        return this.samples.entrySet().stream().map(entry -> new Tuple<>(entry.getKey(), entry.getValue()))
                .collect(Collectors.toUnmodifiableSet());
    }

    /**
     * Retrieves the variant call string associated with a specific sample.
     * <p>
     * This method fetches the variant call string for the given sample identifier from the `samples` map. If the sample identifier does not
     * exist in the map, it returns {@code null}.
     *
     * @param sampleIdentifier The unique identifier of the sample whose relation is to be retrieved.
     * @return The variant call string associated with the given sample identifier, or {@code null} if the sample is not found.
     */
    public String getSampleRelation(String sampleIdentifier) {
        return this.samples.get(sampleIdentifier);
    }

    /**
     * Retrieves a set of allele identifiers associated with a specific feature for this variant.
     * <p>
     * This method fetches the set of allele identifiers for the given feature identifier from the `features` map. If the feature identifier
     * does not exist in the map, it returns an empty set.
     *
     * @param featureIdentifier The unique identifier of the feature whose allele relations are to be retrieved.
     * @return A set of allele identifiers associated with the given feature identifier, or an empty set if the feature is not found.
     */
    public Set<String> getFeatureRelation(String featureIdentifier) {
        return this.features.getOrDefault(featureIdentifier, Collections.emptySet());
    }

    /**
     * Retrieves a set of tuples representing the feature and allele occurrences associated with this variant.
     *
     * @return A set of tuples where each tuple contains:
     * <ul>
     *   <li>The feature identifier as a {@link String}.</li>
     *   <li>The allele identifier as a {@link String}.</li>
     * </ul>
     */
    public Set<Tuple<String, String>> getRelatedAlleles() {
        return this.features.entrySet().stream()
                .flatMap(entry -> entry.getValue().stream().map(allele -> new Tuple<>(entry.getKey(), allele))).collect(Collectors.toUnmodifiableSet());
    }

    /**
     * Associates a sample with this variant by adding the sample identifier and its associated variant calls.
     * <p>
     * This method updates the `samples` map by associating the given sample identifier with a string representation of the provided variant
     * calls. The variant calls are converted to strings using their `toString` method and concatenated with a pipe ('|') delimiter.
     * <p>
     * This will also mark the variant as active by setting the {@code active} property to {@code true}, indicating that the variant has
     * been modified. This will overwrite any existing association for the given sample identifier.
     *
     * @param sampleIdentifier The unique identifier of the sample to associate with this variant.
     * @param variantCalls     A set of {@link VariantCall} objects representing the variant calls to associate with the sample. Each
     *                         variant call is converted to its string representation.
     */
    public void addSampleRelation(String sampleIdentifier, Set<VariantCall> variantCalls) {
        // Convert the set of VariantCall objects to a single string, joined by the pipe ('|') character,
        // and associate it with the given sample identifier in the samples map.
        this.samples.put(sampleIdentifier, variantCalls.stream().map(VariantCall::toString).collect(Collectors.joining(Constants.PIPE)));
        this.active = true;
    }

    /**
     * Associates an allele and its parent feature with this variant.
     *
     * @param featureIdentifier The identifier of the feature to associate with this variant.
     * @param alleleIdentifier  The identifier of the allele to associate with the feature.
     */
    public void addAlleleRelation(String featureIdentifier, String alleleIdentifier) {
        this.features.putIfAbsent(featureIdentifier, new HashSet<>(32));
        this.features.get(featureIdentifier).add(alleleIdentifier);
    }

    /**
     * Removes the association between a sample and this variant.
     * <p>
     * This method removes the specified sample identifier from the `samples` map. After the removal, it checks if the `samples` map is
     * empty and returns the result.
     * <p>
     * This operation does not affect other associations or attributes of the variant.
     *
     * @param identifier The unique identifier of the sample to be disassociated from this variant.
     * @return {@code true} if the `samples` map is empty after the removal; {@code false} otherwise.
     */
    boolean removeSampleRelation(String identifier) {
        samples.remove(identifier);
        return samples.isEmpty();
    }

    /**
     * Removes the association between a specific allele and its parent feature for this variant.
     * <p>
     * This method removes the specified allele identifier from the set of alleles associated with the given feature identifier. If the set
     * of alleles for the feature becomes empty after the removal, the feature itself is removed from the `features` map.
     * <p>
     * In contrast to {@link #removeSampleRelation(String)}, no value is returned, as variants that are not associated with any features can
     * still be valid and exist in the storage.
     *
     * @param featureIdentifier The unique identifier of the feature to disassociate the allele from.
     * @param alleleIdentifier  The unique identifier of the allele to be removed from the feature.
     */
    void removeAlleleRelation(String featureIdentifier, String alleleIdentifier) {
        features.get(featureIdentifier).remove(alleleIdentifier);
        if (features.get(featureIdentifier).isEmpty()) features.remove(featureIdentifier);
    }

    /**
     * Converts this variant to a simplified stub representation.
     * <p>
     * This method creates a {@link Variant.Stub} object that encapsulates the position and alternative allele of this variant. The stub
     * serves as a simplified representation of the variant.
     *
     * @return A {@link Variant.Stub} object containing the position and alternative allele of this variant.
     */
    public Variant.Stub asStub() {
        return new Variant.Stub(this.position, this.alternative);
    }

    /**
     * Converts this variant to a masked stub representation.
     * <p>
     * This method creates a {@link Variant.Stub} object that encapsulates the position of this variant and a masked alternative allele
     * represented by {@link Constants#ANY_NUCLEOTIDE}.
     *
     * @return A {@link Variant.Stub} object containing the position of this variant and a masked alternative allele.
     */
    public Variant.Stub asMaskedStub() {
        return new Variant.Stub(this.position, Constants.ANY_NUCLEOTIDE);
    }

    /**
     * Converts the variant to its string representation.
     * <p>
     * This method generates a string representation of the variant based on its type:
     * <ul>
     *   <li>For {@link Type#SNV}, the format is: "prefix + position + reference + '>' + alternative".</li>
     *   <li>For {@link Type#INSERTION}, the format is: "prefix + position + '_' + (position + 1) + 'ins' + alternative".</li>
     *   <li>For {@link Type#DELETION}, the format is: "prefix + position + '_' + (position + reference.length() - 1) + 'del'".</li>
     * </ul>
     *
     * @return A {@link String} representing the variant.
     */
    public String toString() {
        StringBuilder sb = new StringBuilder();
        switch (type) {
            case SNV -> sb.append(position).append(reference).append(Constants.GREATER_THAN).append(alternative);
            case INSERTION -> sb.append(position).append(Constants.UNDER_SCORE).append(position + 1).append("ins").append(alternative);
            case DELETION -> sb.append(position).append(Constants.UNDER_SCORE).append(position + reference.length() - 1).append("del");
        }
        return sb.toString();
    }

    /**
     * Computes the hash code for this variant.
     * <p>
     * This method calculates the hash code of the variant based on its string representation.
     *
     * @return The hash code of the variant.
     */
    public int hashCode() {
        return this.toString().hashCode();
    }

    /**
     * Compares this variant to another object for equality.
     * <p>
     * This method checks if the provided object is the same instance as this object. If not, it verifies that the object is of the same
     * class and compares their string representations for equality.
     *
     * @param obj The object to compare with this {@link Variant} instance.
     * @return {@code true} if the objects are the same instance or if their string representations are equal; {@code false} otherwise.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        Variant that = (Variant) obj;
        return Objects.equals(this.toString(), that.toString());
    }

}