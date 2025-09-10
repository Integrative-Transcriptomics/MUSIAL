package model;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


/**
 * Represents an allele, which is a specific sequence type associated with a set of samples and a proteoform.
 * <p>
 * This class extends {@link SequenceType} and provides functionality to manage relationships between the allele and associated samples, as
 * well as the proteoform identifier.
 */
public class Allele extends SequenceType {

    /**
     * The unique identifier of the proteoform associated with this allele.
     */
    private String proteoform;

    /**
     * A set of sample identifiers associated with this allele.
     */
    private final Set<String> samples = new HashSet<>();

    /**
     * Constructs an {@link Allele} instance with the specified variants.
     *
     * @param variants A list of {@link Variant.Stub} objects representing the variants that define this allele.
     */
    public Allele(List<Variant.Stub> variants) {
        super(variants);
    }

    /**
     * Checks if this allele is associated with a specific sample by its identifier.
     *
     * @param identifier The unique identifier of the sample to check.
     * @return {@code true} if the sample is associated with this allele, {@code false} otherwise.
     */
    public boolean hasRelation(String identifier) {
        return this.samples.contains(identifier);
    }

    /**
     * Retrieves a collection of sample identifiers associated with this allele.
     * <p>
     * This method returns an unmodifiable view of the set of sample identifiers associated with the allele. The unmodifiable set ensures
     * that the original set cannot be modified externally, preserving data integrity.
     *
     * @return A {@link Set} of sample identifiers related to this allele.
     */
    public Set<String> getRelatedSamples() {
        return Collections.unmodifiableSet(this.samples);
    }

    /**
     * Retrieves the number of unique samples associated with this allele.
     *
     * @return The count of sample identifiers related to this allele.
     */
    public int getRelatedSamplesCount() {
        return this.samples.size();
    }

    /**
     * Retrieves the proteoform identifier associated with this allele.
     *
     * @return The proteoform identifier as a {@link String}, or {@code null} if no proteoform is associated.
     */
    public String getProteoform() {
        return this.proteoform;
    }

    /**
     * Associates a sample with this allele.
     *
     * @param sampleIdentifier The unique identifier of the sample to associate with this allele.
     */
    public void addRelation(String sampleIdentifier) {
        this.samples.add(sampleIdentifier);
    }

    /**
     * Sets the proteoform identifier for this allele.
     *
     * @param identifier The unique identifier of the proteoform to associate with this allele.
     */
    public void setProteoform(String identifier) {
        this.proteoform = identifier;
    }

    /**
     * Converts the allele to its string representation.
     * <p>
     * This method returns the unique identifier of the allele as its string representation.
     *
     * @return A {@link String} representing the unique identifier of the allele.
     */
    public String toString() {
        return this._id;
    }

    /**
     * Computes the hash code for this allele.
     * <p>
     * This method calculates the hash code of the allele based on its unique identifier.
     *
     * @return The hash code of the allele.
     */
    public int hashCode() {
        return this._id.hashCode();
    }

    /**
     * Compares this allele to another object for equality.
     * <p>
     * This method checks if the provided object is the same instance as this allele. If not, it verifies that the object is of the same
     * class and compares their unique identifiers for equality.
     *
     * @param obj The object to compare with this {@link Allele} instance.
     * @return {@code true} if the objects are the same instance or if their unique identifiers are equal; {@code false} otherwise.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        Allele that = (Allele) obj;
        return this._id.equals(that._id);
    }

}
