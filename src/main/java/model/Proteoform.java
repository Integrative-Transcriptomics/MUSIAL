package model;

import util.Constants;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Represents a proteoform, which is a specific sequence type associated with a set of alleles.
 * <p>
 * This class extends {@link SequenceType} and provides functionality to manage relationships between the proteoform and associated
 * alleles.
 */
public class Proteoform extends SequenceType {

    /**
     * A set of allele identifiers associated with this proteoform.
     */
    private final Set<String> alleles = new HashSet<>();

    /**
     * Constructs a {@link Proteoform} instance with the specified variants.
     * <p>
     * This constructor initializes the proteoform by processing a list of {@link Variant.Stub} objects.
     *
     * @param variants A {@link List} of {@link Variant.Stub} objects representing the variants that define this proteoform.
     */
    public Proteoform(List<Variant.Stub> variants) {
        super(variants);
    }

    /**
     * Checks if this proteoform is associated with a specific allele by its identifier.
     *
     * @param identifier The unique identifier of the allele to check.
     * @return {@code true} if the allele is associated with this proteoform, {@code false} otherwise.
     */
    public boolean hasRelation(String identifier) {
        return this.alleles.contains(identifier);
    }

    /**
     * Retrieves a collection of allele identifiers associated with this proteoform.
     * <p>
     * This method returns an unmodifiable view of the set of allele identifiers associated with the proteoform. The unmodifiable set
     * ensures that the original set cannot be modified externally, preserving data integrity.
     *
     * @return A {@link Set} of allele identifiers related to this proteoform.
     */
    public Set<String> getRelatedAlleles() {
        return Collections.unmodifiableSet(this.alleles);
    }

    /**
     * Associates an allele with this proteoform.
     *
     * @param alleleIdentifier The unique identifier of the allele to associate with this proteoform.
     */
    public void addRelation(String alleleIdentifier) {
        this.alleles.add(alleleIdentifier);
    }

    /**
     * Determines if the proteoform is disrupted based on its sequence ontology (SO) effects.
     *
     * @return {@code true} if the proteoform is disrupted (i.e., has "start_lost" or "stop_gained" SO effects), {@code false} otherwise.
     */
    public boolean isDisrupted() {
        Set<String> soEffects = getAttributeSet(Constants.AttributesKeys.SO_EFFECTS);
        return soEffects.contains("start_lost") || soEffects.contains("stop_gained");
    }

    /**
     * Converts the proteoform to its string representation.
     * <p>
     * This method returns the unique identifier of the proteoform as its string representation.
     *
     * @return A {@link String} representing the unique identifier of the proteoform.
     */
    public String toString() {
        return this._id;
    }

    /**
     * Computes the hash code for this proteoform.
     * <p>
     * This method calculates the hash code of the proteoform based on its string representation.
     *
     * @return The hash code of the proteoform.
     */
    public int hashCode() {
        return this.toString().hashCode();
    }

    /**
     * Compares this proteoform to another object for equality.
     * <p>
     * This method checks if the provided object is the same instance as this proteoform. If not, it verifies that the object is of the same
     * class and compares their string representations for equality.
     *
     * @param obj The object to compare with this {@link Proteoform} instance.
     * @return {@code true} if the objects are the same instance or if their string representations are equal; {@code false} otherwise.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        Proteoform that = (Proteoform) obj;
        return this.toString().equals(that.toString());
    }

}
