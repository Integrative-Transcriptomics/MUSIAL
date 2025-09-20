package model;

import htsjdk.samtools.util.Tuple;
import util.Constants;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Represents a single biological sample.
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
     * Indicates whether the sample is active.
     * <p>
     * This boolean flag is used to mark samples that are new or were modified. It is set to {@code true} for newly created or touched
     * entries, but will not be serialized. During deserialization, it is assumed that all samples are idle and should be set to
     * {@code false} (see {@link Storage#typeAdapter()}).
     */
    transient boolean active;

    /**
     * Constructs a new {@link Sample} instance with the specified id and initial capacity for the alleles map.
     * <p>
     * This constructor initializes a {@link Sample} object with the given identifier and allocates a {@link HashMap} instance for the
     * {@link #alleles} field with the specified initial capacity. The {@link #_id} field is set to the provided identifier, and the
     * superclass constructor is invoked to initialize inherited properties.
     * <p>
     * Instances of this class are not related to variants from within this constructor. Instead, variants are linked via the
     * {@link Variant#samples} and {@link Allele} classes.
     *
     * @param identifier       The unique identifier of the sample, used as its unique key.
     * @param capacityFeatures The initial capacity for the {@link #alleles} map, which stores feature-allele associations.
     */
    public Sample(String identifier, int capacityFeatures) {
        super(); // Call the constructor of the superclass to initialize inherited properties.
        this._id = identifier; // Assign the unique identifier to the _id field.
        this.alleles = new HashMap<>(capacityFeatures); // Initialize the alleles map with the specified capacity.
        this.active = true; // Mark the sample as active upon creation.
    }

    /**
     * Associates a specific allele of a feature with this sample.
     * <p>
     * This method updates the {@link #alleles} map by setting the sequence type (allele) for the specified feature. The feature and allele
     * are identified by their id.
     *
     * @param featureIdentifier The id of the feature ({@link Feature#_id}) to associate with the allele.
     * @param alleleIdentifier  The unique identifier of the allele ({@link SequenceType#_id}) to set for the feature.
     */
    public void addRelation(String featureIdentifier, String alleleIdentifier) {
        this.alleles.put(featureIdentifier, alleleIdentifier);
    }

    /**
     * Retrieves the allele associated with a specific feature in this sample.
     * <p>
     * This method looks up the {@link #alleles} map to find the allele associated with the given feature identifier. If no association
     * exists, it returns the default reference allele defined in {@link Constants#REFERENCE}.
     *
     * @param featureIdentifier The unique identifier of the feature ({@link Feature#name}) to retrieve the associated allele for.
     * @return The unique identifier of the allele ({@link SequenceType#_id}) associated with the feature, or the default reference allele
     * if no association exists.
     */
    public String getRelatedAllele(String featureIdentifier) {
        return this.alleles.getOrDefault(featureIdentifier, Constants.REFERENCE);
    }

    /**
     * Retrieves all feature-allele associations in this sample.
     * <p>
     * This method converts the {@link #alleles} map, which stores feature names as keys and their associated allele identifiers as values,
     * into a collection of {@link Tuple} objects. Each tuple contains a feature name and its corresponding allele identifier.
     *
     * @return A {@link Set} of {@link Tuple} objects, where each tuple represents a feature-allele association.
     */
    public Set<Tuple<String, String>> getRelatedAlleles() {
        return this.alleles.entrySet().stream()
                .map(entry -> new Tuple<>(entry.getKey(), entry.getValue()))
                .collect(Collectors.toUnmodifiableSet());
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
     * Converts the sample to its string representation.
     * <p>
     * This method returns the unique identifier of the sample as its string representation.
     *
     * @return A {@link String} representing the unique identifier of the sample.
     */
    public String toString() {
        return this._id;
    }

    /**
     * Computes the hash code for this sample.
     * <p>
     * This method calculates the hash code of the sample based on its unique identifier.
     *
     * @return The hash code of the sample.
     */
    public int hashCode() {
        return this._id.hashCode();
    }

    /**
     * Compares this sample to another object for equality.
     * <p>
     * This method checks if the provided object is the same instance as this sample. If not, it verifies that the object is of the same
     * class and compares their unique identifiers for equality.
     *
     * @param obj The object to compare with this {@link Sample} instance.
     * @return {@code true} if the objects are the same instance or if their unique identifiers are equal; {@code false} otherwise.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        Sample that = (Sample) obj;
        return this._id.equals(that._id);
    }

}
