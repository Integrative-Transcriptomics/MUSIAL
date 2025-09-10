package model;

import util.Bio;
import util.Constants;
import util.IO;

import java.util.*;

/**
 * Represents a sequence type with associated variants and attributes.
 * <p>
 * It extends the {@link Attributes} class to inherit functionality for managing attributes associated with the sequence type. This class is
 * extended by the {@link Allele} and {@link Proteoform} classes.
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
    private final NavigableMap<Integer, String> variants = new TreeMap<>();

    /**
     * Constructs a new {@link SequenceType} instance with the specified variants.
     * <p>
     * This constructor initializes the sequence type by processing a list of {@link Variant.Stub} objects. It calculates the sequence
     * length deviation, generates a unique identifier for the sequence type, and populates the {@link #variants} map with the provided
     * variants.
     *
     * @param variants A {@link List} of {@link Variant.Stub} objects representing the variants that define this sequence type.
     */
    SequenceType(List<Variant.Stub> variants) {
        super();
        int lengthDelta = 0;

        // Builds a unique identifier for the sequence type based on the variants.
        StringBuilder identifierBuilder = new StringBuilder(variants.size() * String.valueOf(variants.get(0).position()).length());

        // Sorts the variants by their position in ascending order.
        variants.sort(Comparator.comparingInt(Variant.Stub::position));

        // Processes each variant to populate the variants map and calculate the length deviation.
        for (Variant.Stub stub : variants) {
            this.variants.put(stub.position(), stub.alternative());
            identifierBuilder.append(stub.position()).append(stub.alternative());
            int variantLength = stub.alternative().length();

            // Updates the length deviation based on the type of variant (insertion or deletion).
            lengthDelta += Bio.isInsertion(stub.alternative()) ? variantLength - 1 :
                    Bio.isDeletion(stub.alternative()) ? -(variantLength - 1) : 0;
        }

        // Appends each variant's position and alternative allele to the identifier builder.
        variants.forEach(variant -> {
            this.variants.put(variant.position(), variant.alternative());
            identifierBuilder.append(variant.position()).append(variant.alternative());
        });

        // Generates a unique identifier for the sequence type using an MD5 hash of the identifier string.
        this._id = IO.md5Hash(identifierBuilder.toString());

        // Adds an attribute for the sequence length deviation.
        addAttribute(Constants.AttributesKeys.SEQUENCE_LENGTH_DEVIATION, String.valueOf(lengthDelta));
    }

    /**
     * Checks if a specific variant exists at the given position.
     * <p>
     * This method verifies whether a variant with the specified position and alternative allele is present in the {@link #variants} map. It
     * first checks if the position exists as a key in the map and then compares the associated value (alternative allele) with the provided
     * alternative allele.
     *
     * @param position    The position of the variant to check.
     * @param alternative The alternative allele to check for at the specified position.
     * @return {@code true} if the variant exists at the given position with the specified alternative allele; {@code false} otherwise.
     */
    public boolean hasVariant(int position, String alternative) {
        return this.variants.containsKey(position) && this.variants.get(position).equals(alternative);
    }

    /**
     * Checks if this sequence type has a variant at the specified position.
     *
     * @param position The position to check for a variant.
     * @return {@code true} if a variant exists at the specified position, {@code false} otherwise.
     */
    public boolean hasVariantAt(int position) {
        return this.variants.containsKey(position);
    }

    /**
     * Retrieves the alternative allele at the specified position.
     * <p>
     * This method fetches the alternative allele for a given position from the {@link #variants} map. If no variant exists at the specified
     * position, it returns {@code null}.
     *
     * @param position The position of the variant to retrieve.
     * @return The alternative allele as a {@link String}, or {@code null} if no variant exists at the position.
     */
    public String getVariant(int position) {
        return this.variants.getOrDefault(position, null);
    }

    /**
     * Retrieves a variant object for a specific position and contig.
     * <p>
     * This method fetches the alternative allele at the specified position from the {@link #variants} map. If the alternative allele
     * exists, it retrieves the corresponding {@link Variant} object from the provided {@link Contig}. If no alternative allele exists at
     * the position, it returns {@code null}.
     *
     * @param position The position of the variant to retrieve.
     * @param contig   The {@link Contig} object to retrieve the variant from.
     * @return The {@link Variant} object, or {@code null} if no variant exists at the position.
     */
    Variant getVariant(int position, Contig contig) {
        String alternative = this.variants.getOrDefault(position, null);
        if (alternative == null) return null;
        return contig.getVariant(position, alternative);
    }

    /**
     * Retrieves all variants as stubs.
     * <p>
     * Converts the {@link #variants} map into a list of {@link Variant.Stub} objects, where each stub contains the position and alternative
     * allele. Returns an empty list if no variants exist. The resulting list is unmodifiable.
     *
     * @return An unmodifiable {@link List} of {@link Variant.Stub} objects.
     */
    public List<Variant.Stub> getVariants() {
        if (this.variants.isEmpty()) {
            return Collections.emptyList();
        }

        // Convert map entries to Variant.Stub objects and collect into a list
        return Collections.unmodifiableList(
                this.variants.entrySet().stream()
                        .map(entry -> new Variant.Stub(entry.getKey(), entry.getValue()))
                        .toList()
        );
    }

    /**
     * Retrieves all variants as full objects for a specific contig.
     * <p>
     * Converts the {@link #variants} map into a list of {@link Variant} objects by fetching the corresponding variant from the provided
     * {@link Contig}. Returns an unmodifiable list to ensure immutability.
     *
     * @param contig The {@link Contig} object to retrieve the variants from.
     * @return A {@link List} of {@link Variant} objects representing all variants for the contig.
     */
    List<Variant> getVariants(Contig contig) {
        if (this.variants.isEmpty()) {
            return Collections.emptyList();
        }

        // Map entries to Variant objects, filter non-null, and collect into a list
        return Collections.unmodifiableList(
                this.variants.entrySet().stream()
                        .map(entry -> contig.getVariant(entry.getKey(), entry.getValue()))
                        .filter(Objects::nonNull)
                        .toList()
        );
    }

}