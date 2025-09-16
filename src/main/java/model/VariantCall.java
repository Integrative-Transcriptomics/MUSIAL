package model;

import util.Constants;

import java.util.List;

/**
 * Represents a variant call in a biological sample.
 * <p>
 * This record encapsulates the details of a variant call, including:
 * <ul>
 *   <li>A flag indicating the state of this call (e.g., filtered due to low frequency, reference call).</li>
 *   <li>The total read depth at the variant site.</li>
 *   <li>The normalized entropy of the call, representing the uncertainty in the variant call.</li>
 *   <li>A list of alternative alleles associated with the variant call, represented as {@link CallAlternative} objects.</li>
 * </ul>
 * <p>
 * Variant calls are stored in the {@link Sample#variantCalls} property of the model.
 *
 * @param flag         Flag indicating the state of this variant call wrt. filters.
 * @param totalDepth   The total read depth at the variant site.
 * @param callEntropy  The normalized entropy of the call.
 * @param alternatives A list of alternative alleles associated with the variant call.
 */
public record VariantCall(Flag flag, short totalDepth, float callEntropy, List<CallAlternative> alternatives) {

    /**
     * Enumeration of flags representing the status of a variant call.
     * <p>
     * This enum defines various flags that describe the state or quality of a variant call:
     * <ul>
     *   <li>{@link #LOW_FREQUENCY} - Indicates that the variant call has a frequency below the minimum threshold.</li>
     *   <li>{@link #LOW_COVERAGE} - Indicates that the variant call has a coverage below the minimum threshold.</li>
     *   <li>{@link #MISSING_UPSTREAM_DELETION} - Indicates that the variant call is missing due to an upstream deletion.</li>
     *   <li>{@link #REFERENCE_CALL} - Indicates that the variant call corresponds to the reference allele.</li>
     *   <li>{@link #PASS} - Indicates that the variant call has passed all filters.</li>
     * </ul>
     */
    public enum Flag {
        /**
         * Indicates that the variant call has a frequency below the minimum threshold.
         */
        LOW_FREQUENCY,
        /**
         * Indicates that the variant call has a coverage below the minimum threshold.
         */
        LOW_COVERAGE,
        /**
         * Indicates that the variant call is declared missing due to an upstream deletion, but the upstream deletion is not detected.
         */
        MISSING_UPSTREAM_DELETION,
        /**
         * Indicates that the variant call corresponds to the reference allele.
         */
        REFERENCE_CALL,
        /**
         * Indicates that the variant call has passed all filters.
         */
        PASS
    }

    /**
     * Represents one alternative in the context of a variant call.
     * <p>
     * This record encapsulates the details of an alternative allele, including:
     * <ul>
     *   <li>The reference nucleotide sequence as a {@link String}.</li>
     *   <li>The alternative nucleotide sequence as a {@link String}.</li>
     *   <li>The allelic depth (number of reads supporting the alternative allele) as an {@link Integer}.</li>
     * </ul>
     * <p>
     * Call alternatives are stored in the {@link VariantCall#alternatives} property of the model.
     *
     * @param reference    The reference allele.
     * @param alternative  The alternative allele.
     * @param allelicDepth The number of reads supporting the alternative allele.
     */
    public record CallAlternative(String reference, String alternative, short allelicDepth) {

        /**
         * Converts the alternative allele to its string representation.
         * <p>
         * This method generates a string representation of the alternative allele by concatenating the reference allele, the alternative
         * allele, and the allelic depth, separated by colons.
         *
         * @return A {@link String} representing the alternative allele.
         */
        public String toString() {
            return reference + Constants.GREATER_THAN + alternative;
        }

        /**
         * Computes the hash code for this alternative allele.
         * <p>
         * This method calculates the hash code of the alternative allele based on its string representation.
         *
         * @return The hash code of the alternative allele.
         */
        public int hashCode() {
            return this.toString().hashCode();
        }

        /**
         * Compares this alternative allele to another object for equality.
         * <p>
         * This method checks if the provided object is the same instance as this object. If not, it verifies that the object is of the same
         * class and compares their string representations for equality.
         *
         * @param obj The object to compare with this {@link CallAlternative} instance.
         * @return {@code true} if the objects are the same instance or if their string representations are equal; {@code false} otherwise.
         */
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (obj == null || getClass() != obj.getClass()) return false;
            CallAlternative that = (CallAlternative) obj;
            return this.toString().equals(that.toString());
        }

    }

    /**
     * Checks if the variant call is filtered.
     * <p>
     * This method determines whether the variant call has been flagged as filtered based on its status. A variant call is considered
     * filtered if its flag is set to one of the following:
     * <ul>
     *   <li>{@link Flag#LOW_FREQUENCY} - The variant call has a frequency below the minimum threshold.</li>
     *   <li>{@link Flag#LOW_COVERAGE} - The variant call has a coverage below the minimum threshold.</li>
     *   <li>{@link Flag#MISSING_UPSTREAM_DELETION} - The variant call is missing due to an upstream deletion.</li>
     * </ul>
     *
     * @return {@code true} if the variant call is filtered; {@code false} otherwise.
     */
    public boolean isFiltered() {
        return (flag.equals(Flag.LOW_FREQUENCY) || flag.equals(Flag.LOW_COVERAGE) || flag.equals(Flag.MISSING_UPSTREAM_DELETION));
    }

    /**
     * Checks if the variant call corresponds to the reference allele.
     * <p>
     * This method determines whether the variant call's flag is set to {@link Flag#REFERENCE_CALL}, indicating that the variant call
     * matches the reference allele.
     *
     * @return {@code true} if the variant call corresponds to the reference allele; {@code false} otherwise.
     */
    public boolean isReference() {
        return flag.equals(Flag.REFERENCE_CALL);
    }

    /**
     * Retrieves the reference allele from the first alternative in the list.
     * <p>
     * This method assumes that the list of alternatives is not empty and returns the reference allele of the first {@link CallAlternative}
     * object in the list.
     *
     * @return The reference allele as a {@link String}.
     */
    public String getCalledReference() {
        return alternatives.get(0).reference;
    }

    /**
     * Retrieves the alternative allele from the first alternative in the list.
     * <p>
     * This method assumes that the list of alternatives is not empty and returns the alternative allele of the first
     * {@link CallAlternative} object in the list.
     *
     * @return The alternative allele as a {@link String}.
     */
    public String getCalledAlternative() {
        return alternatives.get(0).alternative;
    }

    /**
     * Converts the variant call to its string representation.
     * <p>
     * This method generates a string representation of the variant call, including its flag, total depth, entropy, and a formatted list of
     * alternative alleles. Each alternative allele is represented by its string representation, and alternatives are separated by commas.
     *
     * @return A {@link String} representing the variant call.
     */
    public String toString() {
        StringBuilder sb =
                new StringBuilder(flag.name().toLowerCase()).append(Constants.SEMICOLON).append(totalDepth).append(Constants.SEMICOLON)
                        .append(String.format("%.3f", callEntropy)).append(Constants.SEMICOLON);
        for (int i = 0; i < alternatives.size(); i++) {
            sb.append(alternatives.get(i).toString());
            if (i < alternatives.size() - 1) sb.append(Constants.COMMA);
        }
        return sb.toString();
    }

    /**
     * Computes the hash code for this variant call.
     * <p>
     * This method calculates the hash code of the variant call based on its string representation.
     *
     * @return The hash code of the variant call.
     */
    public int hashCode() {
        return this.toString().hashCode();
    }

    /**
     * Compares this variant call to another object for equality.
     * <p>
     * This method checks if the provided object is the same instance as this object. If not, it verifies that the object is of the same
     * class and compares their string representations for equality.
     *
     * @param obj The object to compare with this {@link VariantCall} instance.
     * @return {@code true} if the objects are the same instance or if their string representations are equal; {@code false} otherwise.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        VariantCall that = (VariantCall) obj;
        return this.toString().equals(that.toString());
    }

}