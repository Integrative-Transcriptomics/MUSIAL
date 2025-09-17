package model;

import org.ehcache.spi.serialization.Serializer;
import org.ehcache.spi.serialization.SerializerException;
import util.Constants;

import java.nio.ByteBuffer;
import java.util.ArrayList;
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
         * Converts the alternative allele to a string representation, i.e., the reference and alternative allele, and allelic depth
         * separated by {@link Constants#COLON}.
         *
         * @return A {@link String} representing the alternative allele with allelic depth.
         */
        private String asString() {
            return reference + Constants.COLON + alternative + Constants.COLON + allelicDepth;
        }

        /**
         * Converts the alternative allele to its string representation.
         * <p>
         * This method generates a string representation of the alternative allele by concatenating the reference and alternative allele
         * separated by {@link Constants#GREATER_THAN}.
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
     * Checks if the variant call is buried.
     * <p>
     * A variant call is considered "buried" if it is not filtered, but the called alternative is a missing allele due to an upstream
     * deletion (*).
     *
     * @return {@code true} if the variant call is buried; {@code false} otherwise.
     */
    public boolean isBuried() {
        return !isFiltered() && getCalledAlternative().equals("*");
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
     * Serializer for the {@link VariantCall} class.
     * <p>
     * This class implements the {@link Serializer} interface to provide custom serialization and deserialization logic for
     * {@link VariantCall} objects. It handles the conversion of {@link VariantCall} instances to and from {@link ByteBuffer} format,
     * allowing for efficient storage and retrieval in caching systems.
     */
    public static class VariantCallSerializer implements Serializer<VariantCall> {

        /**
         * Constructs a new {@link VariantCallSerializer} instance.
         *
         * @param object The instance to serialize.
         * @return A {@link ByteBuffer} containing the serialized data of the {@link VariantCall} object.
         * @throws SerializerException If an error occurs during serialization.
         */
        @Override
        public ByteBuffer serialize(VariantCall object) throws SerializerException {
            try {
                // Estimate the buffer size
                int size = Short.BYTES + Float.BYTES + Integer.BYTES; // totalDepth, callEntropy, alternatives size
                size += object.flag().name().getBytes().length + Integer.BYTES; // flag string length
                for (VariantCall.CallAlternative alternative : object.alternatives()) {
                    size += alternative.reference().getBytes().length + Integer.BYTES; // reference string length
                    size += alternative.alternative().getBytes().length + Integer.BYTES; // alternative string length
                    size += Short.BYTES; // allelicDepth
                }

                ByteBuffer buffer = ByteBuffer.allocate(size);

                // Serialize the flag
                byte[] flagBytes = object.flag().name().getBytes();
                buffer.putInt(flagBytes.length);
                buffer.put(flagBytes);

                // Serialize the totalDepth
                buffer.putShort(object.totalDepth());

                // Serialize the callEntropy
                buffer.putFloat(object.callEntropy());

                // Serialize the alternatives
                buffer.putInt(object.alternatives().size());
                for (VariantCall.CallAlternative alternative : object.alternatives()) {
                    byte[] referenceBytes = alternative.reference().getBytes();
                    buffer.putInt(referenceBytes.length);
                    buffer.put(referenceBytes);

                    byte[] alternativeBytes = alternative.alternative().getBytes();
                    buffer.putInt(alternativeBytes.length);
                    buffer.put(alternativeBytes);

                    buffer.putShort(alternative.allelicDepth());
                }

                buffer.flip();
                return buffer;
            } catch (Exception e) {
                throw new SerializerException("Error during serialization", e);
            }
        }

        /**
         * Deserializes a {@link VariantCall} object from a {@link ByteBuffer}.
         *
         * @param binary The binary representation to deserialize.
         * @return The deserialized {@link VariantCall} object.
         * @throws SerializerException If an error occurs during deserialization.
         */
        @Override
        public VariantCall read(ByteBuffer binary) throws SerializerException {
            try {
                // Read the flag
                int flagLength = binary.getInt();
                byte[] flagBytes = new byte[flagLength];
                binary.get(flagBytes);
                VariantCall.Flag flag = VariantCall.Flag.valueOf(new String(flagBytes));

                // Read the totalDepth
                short totalDepth = binary.getShort();

                // Read the callEntropy
                float callEntropy = binary.getFloat();

                // Read the alternatives
                int alternativesSize = binary.getInt();
                List<VariantCall.CallAlternative> alternatives = new ArrayList<>(alternativesSize);
                for (int i = 0; i < alternativesSize; i++) {
                    // Read reference
                    int referenceLength = binary.getInt();
                    byte[] referenceBytes = new byte[referenceLength];
                    binary.get(referenceBytes);
                    String reference = new String(referenceBytes);

                    // Read alternative
                    int alternativeLength = binary.getInt();
                    byte[] alternativeBytes = new byte[alternativeLength];
                    binary.get(alternativeBytes);
                    String alternative = new String(alternativeBytes);

                    // Read allelicDepth
                    short allelicDepth = binary.getShort();

                    // Create CallAlternative and add to list
                    alternatives.add(new VariantCall.CallAlternative(reference, alternative, allelicDepth));
                }

                // Construct and return the VariantCall instance
                return new VariantCall(flag, totalDepth, callEntropy, alternatives);
            } catch (Exception e) {
                throw new SerializerException("Error during deserialization", e);
            }
        }

        /**
         * Compares a {@link VariantCall} object with its binary representation for equality.
         *
         * @param object The instance to check for equality.
         * @param binary The serialized form to check against.
         * @return {@code true} if the object and its binary representation are equal; {@code false} otherwise.
         * @throws SerializerException If an error occurs during comparison.
         */
        @Override
        public boolean equals(VariantCall object, ByteBuffer binary) throws SerializerException {
            try {
                // Deserialize the ByteBuffer into a VariantCall object
                VariantCall deserializedObject = read(binary);

                // Compare the provided object with the deserialized object
                return object.equals(deserializedObject);
            } catch (Exception e) {
                throw new SerializerException("Error during comparison", e);
            }
        }

    }

    public static VariantCall fromString(String s) {
        String[] parts = s.split(Constants.SEMICOLON);
        Flag flag = Flag.valueOf(parts[0].toUpperCase());
        short totalDepth = Short.parseShort(parts[1]);
        float callEntropy = Float.parseFloat(parts[2]);
        List<CallAlternative> alternatives = new ArrayList<>();
        if (parts.length > 3) {
            String[] alts = parts[3].split(Constants.COMMA);
            for (String alt : alts) {
                String[] altParts = alt.split(Constants.COLON);
                String reference = altParts[0];
                String alternative = altParts[1];
                short allelicDepth = Short.parseShort(altParts[2]);
                alternatives.add(new CallAlternative(reference, alternative, allelicDepth));
            }
        }
        return new VariantCall(flag, totalDepth, callEntropy, alternatives);
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
            sb.append(alternatives.get(i).asString());
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