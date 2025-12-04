package op;

import exceptions.MusialException;

/**
 * Interface for generating sequences based on a sample identifier.
 */
public interface SequenceGenerator {

    /**
     * Generates a sequence based on the provided sample identifier.
     *
     * @param sampleIdentifier the identifier for the sample
     * @return the generated sequence as a String
     * @throws MusialException If an error occurs during sequence generation
     */
    String getSequence(String sampleIdentifier) throws MusialException;

    /**
     * Retrieves the name of the sequence generator.
     *
     * @param forFile true to return a file-friendly name, false for the standard name.
     * @return the name of the sequence generator.
     */
    String getName(boolean forFile);

    /**
     * Returns the length of the sequences that can be generated.
     *
     * @return the length of the generated sequence as an integer.
     */
    int getSize();

    /**
     * Indicates whether this generator has a specific feature.
     *
     * @return true if the generator has the feature, false otherwise.
     */
    boolean hasFeature();

}
