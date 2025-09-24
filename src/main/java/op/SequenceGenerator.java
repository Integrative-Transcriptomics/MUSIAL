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

}
