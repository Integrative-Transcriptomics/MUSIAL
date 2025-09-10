package main;

/**
 * {@link Enum} specifying MUSIAL tasks.
 */
public enum MusialTask {
    /**
     * Task to build a MUSIAL storage file.
     */
    BUILD,
    /**
     * Task to add (sample) data from VCF files to a MUSIAL storage file.
     */
    EXPAND,
    /**
     * Task to generate tables of various content from a MUSIAL storage file.
     */
    VIEW,
    /**
     * Task to export sequence data in FASTA format from a MUSIAL storage file.
     */
    SEQUENCE,
    /**
     * Task is undefined.
     */
    UNDEFINED
}
