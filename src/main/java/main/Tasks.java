package main;

/**
 * {@link Enum} specifying tasks of MUSIAL to choose from.
 *
 * @author Simon Hackl
 */
public enum Tasks {
    /**
     * Task to build a MUSIAL storage file.
     */
    BUILD,
    /**
     * Task to view features from a MUSIAL storage file.
     */
    VIEW_FEATURES,
    /**
     * Task to view samples from a MUSIAL storage file.
     */
    VIEW_SAMPLES,
    /**
     * Task to view variants from a MUSIAL storage file.
     */
    VIEW_VARIANTS,
    /**
     * Task to export a table/matrix of variants wrt. one feature/gene.
     */
    EXPORT_TABLE,
    /**
     * Task to export a sequence in FASTA format wrt. one feature/gene.
     */
    EXPORT_SEQUENCE
}
