package main;

/**
 * {@link Enum} specifying modules of MUSIAL to choose from.
 * <p>
 * Please add respective description to {@link Musial#printInfo()} documentation, if
 * new modules are added.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public enum MusialModules {
    /**
     * Generate a new variants dictionary JSON file.
     */
    BUILD,
    /**
     * Extracts nucleotide variants as a table (.tsv).
     */
    EXTRACT_NUC_VARIANTS_TABLE
}
