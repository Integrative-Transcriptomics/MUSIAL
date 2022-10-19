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
     * Updates an existing or generates a new variants dictionary JSON file.
     */
    BUILD,
    /**
     * Infer nucleotide and/or amino-acid sequences from an existing variants dictionary JSON file.
     */
    inferSequences,
    /**
     * Run SnpEff to annotate nucleotide variants and compute various statistics.
     */
    statistics
}
