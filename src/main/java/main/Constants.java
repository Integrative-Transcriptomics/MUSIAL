package main;

/**
 * {@link Enum} to store constant project-wide values, such as names of MUSIAL modules or attribute values.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public enum Constants {
    /**
     * MUSIAL module name; Updates an existing or generates a new variants dictionary JSON file.
     */
    updateVDict,
    /**
     * MUSIAL module name; Infer nucleotide and/or amino-acid sequences from an existing variants dictionary JSON file.
     */
    inferSequences,
    /**
     * MUSIAL module name; Run SnpEff to annotate nucleotide variants and compute various statistics.
     */
    statistics
}
