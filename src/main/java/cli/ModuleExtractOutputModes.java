package cli;

/**
 * {@link Enum} storing categories of output modes, i.e. output format, used for the {@link main.MusialModules#EXTRACT}
 * module.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public enum ModuleExtractOutputModes {
    /**
     * Output mode for writing plain unaligned sequences (.fasta).
     */
    SEQUENCE,
    /**
     * Output mode for writing aligned sequences (.fasta).
     */
    SEQUENCE_ALIGNED,
    /**
     * Output mode for writing tables (.tsv).
     */
    TABLE
}
