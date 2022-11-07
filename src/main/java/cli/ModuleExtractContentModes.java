package cli;

/**
 * {@link Enum} storing categories of content modes, i.e. aminoacid or nucleotide variants, used for the
 * {@link main.MusialModules#EXTRACT} module.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public enum ModuleExtractContentModes {
    /**
     * Content mode for writing nucleotide variants.
     */
    NUCLEOTIDE,
    /**
     * Content mode for writing aminoacid variants.
     */
    AMINOACID
}
