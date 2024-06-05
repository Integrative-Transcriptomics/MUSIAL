package datastructure;

/**
 * Container to store a biological variation of a feature, i.e., an allele of a gene or proteoform of a protein.
 *
 * @author Simon Hackl
 */
public class Form extends OccurrenceContainer {

    /**
     * The internal name to use for this object.
     */
    public final String name;

    /**
     * Variants that occur for this form in the format {@code POS:ALT;...}.
     */
    public final String variants;

    /**
     * Constructor of {@link Form}.
     *
     * @param name     Name to use for this object.
     * @param variants Variants fingerprint for this object, i.e., sequence of {@code POS:ALT;...}.
     */
    public Form(String name, String variants) {
        super();
        this.name = name;
        this.variants = variants;
    }

}
