package datastructure;

/**
 * Container to store annotations and information of a single variant on a feature.
 *
 * @author Simon Hackl
 */
public class VariantInformation extends OccurrenceContainer {

    /**
     * The reference base content of this variant.
     */
    public final String referenceContent;

    /**
     * Constructor of {@link VariantInformation}.
     */
    public VariantInformation(String referenceContent) {
        super();
        this.referenceContent = referenceContent;
    }

}