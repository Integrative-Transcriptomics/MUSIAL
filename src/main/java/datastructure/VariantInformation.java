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
     *
     * @param referenceContent The {@link String} to set as the reference content of this variant.
     */
    public VariantInformation(String referenceContent) {
        super();
        this.referenceContent = referenceContent;
    }

}