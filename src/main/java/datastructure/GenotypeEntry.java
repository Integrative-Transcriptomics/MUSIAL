package datastructure;

import java.util.HashMap;
import java.util.TreeSet;

/**
 * Represents a proteoform, i.e. a set of single-position amino-acid variants (may include inserted and deleted
 * positions), wrt. a {@link FeatureEntry}.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class GenotypeEntry {

    /**
     * The internal name to use for this proteoform. Should be derived with {@link GenotypeEntry#generateGenotypeName(String)}
     */
    public final String name;
    /**
     * A {@link HashMap} yielding any meta-information {@link String} key-value pairs about this {@link GenotypeEntry}.
     */
    public final HashMap<String, String> annotations = new HashMap<>();
    /**
     * {@link TreeSet} of {@link String} values matching exactly one {@link SampleEntry#name} value. Represents the list of
     * all samples yielding this proteoform.
     */
    public final TreeSet<String> samples = new TreeSet<>();

    public final static String PROPERTY_NAME_REFERENCE_ID = "GT00000000WT";

    public final static String PROPERTY_NAME_VARIANTS = "VARIANTS";

    /**
     * Constructor of {@link GenotypeEntry}.
     *
     * @param name           The internal name to use for the proteoform.
     * @param sId            The sample ID/name for which the proteoform was derived.
     * @param concatVariants Mandatory variants annotation to use for this proteoform.
     */
    public GenotypeEntry(String name, String sId, String concatVariants) {
        this.name = name;
        this.annotations.put(PROPERTY_NAME_VARIANTS, concatVariants);
        this.samples.add(sId);
    }

    /**
     * Returns a {@link String} that is used as internal genotype name.
     * <p>
     * See {@link FeatureEntry#generateEntryName(String, String)} for more details.
     *
     * @param concatVariants {@link String} representation of the variants of one genotype in the format <ALT>#<POS>.
     * @return {@link String} intended to be used as internal genotype name.
     */
    public static String generateGenotypeName(String concatVariants) {
        if (concatVariants.equals("")) {
            return GenotypeEntry.PROPERTY_NAME_REFERENCE_ID;
        } else {
            return FeatureEntry.generateEntryName(concatVariants, "GT");
        }
    }

}
