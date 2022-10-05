package datastructure;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;

/**
 * Represents a proteoform, i.e. a set of single-position amino-acid variants (may include inserted and deleted
 * positions), wrt. a {@link FeatureEntry}.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class ProteoformEntry {

    /**
     * The internal name to use for this proteoform. Should be derived with {@link ProteoformEntry#generateProteoformName(String)}
     */
    public final String name;
    /**
     * A {@link HashMap} yielding any meta-information {@link String} key-value pairs about this {@link ProteoformEntry}.
     */
    public final HashMap<String, String> annotations = new HashMap<>();
    /**
     * {@link TreeSet} of {@link String} values matching exactly one {@link SampleEntry#name} value. Represents the list of
     * all samples yielding this proteoform.
     */
    public final TreeSet<String> samples = new TreeSet<>();
    /**
     * {@link String} identifier used to store the wild type proteoform.
     */
    public final static String PROPERTY_NAME_REFERENCE_ID = "PF00000000WT";
    public final static String PROPERTY_NAME_VARIANTS = "VARIANTS";
    public final static String PROPERTY_NAME_DIVERGING_TERMINATION_POSITION = "DIV_TERM_POS";
    public final static String PROPERTY_NAME_DIVERGING_TERMINATION_TRUNCATED_PERCENTAGE = "DIV_TERM_TRC_PRC";
    public final static String PROPERTY_NAME_VARIABLE_POSITIONS_TOTAL = "VAR_POS";
    public final static String PROPERTY_NAME_CONGLOMERATION_INDEX = "CONG_IDX";

    /**
     * Constructor of {@link ProteoformEntry}.
     *
     * @param name           The internal name to use for the proteoform.
     * @param sampleId       The sample ID/name for which the proteoform was derived.
     * @param concatVariants Mandatory variants annotation to use for this proteoform.
     */
    public ProteoformEntry(String name, String sampleId, String concatVariants) {
        this.name = name;
        this.annotations.put(PROPERTY_NAME_VARIANTS, concatVariants);
        this.samples.add(sampleId);
    }

    /**
     * Returns a {@link String} that is used as internal proteoform name.
     * <p>
     * See {@link FeatureEntry#generateEntryName(String, String)} for more details.
     *
     * @param concatVariants {@link String} representation of the variants of one proteoform in the format <ALT>#<POS>.
     * @return {@link String} intended to be used as internal proteoform name.
     */
    public static String generateProteoformName(String concatVariants) {
        if (concatVariants.equals("")) {
            return ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
        } else {
            return FeatureEntry.generateEntryName(concatVariants, "PF");
        }
    }

}
