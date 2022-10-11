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
public class AlleleEntry {

    /**
     * The internal name to use for this allele. Should be derived with {@link AlleleEntry#generateAlleleName(String)}
     */
    public final String name;
    /**
     * A {@link HashMap} yielding any meta-information {@link String} key-value pairs about this {@link AlleleEntry}.
     */
    public final HashMap<String, String> annotations = new HashMap<>();
    /**
     * {@link TreeSet} of {@link String} values matching exactly one {@link SampleEntry#name} value. Represents the list of
     * all samples yielding this allele.
     */
    public final TreeSet<String> samples = new TreeSet<>();

    public final static String PROPERTY_NAME_REFERENCE_ID = "AL00000000WT";

    public final static String PROPERTY_NAME_VARIANTS = "VARIANTS";
    public final static String PROPERTY_NAME_VARIABLE_POSITIONS = "VAR_POS";
    public final static String PROPERTY_NAME_NUMBER_OF_SUBSTITUTIONS = "NO_SUB";
    public final static String PROPERTY_NAME_NUMBER_OF_INSERTIONS = "NO_INS";
    public final static String PROPERTY_NAME_NUMBER_OF_DELETIONS = "NO_DEL";
    public final static String PROPERTY_NAME_CONGLOMERATION_INDEX = "CONG_IDX";
    public final static String PROPERTY_NAME_FREQUENCY = "FRQ";

    /**
     * Constructor of {@link AlleleEntry}.
     *
     * @param name           The internal name to use for the allele.
     * @param sId            The sample ID/name for which the allele was derived.
     * @param concatVariants Mandatory annotation to use for this allele.
     */
    public AlleleEntry(String name, String sId, String concatVariants) {
        this.name = name;
        this.annotations.put(PROPERTY_NAME_VARIANTS, concatVariants);
        this.samples.add(sId);
    }

    /**
     * Returns a {@link String} that is used as internal genotype name.
     * <p>
     * See {@link FeatureEntry#generateEntryName(String, String)} for more details.
     *
     * @param concatVariants {@link String} representation of the variants of one allele in the format
     *                       <POS0>_<ALT0>;...;<POSn>_<ALTn>.
     * @return {@link String} intended to be used as internal allele name.
     */
    public static String generateAlleleName(String concatVariants) {
        if (concatVariants.equals("")) {
            return AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
        } else {
            return FeatureEntry.generateEntryName(concatVariants, "AL");
        }
    }

}
