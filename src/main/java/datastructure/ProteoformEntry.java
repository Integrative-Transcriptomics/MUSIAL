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

    /**
     * Constructor of {@link ProteoformEntry}.
     *
     * @param name         The internal name to use for the proteoform.
     * @param sampleId     The sample ID/name for which the proteoform was derived.
     * @param variantsSwab Mandatory VARIANTS annotation to use for this proteoform.
     */
    public ProteoformEntry(String name, String sampleId, String variantsSwab) {
        this.name = name;
        this.annotations.put(ProteoformEntry.PROPERTY_NAME_VARIANTS, variantsSwab);
        this.samples.add(sampleId);
    }

    /**
     * Returns a {@link String} that is used as internal proteoform name.
     * <p>
     * These names are generated by concatenating the prefix 'PF' with the HashCode of the {@link ProteoformEntry}
     * mandatory 'VARIANTS' annotation value converted to a {@link String} for which a negative sign is converted to a one
     * and which in total was padded with prepending 0s to a length of 10.
     *
     * @param variantsSwab The {@link String} value of the 'VSWAB' annotation of a {@link ProteoformEntry}.
     * @return {@link String} intended to be used as internal proteoform name.
     */
    public static String generateProteoformName(String variantsSwab) {
        if (variantsSwab.equals("")) {
            return ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
        } else {
            StringBuilder proteoformNameBuilder = new StringBuilder();
            proteoformNameBuilder.append("PF");
            String hashCodeString = String.valueOf(variantsSwab.hashCode());
            if (hashCodeString.startsWith("-")) {
                proteoformNameBuilder.append("1");
                hashCodeString = hashCodeString.replace("-", "");
            } else {
                proteoformNameBuilder.append("0");
            }
            proteoformNameBuilder.append("0".repeat(10 - hashCodeString.length()));
            proteoformNameBuilder.append(hashCodeString);
            return proteoformNameBuilder.toString();
        }
    }

}
