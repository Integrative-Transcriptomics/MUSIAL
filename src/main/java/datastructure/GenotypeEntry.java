package datastructure;

import java.util.HashMap;
import java.util.TreeSet;

/**
 * Represents a proteoform, i.e. a set of single-position amino-acid variants (may include inserted and deleted
 * positions), wrt. a {@link FeatureEntry#allocatedProtein}.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class GenotypeEntry {

    /**
     * The internal name to use for this proteoform. Should be derived with {@link AllocatedProteinEntry#generateProteoformName(String)}
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

    /**
     * Constructor of {@link GenotypeEntry}.
     *
     * @param name            The internal name to use for the proteoform.
     * @param sId             The sample ID/name for which the proteoform was derived.
     * @param proteoformVSwab The mandatory VSWAB annotation to use for this proteoform.
     */
    public GenotypeEntry(String name, String sId, String proteoformVSwab) {
        this.name = name;
        this.annotations.put("VSWAB", proteoformVSwab);
        this.samples.add(sId);
    }

}
