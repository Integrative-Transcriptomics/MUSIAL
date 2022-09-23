package datastructure;

import java.util.HashMap;
import java.util.concurrent.ConcurrentSkipListMap;

/**
 * Object used to store annotations for a single nucleotide variant.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class NucleotideVariantEntry {

    /**
     * {@link HashMap} of {@link String} key/value pairs specifying any information related to the variant.
     */
    public final HashMap<String, String> annotations = new HashMap<>();

    /**
     * {@link ConcurrentSkipListMap} of {@link String} key/value pairs indicating the occurrence of this variant.
     * <p>
     * - Keys correspond to sample names.
     * - Values correspond to return values of the {@link NucleotideVariantEntry#constructSampleSpecificAnnotation} method.
     */
    public final ConcurrentSkipListMap<String, String> occurrence = new ConcurrentSkipListMap<>();

    /**
     * Constructor of {@link AminoacidVariantEntry}.
     */
    public NucleotideVariantEntry() {

    }

    /**
     * Constructs a {@link String} yielding information about the following properties wrt. a single sample, separated by a `|` symbol:
     * <p>
     * - If the variant was rejected, i.e. failed any filter criteria.
     * - If the variant is primary, i.e. has the highest frequency in the case of a het. variant.
     * - The quality of the variant call.
     * - The frequency wrt. coverage of the variant call.
     * - The total coverage at the variant site.
     *
     * @param rejected  {@link Boolean} whether the variant was rejected.
     * @param primary   {@link Boolean} whether the variant is a primary variant.
     * @param quality   {@link Double}; The quality of the variant call.
     * @param frequency {@link Double}; The frequency of the variant call wrt. coverage.
     * @param coverage  {@link Double}; The total coverage at the variant site.
     * @return {@link String}; The passed parameters separated by `|` symbols.
     */
    public static String constructSampleSpecificAnnotation(boolean rejected, boolean primary, double quality,
                                                           double frequency, double coverage) {
        return rejected + "|" + primary + "|" + quality + "|" + frequency + "|" + coverage;
    }

}