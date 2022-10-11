package datastructure;

import java.util.HashMap;
import java.util.concurrent.ConcurrentSkipListSet;

/**
 * Object used to store annotations for a single amino-acid variant.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class AminoacidVariantEntry {

    /**
     * {@link HashMap} of {@link String} key/value pairs specifying any information related to the variant.
     */
    public final HashMap<String, String> annotations = new HashMap<>();

    /**
     * {@link ConcurrentSkipListSet} of {@link String}s indicating the occurrence of this variant.
     */
    public final ConcurrentSkipListSet<String> occurrence = new ConcurrentSkipListSet<>();

    public final static String PROPERTY_NAME_FREQUENCY = "FRQ";

    /**
     * Constructor of {@link AminoacidVariantEntry}.
     */
    public AminoacidVariantEntry() {

    }

}
