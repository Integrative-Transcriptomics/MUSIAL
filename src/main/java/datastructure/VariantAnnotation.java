package datastructure;

import main.MusialConstants;

import java.util.NavigableSet;
import java.util.concurrent.ConcurrentSkipListMap;

/**
 * Container to store annotations of a single variant.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.1
 */
@SuppressWarnings("unused")
public class VariantAnnotation {

    /**
     * {@link ConcurrentSkipListMap} of {@link String} key/value pairs specifying any information related to the variant.
     */
    private final ConcurrentSkipListMap<String, String> properties = new ConcurrentSkipListMap<>();

    /**
     * Constructor of {@link VariantAnnotation}.
     */
    public VariantAnnotation() {
    }

    /**
     * Puts a meta-information key/value pair into this {@link VariantAnnotation#properties}.
     *
     * @param key   {@link String}
     * @param value {@link String}
     */
    public void addProperty(String key, String value) {
        this.properties.put(key, value);
    }

    /**
     * Appends an existing or puts a meta-information key/value pair into this {@link VariantAnnotation#properties}.
     * <p>
     * For appending, the passed value is attached to the currently stored value separated by `:`.
     *
     * @param key   {@link String}
     * @param value {@link String}
     */
    public void appendProperty(String key, String value) {
        if (this.properties.containsKey(key)) {
            this.properties.put(key, this.properties.get(key) + MusialConstants.FIELD_SEPARATOR_1 + value);
        } else {
            addProperty(key, value);
        }
    }

    /**
     * Removes a meta-information key/value pair into this {@link VariantAnnotation#properties}. The stored value is returned.
     *
     * @param key {@link String}
     */
    public void removeProperty(String key) {
        this.properties.remove(key);
    }

    /**
     * Retrieves the value stored under key in this {@link VariantAnnotation#properties}.
     *
     * @param key {@link String}
     * @return {@link String} value stored under key.
     */
    public String getProperty(String key) {
        return this.properties.get(key);
    }

    /**
     * Retrieves the keys stored in this {@link VariantAnnotation#properties}.
     *
     * @return {@link NavigableSet} of stored property keys.
     */
    public NavigableSet<String> getPropertyKeys() {
        return this.properties.keySet();
    }

    /**
     * Check if element is contained in this {@link VariantAnnotation#properties}.
     *
     * @param key {@link String}
     * @return {@link Boolean}
     */
    public boolean hasProperty(String key) {
        return this.properties.containsKey(key);
    }


}