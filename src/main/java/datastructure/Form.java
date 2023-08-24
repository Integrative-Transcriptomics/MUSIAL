package datastructure;

import java.util.HashMap;
import java.util.TreeSet;

/**
 * Container to store a biological variation of a feature, i.e., an allele of a gene or proteoform of a protein.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.2
 */
@SuppressWarnings("unused")
public class Form {

    /**
     * The internal name to use for this form. This attribute is inferred from the variants that define this instance.
     */
    public final String name;
    /**
     * Stores meta-information as key-value pairs about this instance.
     */
    private final HashMap<String, String> annotations = new HashMap<>();
    /**
     * Stores all occurrences of this form in, e.g., variant call sets (samples) or alleles.
     */
    private final TreeSet<String> occurrence = new TreeSet<>();

    /**
     * Constructor of {@link Form}.
     *
     * @param name Name to use for this object.
     */
    public Form(String name) {
        this.name = name;
    }

    /**
     * Add an annotation to this {@link #annotations}.
     *
     * @param key   Key of the annotation to add.
     * @param value Value of the annotation to add.
     */
    public void addAnnotation(String key, String value) {
        this.annotations.put(key, value);
    }

    /**
     * Remove an annotation from {@link #annotations}, if it exists.
     *
     * @param key Key of the annotation to remove.
     */
    public void removeAnnotation(String key) {
        this.annotations.remove(key);
    }

    /**
     * Query an annotation value from {@link #annotations} or null, if it does not exist.
     *
     * @param key Key of the annotation to query.
     * @return Value of the annotation to query or null.
     */
    public String getAnnotation(String key) {
        return this.annotations.get(key);
    }

    /**
     * Check if an annotation in {@link #annotations} exists.
     *
     * @param key Key of the annotation to query.
     * @return True, if an annotation is stored with the specified key.
     */
    public boolean hasAnnotation(String key) {
        return this.annotations.containsKey(key);
    }

    /**
     * Add an entry to this {@link #occurrence}.
     *
     * @param occurrenceName Name of the entry to add.
     */
    public void addOccurrence(String occurrenceName) {
        this.occurrence.add(occurrenceName);
    }

    /**
     * Remove an entry from {@link #occurrence}, if it exists.
     *
     * @param occurrenceName Name of the entry to remove.
     */
    public void removeOccurrence(String occurrenceName) {
        this.occurrence.remove(occurrenceName);
    }

    /**
     * @return The {@link #occurrence} of this instance.
     */
    public TreeSet<String> getOccurrence() {
        return this.occurrence;
    }

    /**
     * Check if an entry in {@link #occurrence} exists.
     *
     * @param occurrenceName Name of the entry to query.
     * @return True, if an entry is stored with the specified name.
     */
    public boolean hasOccurrence(String occurrenceName) {
        return this.occurrence.contains(occurrenceName);
    }
}
