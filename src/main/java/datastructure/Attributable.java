package datastructure;

import utility.Constants;

import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

/**
 * Base class for entities that can have attributes.
 * <p>
 * This class provides methods to manage attributes associated with an entity. Attributes are stored
 * as key-value pairs in a {@link TreeMap}, allowing efficient retrieval, addition, extension, and removal
 * of attributes. It also supports operations like checking for the existence of attributes and converting
 * attributes to a string representation.
 */
public class Attributable {

    /**
     * Key used to represent the occurrence of samples in classes that extend {@link Attributable}.
     */
    public final static String sampleOccurrence = "sample";
    /**
     * Attributes associated with this entity, stored as key-value pairs.
     */
    private final TreeMap<String, String> attributes = new TreeMap<>();

    /**
     * Constructor of {@link Attributable}.
     * <p>
     * Initializes an empty attributes map for the entity.
     */
    protected Attributable() {
    }

    /**
     * Adds an attribute to this entity. If an attribute with the same key already exists, it will be overwritten.
     *
     * @param key   The key of the attribute.
     * @param value The value of the attribute.
     */
    public void setAttribute(String key, String value) {
        this.attributes.put(key, value);
    }

    /**
     * Adds multiple attributes to this entity. Existing attributes with the same keys will be overwritten.
     *
     * @param attributes A map of attributes to associate with this entity.
     */
    public void setAttributes(Map<String, String> attributes) {
        attributes.forEach(this::setAttribute);
    }

    /**
     * Adds an attribute to this entity only if it does not already exist.
     *
     * @param key   The key of the attribute.
     * @param value The value of the attribute.
     */
    public void addAttributeIfAbsent(String key, String value) {
        this.attributes.putIfAbsent(key, value);
    }

    /**
     * Adds multiple attributes to this entity only if they do not already exist.
     *
     * @param attributes A map of attributes to associate with this entity.
     */
    public void addAttributesIfAbsent(Map<String, String> attributes) {
        attributes.forEach(this::addAttributeIfAbsent);
    }

    /**
     * Extends an attribute by appending a value to the existing value. Commas are used to separate values.
     * If the value already exists, it will not be added again.
     *
     * @param key   The key of the attribute.
     * @param value The value to append to the attribute.
     */
    public void extendAttribute(String key, String value) {
        String currentValue = this.attributes.get(key);
        if (currentValue != null) {
            if (!currentValue.contains(value)) {
                this.attributes.put(key, currentValue + Constants.COMMA + value);
            }
        } else {
            this.attributes.put(key, value);
        }
    }

    /**
     * Extends multiple attributes by appending values to the existing values. Commas are used to separate values.
     *
     * @param attributes A map of attributes to extend.
     */
    public void extendAttributes(Map<String, String> attributes) {
        attributes.forEach(this::extendAttribute);
    }

    /**
     * Retrieves the value of an attribute associated with this entity. If the attribute does not exist,
     * {@link Constants#EMPTY} is returned.
     *
     * @param key The key of the attribute to retrieve.
     * @return The value of the attribute, or {@link Constants#EMPTY} if the attribute does not exist.
     */
    public String getAttribute(String key) {
        return this.attributes.getOrDefault(key, Constants.EMPTY);
    }

    /**
     * Retrieves the value of an attribute as a collection of strings.
     * <p>
     * The attribute value is split into individual elements using the comma (`,`) as a delimiter.
     * If the attribute does not exist, an empty collection is returned.
     *
     * @param key The key of the attribute to retrieve.
     * @return A collection of strings representing the split values of the attribute, or an empty collection if the attribute does not exist.
     */
    public Collection<String> getAttributeAsCollection(String key) {
        return Arrays.stream(this.attributes.getOrDefault(key, Constants.EMPTY).split(Constants.COMMA))
                .collect(Collectors.toSet());
    }

    /**
     * Retrieves all attributes associated with this entity.
     *
     * @return A map of all attributes.
     */
    public Map<String, String> getAttributes() {
        return this.attributes;
    }

    /**
     * Checks if an attribute with the specified key exists in this entity.
     *
     * @param key The key of the attribute to query.
     * @return {@code true} if the attribute exists, {@code false} otherwise.
     */
    public boolean hasAttribute(String key) {
        return this.attributes.containsKey(key);
    }

    /**
     * Checks if this entity has any attributes.
     *
     * @return {@code true} if at least one attribute exists, {@code false} otherwise.
     */
    public boolean hasAttributes() {
        return !this.attributes.isEmpty();
    }

    /**
     * Removes an attribute with the specified key from this entity.
     *
     * @param key The key of the attribute to remove.
     */
    public void removeAttribute(String key) {
        this.attributes.remove(key);
    }

    /**
     * Removes all attributes from this entity.
     */
    public void clearAttributes() {
        this.attributes.clear();
    }

    /**
     * Converts the attributes of this entity to a string representation in the format {@code KEY=VALUE;...}.
     *
     * @return A string representation of the attributes.
     */
    public String attributesAsString() {
        return this.attributes.entrySet().stream()
                .map(entry -> entry.getKey() + Constants.EQUAL + entry.getValue())
                .collect(Collectors.joining(Constants.SEMICOLON));
    }

    /**
     * Converts the attributes of this entity to a string representation in the format {@code KEY=VALUE;...},
     * excluding attributes with keys in the specified collection.
     *
     * @param except A collection of keys to exclude from the string representation.
     * @return A string representation of the attributes, excluding the specified keys.
     */
    public String attributesAsString(Collection<String> except) {
        return this.attributes.entrySet().stream()
                .filter(entry -> !except.contains(entry.getKey()))
                .map(entry -> entry.getKey() + Constants.EQUAL + entry.getValue())
                .collect(Collectors.joining(Constants.SEMICOLON));
    }
}