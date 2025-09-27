package model;

import util.Constants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Base class for entities to store arbitrary attributes as {@link String}s.
 * <p>
 * This class provides methods to manage attributes associated with an entity. Attributes are stored as key-value pairs in a {@link Map},
 * allowing efficient retrieval, addition, extension, and removal of such. It also supports operations like checking for the existence of
 * attributes and converting attributes to a string representation.
 */
public class Attributes {

    /**
     * Attributes associated with this entity, stored as key-value pairs.
     */
    private final Map<String, String> attributes = new HashMap<>();

    /**
     * Constructor of {@link Attributes}.
     * <p>
     * Initializes an empty attributes map for the entity.
     */
    Attributes() {
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
     * Checks if this entity has any attribute.
     *
     * @return {@code true} if at least one attribute exists, {@code false} otherwise.
     */
    public boolean hasAnyAttribute() {
        return !this.attributes.isEmpty();
    }

    /**
     * Retrieves the value of an attribute associated with this entity. If the attribute does not exist, {@link Constants#EMPTY} is
     * returned.
     *
     * @param key The key of the attribute to retrieve.
     * @return The value of the attribute, or {@link Constants#EMPTY} if the attribute does not exist.
     */
    public String getAttribute(String key) {
        return this.attributes.getOrDefault(key, Constants.EMPTY);
    }

    /**
     * Retrieves the value of an attribute associated with this entity. If the attribute does not exist, the specified default value is
     * returned.
     *
     * @param key          The key of the attribute to retrieve.
     * @param defaultValue The default value to return if the attribute does not exist.
     * @return The value of the attribute, or the specified default value if the attribute does not exist.
     */
    public String getAttributeOrDefault(String key, String defaultValue) {
        return this.attributes.getOrDefault(key, defaultValue);
    }

    /**
     * Retrieves the value of an attribute as a collection of strings.
     * <p>
     * This method fetches the value of the specified attribute key from the attributes map. The value is split into individual elements
     * using a comma (`,`) as the delimiter. If the attribute does not exist, an empty collection is returned.
     *
     * @param key The key of the attribute to retrieve.
     * @return A {@link Set} of strings representing the split values of the attribute, or an empty set if the attribute does not exist.
     */
    public Set<String> getAttributeSet(String key) {
        return Arrays.stream(this.attributes.getOrDefault(key, Constants.EMPTY).split(Constants.COMMA))
                .collect(Collectors.toUnmodifiableSet());
    }

    /**
     * Retrieves all attributes associated with this entity.
     * <p>
     * This method provides an unmodifiable view of the attributes map, ensuring that the original map cannot be modified externally. The
     * attributes are stored as key-value pairs, where both the key and value are {@link String}.
     *
     * @return An unmodifiable {@link Map} containing all attributes associated with this entity.
     */
    public Map<String, String> getAttributes() {
        return Collections.unmodifiableMap(this.attributes);
    }

    /**
     * Converts the attributes of this entity to a string representation..
     *
     * @param separator The separator to use between key-value pairs.
     * @return A string representation of the attributes.
     */
    public String attributesAsString(String separator) {
        return this.attributes.entrySet().stream()
                .map(entry -> entry.getKey() + Constants.EQUAL + entry.getValue())
                .collect(Collectors.joining(separator));
    }

    /**
     * Converts the attributes of this entity to a string representation, excluding attributes with keys in the specified collection.
     *
     * @param except    A collection of keys to exclude from the string representation.
     * @param separator The separator to use between key-value pairs.
     * @return A string representation of the attributes, excluding the specified keys.
     */
    public String attributesAsString(Collection<String> except, String separator) {
        return this.attributes.entrySet().stream()
                .filter(entry -> !except.contains(entry.getKey()))
                .map(entry -> entry.getKey() + Constants.EQUAL + entry.getValue())
                .collect(Collectors.joining(separator));
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
    public void setAttributeIfAbsent(String key, String value) {
        String result = this.attributes.putIfAbsent(key, value);
    }

    /**
     * Adds multiple attributes to this entity only if they do not already exist.
     *
     * @param attributes A map of attributes to associate with this entity.
     */
    public void setAttributesIfAbsent(Map<String, String> attributes) {
        attributes.forEach(this::setAttributeIfAbsent);
    }

    /**
     * Extends an attribute by appending a value to the existing value. Commas are used to separate values. If the value already exists, it
     * will not be added again.
     *
     * @param key   The key of the attribute.
     * @param value The value to append to the attribute.
     */
    public void extendAttribute(String key, String value) {
        String currentValue = this.attributes.get(key);
        if (currentValue != null) {
            if (!currentValue.contains(value)) {
                this.setAttribute(key, currentValue + Constants.COMMA + value);
            }
        } else {
            this.setAttribute(key, value);
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

}