package model;

import exceptions.MusialException;
import util.Constants;

import java.util.*;

/**
 * Representation of a genomic feature.
 * <p>
 * This class models a genomic feature, such as a gene, exon, or coding sequence (CDS), that is analyzed in the context of genomic data
 * processing. It extends the {@link Attributes} class to inherit functionality for managing attributes associated with the feature. Each
 * instance contains information about its type and location on the reference genome - i.e., a {@link Contig}. In addition, an inner map is
 * used to store {@link SequenceType}s associated with the feature.
 * <p>
 * Features are stored in the {@link Storage#features} property of the model.
 */
public class Feature extends Attributes {

    /**
     * Unique identifier of this feature.
     * <p>
     * This field serves as the unique identifier for the feature and is used to reference it in the model.
     */
    public final String _id;

    /**
     * The type of this genomic feature.
     * <p>
     * This field specifies the type of the feature, such as "gene", "exon", or "CDS". The type is defined according to the GFF3
     * specification and provides information about the biological or functional classification of the feature.
     * <p>
     * For more details, refer to the GFF3 specification:
     * <a href="https://gmod.org/wiki/GFF3">https://gmod.org/wiki/GFF3</a>.
     */
    public final String type;

    /**
     * An additional (human-readable) name for this feature.
     * <p>
     * This should be at best a database identifier or common gene name.
     */
    public final String name;

    /**
     * The parent genomic location of this feature.
     */
    public final String contig;

    /**
     * 1-based starting position of the feature.
     */
    public final int start;

    /**
     * 1-based end position of the feature.
     */
    public final int end;

    /**
     * Strand of the feature.
     */
    public final char strand;

    /**
     * A list of sub-features associated with this genomic feature.
     * <p>
     * This list stores {@link SubFeature} objects, each representing a sub-feature of the genomic feature, such as exons or coding
     * sequences (CDS). Sub-features are defined by their type and genomic location (start and end positions).
     */
    private final List<SubFeature> subFeatures = new ArrayList<>();

    /**
     * Represents a sub-feature of a genomic feature.
     * <p>
     * This record encapsulates the type and genomic location of a sub-feature, such as an exon or CDS, within a parent genomic feature. It
     * includes:
     * <ul>
     *   <li>The type of the sub-feature (e.g., "exon", "CDS").</li>
     *   <li>The 1-based start position of the sub-feature on the reference genome.</li>
     *   <li>The 1-based end position of the sub-feature on the reference genome.</li>
     * </ul>
     *
     * @param type  The type of the sub-feature (e.g., "exon", "CDS").
     * @param start The 1-based start position of the sub-feature.
     * @param end   The 1-based end position of the sub-feature.
     */
    public record SubFeature(String type, int start, int end) {

        /**
         * Generates a GFF3 conform identifier (ID) for a sub-feature based on its type and a parent identifier.
         * <p>
         * This method constructs an ID string for the sub-feature using its type and the provided parent identifier. The format of the ID
         * depends on the type of the sub-feature:
         * <ul>
         *   <li>If the type is "CDS", the ID is formatted as "ID=cds-{parentIdentifier};Parent=transcript-{parentIdentifier}".</li>
         *   <li>If the type is "exon", the ID is formatted as "ID=exon-{parentIdentifier};Parent=transcript-{parentIdentifier}".</li>
         *   <li>If the type contains "RNA", the ID is formatted as "ID=transcript-{parentIdentifier};Parent={parentIdentifier}".</li>
         *   <li>For other types, the ID is formatted as "ID={type}-{parentIdentifier};Parent={parentIdentifier}".</li>
         * </ul>
         *
         * @param parentIdentifier The identifier of the parent feature.
         * @return A formatted string representing the GFF3 ID value of the sub-feature.
         */
        public String ID(String parentIdentifier) {
            return switch (this.type) {
                case "CDS" -> "ID=cds-%s;Parent=transcript-%s".formatted(parentIdentifier, parentIdentifier);
                case "exon" -> "ID=exon-%s;Parent=transcript-%s".formatted(parentIdentifier, parentIdentifier);
                default -> this.type.contains("RNA")
                        ? "ID=transcript-%s;Parent=%s".formatted(parentIdentifier, parentIdentifier)
                        : "ID=%s-%s;Parent=%s".formatted(this.type, parentIdentifier, parentIdentifier);
            };
        }

    }

    /**
     * Alleles ({@link SequenceType} instances) associated with this feature.
     * <p>
     * This map stores alleles that are associated with the feature. Alleles represent specific sequence variations of the feature. The keys
     * in the map are unique identifiers for the alleles, and the values are the corresponding {@link Allele} instances.
     * <p>
     * Alleles are used to track and manage sequence variations resulting from genomic changes. Each allele is linked to its unique
     * identifier and contains information about its sequence and attributes.
     */
    private final Map<String, Allele> alleles = new HashMap<>();

    /**
     * Proteoforms ({@link SequenceType} instances) associated with this feature.
     * <p>
     * This map stores proteoforms that are associated with the feature. Proteoforms represent specific sequence variants of proteins
     * derived from the feature. The keys in the map are unique identifiers for the proteoforms, and the values are the corresponding
     * {@link Proteoform} instances.
     * <p>
     * Proteoforms are only relevant for coding features and are used to track and manage protein sequence variations resulting from genomic
     * changes.
     */
    private final Map<String, Proteoform> proteoforms = new HashMap<>();

    /**
     * Constructs a new {@link Feature} instance with the specified properties.
     * <p>
     * This constructor initializes a genomic feature with its id, location, strand orientation, type, and unique identifier. The feature's
     * start and end positions are converted to integers to ensure proper indexing. The {@link Attributes} superclass is also initialized.
     *
     * @param name       The id of the feature, used as its internal identifier.
     * @param contig     The id of the reference location (e.g., contig, chromosome, plasmid) where the feature is located.
     * @param start      The 1-based indexed starting position of the feature on the reference.
     * @param end        The 1-based indexed end position of the feature on the reference.
     * @param strand     The strand orientation of the feature ('+' for forward strand, '-' for reverse strand).
     * @param type       The type of the feature (e.g., coding, non-coding).
     * @param identifier The unique identifier of the feature.
     */
    public Feature(String name, String contig, Number start, Number end, char strand, String type, String identifier) {
        super();
        this.name = name;
        this.contig = contig;
        this.type = type;
        this.start = start.intValue();
        this.end = end.intValue();
        this.strand = strand;
        this._id = identifier;
    }

    /**
     * Determines if this feature is a coding feature.
     * <p>
     * This method checks whether the feature is of type "CDS" (coding sequence) or if any of its sub-features are of type "CDS". A feature
     * is considered coding if it directly represents a coding sequence or contains sub-features that do.
     *
     * @return {@code true} if the feature is of type "CDS" or has sub-features of type "CDS"; {@code false} otherwise.
     */
    public boolean isCoding() {
        return this.type.equals("CDS") || subFeatures.stream().anyMatch(sf -> sf.type.equals("CDS"));
    }

    /**
     * Returns if this feature is on the reverse strand.
     * <p>
     * This method determines whether the strand orientation of the feature is reverse by checking if the strand character is {@code '-'}.
     *
     * @return {@code true} if this feature is on the reverse strand, {@code false} otherwise.
     */
    public boolean isReverse() {
        return this.strand == '-';
    }

    /**
     * Adds a sub-feature to this genomic feature.
     * <p>
     * This method validates and adds a sub-feature to the list of sub-features associated with this genomic feature. The sub-feature is
     * defined by its type, start position, and end position. Validation ensures that:
     * <ul>
     *   <li>The sub-feature type is recognized in the {@link Storage#SEQUENCE_ONTOLOGY_HIERARCHY} map.</li>
     *   <li>The sub-feature's start and end positions are within the bounds of the parent feature.</li>
     * </ul>
     * If validation fails, an {@link IllegalArgumentException} is thrown.
     *
     * @param type  The type of the sub-feature (e.g., "exon", "CDS").
     * @param start The 1-based start position of the sub-feature.
     * @param end   The 1-based end position of the sub-feature.
     * @throws MusialException if the sub-feature type is unrecognized or if its positions are out of bounds.
     */
    public void addSubFeature(String type, int start, int end) throws MusialException {
        if (!Storage.SEQUENCE_ONTOLOGY_HIERARCHY.containsKey(type))
            throw new MusialException("Sub-features of type '%s' are not recognized.".formatted(type));
        if (start < this.start || end > this.end)
            throw new MusialException("Sub-feature %s (%s:g.%d_%d=) is out of bounds of its parent feature %s (%s:g.%d_%d=)."
                    .formatted(type, this.contig, start, end, this.name, this.contig, this.start, this.end));
        if (start > end)
            throw new MusialException("Sub-feature %s (%s:g.%d_%d=) has an invalid location (start > end)."
                    .formatted(type, this.contig, start, end));
        this.subFeatures.add(new SubFeature(type, start, end));
        // Sort the sub-features based on their Sequence Ontology hierarchy level after adding a new one.
        // Note: The repeated sorting could be optimized if performance becomes an issue.
        this.subFeatures.sort(Comparator.comparingInt(sf -> Storage.SEQUENCE_ONTOLOGY_HIERARCHY.get(sf.type)));
    }

    /**
     * Retrieves all sub-features associated with this genomic feature.
     * <p>
     * This method provides an unmodifiable view of the list of sub-features associated with this genomic feature. Sub-features represent
     * smaller components of the feature, such as exons or coding sequences (CDS), and include their type and genomic location (start and
     * end positions).
     * <p>
     * The unmodifiable list ensures that the original list cannot be modified externally, preserving data integrity.
     *
     * @return An unmodifiable {@link List} of {@link SubFeature} objects representing the sub-features of this genomic feature.
     */
    public List<SubFeature> getSubFeatures() {
        return Collections.unmodifiableList(this.subFeatures);
    }

    /**
     * Clears all sub-features associated with this genomic feature.
     * <p>
     * This method removes all sub-features from the list of sub-features associated with this genomic feature. After calling this method,
     * the list of sub-features will be empty.
     */
    public void clearSubFeatures() {
        this.subFeatures.clear();
    }

    /**
     * Clears all sub-features of a specific Sequence Ontology (SO) hierarchy level associated with this genomic feature.
     * <p>
     * This method removes all sub-features from the list of sub-features that match the specified SO hierarchy level. The hierarchy level
     * is determined using the {@link Storage#SEQUENCE_ONTOLOGY_HIERARCHY} map.
     * <p>
     * After calling this method, only sub-features that do not match the specified level will remain in the list.
     *
     * @param level The SO hierarchy level of the sub-features to remove (e.g., 0 for "region", 1 for "gene").
     */
    public void clearSubFeatures(int level) {
        this.subFeatures.removeIf(sf -> Storage.SEQUENCE_ONTOLOGY_HIERARCHY.get(sf.type) == level);
    }

    /**
     * Checks if an allele with the specified identifier exists in this feature.
     *
     * @param alleleIdentifier The identifier of the allele to check for.
     * @return {@code true} if an allele with the given identifier exists, {@code false} otherwise.
     */
    public boolean hasAllele(String alleleIdentifier) {
        return this.alleles.containsKey(alleleIdentifier);
    }

    /**
     * Adds an allele to this feature.
     * <p>
     * This method adds the specified {@link Allele} object to the internal map of alleles associated with this feature. The allele is
     * stored using its unique identifier as the key.
     * <p>
     * <i>Note: No internal validation is performed based on the coordinates of the feature.</i>
     *
     * @param allele The {@link Allele} object to be added to this feature.
     */
    public void addAllele(Allele allele) {
        this.alleles.put(allele._id, allele);
    }

    /**
     * Retrieves an allele associated with this feature by its unique identifier (_id) or {@code null}.
     *
     * @param alleleIdentifier The identifier of the allele to retrieve.
     * @return The {@link Allele} object associated with the given UID or {@code null} if not found.
     */
    public Allele getAllele(String alleleIdentifier) {
        return alleles.getOrDefault(alleleIdentifier, null);
    }

    /**
     * Removes an allele associated with this feature.
     * <p>
     * This method removes the allele identified by the given identifier from the internal map of alleles associated with this feature. If
     * the identifier does not exist in the map, no action is performed.
     *
     * @param alleleIdentifier The unique identifier of the allele to be removed.
     */
    void removeAllele(String alleleIdentifier) {
        this.alleles.remove(alleleIdentifier);
    }

    /**
     * Retrieves all alleles associated with this feature.
     * <p>
     * This method provides an unmodifiable view of the collection of alleles associated with this feature. The alleles are stored as values
     * in the internal map, ensuring that the original collection cannot be modified externally.
     *
     * @return An unmodifiable {@link Collection} of {@link Allele} objects associated with this feature.
     */
    public Collection<Allele> getAlleles() {
        return Collections.unmodifiableCollection(this.alleles.values());
    }

    /**
     * Retrieves the number of alleles associated with this feature.
     * <p>
     * This method returns the size of the internal map of alleles, which represents the total number of unique alleles associated with this
     * feature.
     *
     * @return The number of alleles associated with this feature.
     */
    public int getAlleleCount() {
        return this.alleles.size();
    }

    /**
     * Checks if a proteoform with the specified identifier exists in this feature.
     *
     * @param proteoformIdentifier The identifier of the proteoform to check for.
     * @return {@code true} if a proteoform with the given identifier exists, {@code false} otherwise.
     */
    public boolean hasProteoform(String proteoformIdentifier) {
        return this.proteoforms.containsKey(proteoformIdentifier);
    }

    /**
     * Adds a proteoform to this feature.
     * <p>
     * This method adds the specified {@link Proteoform} object to the internal map of proteoforms associated with this feature. The
     * proteoform is stored using its unique identifier as the key.
     * <p>
     * <i>Note: No internal validation is performed based on the coordinates of the feature.</i>
     *
     * @param proteoform The {@link Proteoform} object to be added to this feature.
     */
    public void addProteoform(Proteoform proteoform) {
        this.proteoforms.put(proteoform._id, proteoform);
    }

    /**
     * Retrieves a proteoform associated with this feature by its unique identifier or {@code null}.
     *
     * @param proteoformIdentifier The unique identifier of the proteoform to retrieve.
     * @return The {@link Proteoform} object associated with the given UID or {@code null} if not found.
     */
    public Proteoform getProteoform(String proteoformIdentifier) {
        return proteoforms.getOrDefault(proteoformIdentifier, null);
    }

    /**
     * Retrieves all proteoforms associated with this feature.
     * <p>
     * This method provides an unmodifiable view of the collection of proteoforms associated with this feature. Proteoforms represent
     * specific sequence variants of proteins derived from the feature. The unmodifiable collection ensures that the original data cannot be
     * modified externally.
     *
     * @return An unmodifiable {@link Collection} of {@link Proteoform} objects associated with this feature.
     */
    public Collection<Proteoform> getProteoforms() {
        return Collections.unmodifiableCollection(this.proteoforms.values());
    }

    /**
     * Retrieves the number of proteoforms associated with this feature.
     * <p>
     * This method returns the size of the internal map of proteoforms, which represents the total number of unique proteoforms associated
     * with this feature.
     *
     * @return The number of proteoforms associated with this feature.
     */
    public int getProteoformCount() {
        return this.proteoforms.size();
    }

    /**
     * Generates a string representation of this feature.
     * <p>
     * This method constructs a string representation of the feature using its contig, genomic coordinates, and name. The format includes
     * the contig identifier, start and end positions, and the feature name, separated by specific constants.
     *
     * @return A {@link String} representing the feature in the format: {@code contig:g.start_end=featureName}.
     */
    public String toString() {
        return "%s%s%s%d%s%d%s%s".formatted(contig, Constants.COLON, Constants.GENOMIC_COORDINATES_PREFIX, start, Constants.UNDER_SCORE,
                end, Constants.EQUAL, name);
    }

    /**
     * Computes the hash code for this feature.
     * <p>
     * This method calculates the hash code of the feature based on its string representation.
     *
     * @return The hash code of this feature.
     */
    public int hashCode() {
        return this.toString().hashCode();
    }

    /**
     * Compares this feature to another object for equality.
     * <p>
     * This method checks if the provided object is the same instance as this feature. If not, it verifies that the object is of the same
     * class and compares their string representations for equality.
     *
     * @param obj The object to compare with this {@link Feature} instance.
     * @return {@code true} if the objects are the same instance or if their string representations are equal; {@code false} otherwise.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        Feature that = (Feature) obj;
        return this.toString().equals(that.toString());
    }

}