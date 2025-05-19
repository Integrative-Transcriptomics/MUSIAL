package datastructure;

import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import main.Musial;
import utility.Constants;
import utility.IO;
import utility.Logging;
import utility.SequenceOperations;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static utility.Constants.DOT;
import static utility.Constants.TAB;

/**
 * Representation of a genomic feature that is subject to analysis.
 * <p>
 * This class models a genomic feature, such as a gene, exon, or coding sequence (CDS),
 * that is analyzed in the context of genomic data processing. It extends the {@link Attributable}
 * class to inherit functionality for managing attributes associated with the feature.
 * <p>
 * Each instance of this class is uniquely identified by its {@code name} and contains
 * information about its type, location on the reference genome, and other relevant properties.
 */
public class Feature extends Attributable {

    /**
     * The type of this genomic feature.
     * <p>
     * This field specifies the type of the feature, such as "gene", "exon", or "CDS".
     * The type is defined according to the GFF3 specification and provides information
     * about the biological or functional classification of the feature.
     * <p>
     * For more details, refer to the GFF3 specification:
     * <a href="https://gmod.org/wiki/GFF3">https://gmod.org/wiki/GFF3</a>.
     */
    public final String type;

    /**
     * The name or internal identifier of this genomic feature.
     * <p>
     * This field uniquely identifies the feature within the context of the analysis.
     * It is a final field, meaning its value is immutable once assigned during the
     * construction of the {@link Feature} instance.
     */
    public final String name;

    /**
     * The location of the feature on the reference, i.e., the contig/chromosome/plasmid.
     */
    public final String contig;

    /**
     * 1-based indexed starting position of the feature.
     */
    public final int start;

    /**
     * 1-based indexed end position of the feature.
     */
    public final int end;

    /**
     * Strand of the feature.
     */
    public final char strand;

    /**
     * Unique identifier of this feature.
     */
    public final String _uid;

    /**
     * Represents an allele associated with a genomic feature.
     * <p>
     * An allele is a specific variant of a sequence type that is associated with
     * a genomic feature. This class extends {@link SequenceType} and provides
     * functionality for managing alleles, including their unique identifiers,
     * variants, and attributes.
     */
    public class Allele extends SequenceType {

        /**
         * Constructs a new {@link Allele} instance associated with a genomic feature.
         * <p>
         * This constructor initializes an allele with a unique identifier and a list of variants. Default attributes are set
         * in the {@link Feature#updateAllele} method.
         *
         * @param uid      The unique identifier for the allele.
         * @param variants A list of {@link Tuple} objects representing the variants
         *                 associated with the allele. Each tuple contains:
         *                 <ul>
         *                   <li>The position of the variant.</li>
         *                   <li>The alternate base string of the variant.</li>
         *                 </ul>
         */
        private Allele(String uid, List<Tuple<Integer, String>> variants) {
            super(uid, variants);
            // Add the allele to the feature, if it is not already present.
            Feature.this.alleles.putIfAbsent(uid, this);
        }
    }

    /**
     * Represents a proteoform associated with a genomic feature.
     * <p>
     * A proteoform is a specific variant of a protein sequence that is derived
     * from the genomic feature. This class extends {@link SequenceType} and provides
     * functionality for managing proteoforms, including their unique identifiers,
     * variants, and attributes.
     */
    public class Proteoform extends SequenceType {

        /**
         * Constructs a new {@link Proteoform} instance associated with a genomic feature.
         * <p>
         * This constructor initializes a proteoform with a unique identifier and a list of variants. Default attributes are set
         * in the {@link Feature#updateProteoform} method.
         *
         * @param uid      The unique identifier for the proteoform.
         * @param variants A list of {@link Tuple} objects representing the variants
         *                 associated with the proteoform. Each tuple contains:
         *                 <ul>
         *                   <li>The position of the variant.</li>
         *                   <li>The alternate base string of the variant.</li>
         *                 </ul>
         */
        private Proteoform(String uid, List<Tuple<Integer, String>> variants) {
            super(uid, variants);
            // Add the proteoform to the feature.
            Feature.this.proteoforms.putIfAbsent(uid, this);
        }
    }

    /**
     * Alleles ({@link SequenceType} instances) associated with this feature.
     * <p>
     * This map stores alleles that are associated with the feature. Alleles represent
     * specific sequence variations of the feature. The keys in the map are unique
     * identifiers for the alleles, and the values are the corresponding {@link Allele}
     * instances.
     * <p>
     * Alleles are used to track and manage sequence variations resulting from genomic
     * changes. Each allele is linked to its unique identifier and contains information
     * about its sequence and attributes.
     */
    protected final HashMap<String, Allele> alleles = new HashMap<>();

    /**
     * Proteoforms ({@link SequenceType} instances) associated with this feature.
     * <p>
     * This map stores proteoforms that are associated with the feature. Proteoforms represent
     * specific sequence variants of proteins derived from the feature. The keys in the map
     * are unique identifiers for the proteoforms, and the values are the corresponding
     * {@link Proteoform} instances.
     * <p>
     * Proteoforms are only relevant for coding features and are used to track and manage
     * protein sequence variations resulting from genomic changes.
     */
    protected final HashMap<String, Proteoform> proteoforms = new HashMap<>();

    /**
     * Constructs a new {@link Feature} instance with the specified properties.
     * <p>
     * This constructor initializes a genomic feature with its name, location, strand orientation,
     * type, and unique identifier. The feature's start and end positions are converted to integers
     * to ensure proper indexing. The {@link Attributable} superclass is also initialized.
     *
     * @param name   The name of the feature, used as its internal identifier.
     * @param contig The name of the reference location (e.g., contig, chromosome, plasmid) where the feature is located.
     * @param start  The 1-based indexed starting position of the feature on the reference.
     * @param end    The 1-based indexed end position of the feature on the reference.
     * @param strand The strand orientation of the feature ('+' for forward strand, '-' for reverse strand).
     * @param type   The type of the feature (e.g., coding, non-coding).
     * @param uid    The unique identifier of the feature.
     */
    protected Feature(String name, String contig, Number start, Number end, char strand, String type, String uid) {
        super();
        // TODO: Implement checks for the provided parameters.
        this.name = name;
        this.contig = contig;
        this.type = type;
        this.start = start.intValue();
        this.end = end.intValue();
        this.strand = strand;
        this._uid = uid;
    }

    /**
     * Updates or creates an allele associated with a given contig and sample.
     * <p>
     * This method checks if the specified allele is already associated with the feature.
     * If it is, the existing allele is retrieved. Otherwise, a new allele is created
     * based on the provided variants and contig. The method validates the input parameters,
     * adds the sample occurrence to the allele, and updates the contig with the sequence
     * type occurrences for each variant.
     *
     * @param contig   The {@link Contig} object containing the reference sequence.
     * @param variants A list of {@link Tuple} objects representing the variants
     *                 associated with the allele. Each tuple contains:
     *                 <ul>
     *                   <li>The position of the variant.</li>
     *                   <li>The alternate allele sequence.</li>
     *                 </ul>
     * @param sample   The {@link Sample} object representing the sample associated with this feature.
     * @return The unique identifier (UID) of the updated or created allele.
     */
    protected String updateAllele(Contig contig, ArrayList<Tuple<Integer, String>> variants, Sample sample) {
        // Generate a unique identifier (UID) for the allele based on the variants.
        String uid = IO.md5Hash(SequenceType.variantsAsString(variants));
        Allele allele;

        // Check if the allele already exists in the feature.
        if (this.alleles.containsKey(uid)) {
            allele = getAllele(uid);
        } else {
            // Validate that the contig name matches the feature's contig.
            if (!Objects.equals(contig.name, Feature.this.contig))
                throw new IllegalArgumentException("Contig does not match feature contig.");

            // Validate that the variants map is not empty.
            if (variants.isEmpty())
                throw new IllegalArgumentException("Variants cannot be empty.");

            // Create a new allele if it does not already exist.
            allele = new Allele(uid, variants);

            // Add default attributes for effects and net-shift.
            int lengthVariation = SequenceType.computeLengthVariation(variants);
            int netFrameshift = Math.abs(lengthVariation % 3);
            allele.setAttribute(Constants.$SequenceType_sequenceLengthVariation, String.valueOf(SequenceType.computeLengthVariation(variants)));
            Set<String> effects = contig.getVariantsEffects(variants);
            if (netFrameshift != 0) {
                effects.add(lengthVariation > 0 ? "plus_%d_frameshift".formatted(netFrameshift) : "minus_%d_frameshift".formatted(netFrameshift));
            }
            allele.setAttribute(Constants.$SequenceType_effects, String.join(Constants.COMMA, effects));
        }

        // Add the sample occurrence to the allele.
        allele.addOccurrence(sample.name);

        // Add the allele (occurrence) to the sample.
        sample.setAllele(this.name, allele._uid);

        // Add the allele occurrence to each variant.
        variants.forEach(variant -> contig.getVariantInformation(variant.a, variant.b).addAlleleOccurrence(this.name, allele._uid));

        return uid;
    }

    /**
     * Retrieves an allele associated with this feature by its unique identifier.
     * <p>
     * This method searches for an {@link Allele} in the internal map of alleles using the provided
     * unique identifier (UID). If the UID is not found, the method returns {@code null}.
     *
     * @param uid The unique identifier of the allele to retrieve.
     * @return The {@link Allele} object associated with the given UID, or {@code null} if not found.
     */
    public Allele getAllele(String uid) {
        return this.alleles.getOrDefault(uid, null);
    }

    /**
     * Retrieves all alleles associated with this feature.
     * <p>
     * This method returns a collection of {@link Allele} objects that are associated
     * with this feature. The alleles are stored as values in the internal map of alleles.
     *
     * @return A {@link Collection} of {@link Allele} objects associated with this feature.
     */
    public Collection<Allele> getAlleles() {
        return this.alleles.values();
    }

    /**
     * Updates or creates a proteoform associated with a given allele and contig.
     * <p>
     * This method checks if the specified allele is already associated with a proteoform.
     * If it is, the existing proteoform is retrieved. Otherwise, a new proteoform is created
     * based on the allele's variants and the contig's sequence. The method validates the input
     * parameters, computes the proteoform sequence, and associates the allele with the proteoform.
     *
     * @param contig    The {@link Contig} object containing the reference sequence.
     * @param alleleUid The unique identifier of the allele to update or associate with a proteoform.
     * @throws IOException              If an I/O error occurs during sequence operations.
     * @throws MusialException          If an error occurs during sequence alignment or translation.
     * @throws IllegalArgumentException If the allele does not exist, the contig does not match the feature's contig,
     *                                  the contig lacks a sequence, the feature is not coding, or the variants map is empty.
     */
    protected void updateProteoform(Contig contig, String alleleUid) throws IOException, MusialException {
        // Check if the allele UID is valid.
        Allele allele = getAllele(alleleUid);
        if (allele == null)
            throw new IllegalArgumentException("Allele %s does not exist.".formatted(alleleUid));

        // Validate that the contig matches the feature's contig.
        if (!Objects.equals(contig.name, Feature.this.contig))
            throw new IllegalArgumentException("Contig does not match feature contig.");

        // Validate that the contig has a sequence.
        if (!contig.hasSequence())
            throw new IllegalArgumentException("Contig does not have a sequence.");

        // Validate that the feature is coding.
        if (!Feature.this.isCoding())
            throw new IllegalArgumentException("Feature is not coding.");

        // Access the variants of the allele.
        NavigableMap<Integer, String> variants = allele.variants;
        // Validate that the variants map is not empty.
        if (variants.isEmpty())
            throw new IllegalArgumentException("Variants of an allele can't be empty.");

        // Compute the reference and proteoform sequences.
        String referenceSequence =
                SequenceOperations.translateSequence(contig.getSubsequence(start, end), isReverse());
        String proteoformSequence =
                SequenceOperations.translateSequence(SequenceOperations.integrateVariants(contig, this, variants, true), isReverse());

        /* Generate a unique identifier for the proteoform out of its sequence; if no variants are present, use "synonymous".
         * NOTE: In contrast to alleles, the UID of the proteoform is generated from its sequence to avoid the incorporation of possible
         * alignment errors and to avoid unnecessary calculations for the alignment if the sequence is already known.
         */
        String proteoformUid;
        if (Objects.equals(referenceSequence, proteoformSequence)) {
            proteoformUid = Constants.synonymous;
        } else {
            Proteoform proteoform;
            proteoformUid = IO.md5Hash(proteoformSequence);

            // Check if the proteoform already exists in the feature.
            if (this.proteoforms.containsKey(proteoformUid)) {
                proteoform = getProteoform(proteoformUid);
            } else {
                // Perform global protein sequence alignment.
                Tuple<String, String> alignment = SequenceOperations.globalProteinSequenceAlignment(
                        referenceSequence,
                        proteoformSequence,
                        Math.max(referenceSequence.length(), proteoformSequence.length()), // To avoid more gaps than defined by variant's deletions.
                        6,
                        SequenceOperations.MarginalGaps.FORBID,
                        SequenceOperations.MarginalGaps.PENALIZE,
                        Math.abs(referenceSequence.length() - proteoformSequence.length())
                );

                // Extract canonical amino acid variants from the alignment.
                List<Tuple<Integer, String>> aaVariants = SequenceOperations.getCanonicalVariants(alignment.a, alignment.b).stream()
                        .map(entry -> new Tuple<>(entry.getLeft() + 1, entry.getRight()))
                        .collect(Collectors.toList());
                if (aaVariants.isEmpty()) // This should not happen, but in case it does, throw an exception.
                    throw new IllegalArgumentException("Proteoform has no variants.");

                // Create a new proteoform instance.
                proteoform = new Proteoform(proteoformUid, aaVariants);

                // Add default attributes for effects and net-shift.
                int lengthVariation = SequenceType.computeLengthVariation(aaVariants);
                proteoform.setAttribute(Constants.$SequenceType_sequenceLengthVariation, String.valueOf(lengthVariation));

                Set<String> effects = new HashSet<>();
                String alleleEffects = allele.getAttribute(Constants.$SequenceType_effects);

                if (alleleEffects.contains("frameshift")) effects.add("frameshift_sequence_variation");
                if (aaVariants.stream().anyMatch(variant -> VariantInformation.isInsertion(variant.b))) effects.add("amino_acid_insertion");
                if (aaVariants.stream().anyMatch(variant -> VariantInformation.isDeletion(variant.b))) effects.add("amino_acid_deletion");
                if (proteoform.hasVariant(1) && proteoform.getVariant(1).charAt(0) != referenceSequence.charAt(0))
                    effects.add("start_lost");
                aaVariants.stream()
                        .filter(variant -> variant.b.contains(Constants.stopCodon))
                        .findFirst()
                        .ifPresent(stopCodonVariant -> {
                            int stopCodonPosition = stopCodonVariant.a + stopCodonVariant.b.indexOf(Constants.stopCodon);
                            if (stopCodonPosition <= referenceSequence.length()) {
                                effects.add("stop_gained");
                            } else {
                                effects.add("redundant_inserted_stop_gained");
                            }
                        });
                proteoform.setAttribute(Constants.$SequenceType_effects, String.join(Constants.COMMA, effects));

            }

            // Associate the proteoform with the allele.
            proteoform.addOccurrence(alleleUid);
        }

        // Associate the allele with the proteoform.
        allele.setAttribute(Constants.$Allele_proteoform, proteoformUid);
    }

    /**
     * Retrieves a proteoform associated with this feature by its unique identifier.
     * <p>
     * This method searches for a {@link Proteoform} in the internal map of proteoforms using the provided
     * unique identifier (UID). If the UID is not found, the method returns {@code null}.
     *
     * @param uid The unique identifier of the proteoform to retrieve.
     * @return The {@link Proteoform} object associated with the given UID, or {@code null} if not found.
     */
    public Proteoform getProteoform(String uid) {
        return this.proteoforms.getOrDefault(uid, null);
    }

    /**
     * Retrieves all proteoforms associated with this feature.
     * <p>
     * This method returns a collection of {@link Proteoform} objects that are associated
     * with this feature. The proteoforms are stored as values in the internal map of alleles.
     *
     * @return A {@link Collection} of {@link Proteoform} objects associated with this feature.
     */
    public Collection<Proteoform> getProteoforms() {
        return this.proteoforms.values();
    }

    /**
     * Retrieves a sorted map of child features associated with this feature.
     * <p>
     * This method parses the "children" attribute of the feature, if present, and constructs
     * a sorted map of child features. Each child feature is represented by a key (child type)
     * and a list of tuples, where each tuple contains the start and end positions of the child feature.
     * The map is sorted using a custom comparator based on the order defined in {@link Storage#SO}.
     * If a child type is not found in {@link Storage#SO}, it is assigned the maximum possible value.
     *
     * @return A {@link SortedMap} where the keys are child feature types (e.g., "CDS"),
     * and the values are lists of {@link Tuple} objects representing the start
     * and end positions of the child features.
     */
    public SortedMap<String, List<Tuple<Integer, Integer>>> getChildren() {
        // Initialize a sorted map to store child features, sorted by a custom comparator.
        SortedMap<String, List<Tuple<Integer, Integer>>> children =
                new TreeMap<>(Comparator.comparingInt(k -> Storage.SO.getOrDefault(k, Integer.MAX_VALUE)));

        // Check if the "children" attribute exists for this feature.
        if (hasAttribute("children")) {
            // Split the "children" attribute value into individual child entries using a comma as the delimiter.
            for (String child : getAttribute("children").split(Constants.COMMA)) {

                // Split each child entry into parts using a colon as the delimiter.
                String[] parts = child.split(Constants.COLON);

                // Ensure the child entry has exactly three parts: type, start, and end.
                if (parts.length == 3) {
                    // Add the child feature to the map. If the key (child type) does not exist,
                    // create a new list for it. Then, add a tuple containing the start and end positions.
                    children.computeIfAbsent(parts[0], k -> new ArrayList<>())
                            .add(new Tuple<>(Integer.parseInt(parts[1]), Integer.parseInt(parts[2])));
                }
            }
        }

        // Return the sorted map of child features.
        return children;
    }

    /**
     * Stores the child features of this feature as a serialized string in the "children" attribute.
     * <p>
     * This method serializes the child features into a string format where each child is represented
     * as "type:start:end". Multiple child features are separated by commas. The serialized string
     * is then stored in the "children" attribute of this feature.
     *
     * @param children A {@link SortedMap} where the keys are child feature types (e.g., "CDS"),
     *                 and the values are lists of {@link Tuple} objects representing the start
     *                 and end positions of the child features.
     */
    public void setChildren(SortedMap<String, List<Tuple<Integer, Integer>>> children) {
        // Initialize a StringBuilder to construct the serialized "children" attribute value.
        StringBuilder sb = new StringBuilder();

        // Iterate over each child type and its associated locations.
        children.forEach((key, locations) ->
                // For each location, append the type, start, and end positions to the StringBuilder.
                locations.forEach(location ->
                        sb.append(key).append(Constants.COLON) // Append the child type.
                                .append(location.a).append(Constants.COLON) // Append the start position.
                                .append(location.b).append(Constants.COMMA) // Append the end position and a comma.
                )
        );

        // If the StringBuilder is not empty, remove the trailing comma.
        if (sb.length() > 0) {
            sb.setLength(sb.length() - 1);
        }

        // Add the serialized "children" string as an attribute to this feature.
        setAttribute(Constants.$Feature_children, sb.toString());
    }

    /**
     * Determines if this feature is a coding feature.
     * <p>
     * This method checks whether the feature has child elements of type "CDS" (coding sequence)
     * and ensures that the list of such child elements is not empty. A feature is considered
     * coding if it contains at least one "CDS" child.
     *
     * @return {@code true} if the feature is a coding feature, {@code false} otherwise.
     */
    public boolean isCoding() {
        SortedMap<String, List<Tuple<Integer, Integer>>> children = this.getChildren();
        return children.containsKey("CDS") && !children.get("CDS").isEmpty();
    }

    /**
     * Checks if this feature is on the reverse strand.
     * <p>
     * This method determines whether the strand orientation of the feature is reverse
     * by checking if the strand character is {@code '-'}.
     *
     * @return {@code true} if this feature is on the reverse strand, {@code false} otherwise.
     */
    public boolean isReverse() {
        return this.strand == '-';
    }

    /**
     * Converts this feature into a tab-delimited string representation.
     * <p>
     * This method generates a string containing the contig, type, name, start position,
     * end position, strand, and attributes of the feature. Certain attributes, such as
     * "children", are excluded from the string representation.
     *
     * @return A {@link String} representing the feature in a tab-delimited format.
     */
    public String toString() {
        Set<String> exclude = new HashSet<>(List.of(Constants.$Feature_children));
        return contig + "\t" +
                type + "\t" +
                name + "\t" +
                start + "\t" +
                end + "\t" +
                strand + "\t" +
                attributesAsString(exclude);
    }

    /**
     * Converts this feature and its child features into a GFF3 format string.
     * <p>
     * This method generates a GFF3 representation of the feature, including its attributes
     * and child features. The main feature is represented with its type, location, strand,
     * and attributes. Child features are appended with their respective types, locations,
     * and parent-child relationships.
     *
     * @return A {@link String} containing the GFF3 representation of this feature and its children.
     */
    public String toGffString() {
        StringBuilder contentBuilder = new StringBuilder();

        // Determine GFF3 conforming ID attribute.
        String id = type.contains("gene") ? "gene-%s".formatted(_uid) : "%s-%s".formatted(type, _uid);
        if (!type.contains("gene")) {
            Logging.logWarning("Feature %s type is not a gene; This may conflict with the GFF3 definition.".formatted(name));
        }

        // Append feature information.
        contentBuilder.append(String.join(TAB,
                contig, Musial.softwareName, type, String.valueOf(start), String.valueOf(end),
                DOT, String.valueOf(strand), DOT, "ID=%s".formatted(id)));

        // Append attributes if present.
        if (hasAttributes()) {
            contentBuilder.append(";%s".formatted(attributesAsString(Set.of(
                    "children", "reference_proportion", "sequence_types_disrupted",
                    "sequence_types_modified", "sequence_types_synonymous"))));
        }
        contentBuilder.append(Constants.lineSeparator);

        // Append child features.
        getChildren().forEach((childType, locations) -> locations.forEach(location ->
                contentBuilder.append(String.join(TAB,
                                contig, Musial.softwareName, childType, String.valueOf(location.a), String.valueOf(location.b),
                                DOT, String.valueOf(strand), DOT, getChildIdAndParent(childType, id)))
                        .append(Constants.lineSeparator)
        ));

        return contentBuilder.toString();
    }

    /**
     * Generates a GFF3-compliant ID and Parent attribute string for a child feature.
     * <p>
     * This method constructs the ID and Parent attributes for a child feature based on its type.
     * The format of the attributes depends on the child type:
     * <ul>
     *     <li>For "CDS", the ID is prefixed with "cds-" and the Parent is "transcript-".</li>
     *     <li>For "exon", the ID is prefixed with "exon-" and the Parent is "transcript-".</li>
     *     <li>For RNA types (e.g., "mRNA"), the ID is prefixed with "transcript-" and the Parent is the provided parent ID.</li>
     *     <li>For other types, the ID is prefixed with the child type and the Parent is the provided parent ID.</li>
     * </ul>
     *
     * @param childType The type of the child feature (e.g., "CDS", "exon", "mRNA").
     * @param parentId  The ID of the parent feature.
     * @return A formatted string containing the ID and Parent attributes for the child feature.
     */
    private String getChildIdAndParent(String childType, String parentId) {
        return switch (childType) {
            case "CDS" -> "ID=cds-%s;Parent=transcript-%s".formatted(_uid, _uid);
            case "exon" -> "ID=exon-%s;Parent=transcript-%s".formatted(_uid, _uid);
            default -> childType.contains("RNA")
                    ? "ID=transcript-%s;Parent=%s".formatted(_uid, parentId)
                    : "ID=%s-%s;Parent=%s".formatted(childType, _uid, parentId);
        };
    }

}