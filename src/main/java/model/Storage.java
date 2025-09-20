package model;

import com.google.gson.Gson;
import com.google.gson.TypeAdapter;
import com.google.gson.internal.LinkedTreeMap;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonWriter;
import exceptions.MusialException;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import util.Bio;
import util.IO;
import util.Logging;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * Central component of the MUSIAL model, designed to manage genomic data, including contigs, features, and samples.
 * <p>
 * It provides methods for adding, retrieving, and processing genomic information, as well as handling variant calls and annotations. The
 * class is structured to support efficient storage and manipulation of data, leveraging Java collections and utility classes. Operations
 * executed on storage instances are implemented in separate classes within the {@code op} package.
 * <p>
 * <b>Core Data Structures:</b>
 * <ul>
 * <li>Contigs: Stored in a `Map (String, {@link Contig})`, contigs represent chromosomes or plasmids. Each contig can store its sequence
 * and associated variants.</li>
 * <li>Features: Stored in a `Map (String, {@link Feature})`, features represent genomic elements like genes or mRNA. These are validated
 * and processed using Sequence Ontology (SO) terms.</li>
 * <li>Samples: Stored in a `Map (String, {@link Sample}), samples represent variant calls from distinct biological samples. Metadata and
 * variant calls are associated with each sample.</li>
 * </ul>
 */
public class Storage {

    /**
     * A static map defining the hierarchy levels of Sequence Ontology (SO) terms used in the model.
     * <p>
     * This map associates various SO terms with their respective hierarchy levels, which are used to categorize genomic features. The
     * hierarchy levels are represented as integers, where a lower number indicates a higher-level feature (e.g., "region" at level 0) and a
     * higher number indicates a more specific feature (e.g., "CDS" at level 3).
     * <p>
     * The map includes common SO terms such as "gene", "mRNA", "CDS", and others. It is designed to support the processing and validation
     * of genomic features in the model. Additional terms like UTRs can be optionally added in the future.
     * <p>
     * <b>Example:</b>
     * <ul>
     *   <li>"region" is assigned level 0, representing the highest-level feature.</li>
     *   <li>"gene" and "pseudogene" are assigned level 1, representing primary genomic elements.</li>
     *   <li>"mRNA" and other RNA types are assigned level 2, representing transcripts.</li>
     *   <li>"CDS" and "exon" are assigned level 3, representing coding sequences and exons.</li>
     * </ul>
     */
    public static final Map<String, Integer> SEQUENCE_ONTOLOGY_HIERARCHY = Map.ofEntries(
            Map.entry("region", 0),
            Map.entry("gene", 1),
            Map.entry("pseudogene", 1),
            Map.entry("mRNA", 2),
            Map.entry("tRNA", 2),
            Map.entry("rRNA", 2),
            Map.entry("tmRNA", 2),
            Map.entry("ncRNA", 2),
            Map.entry("SRP_RNA", 2),
            Map.entry("RNase_P_RNA", 2),
            Map.entry("CDS", 3),
            Map.entry("exon", 3)
    );

    /**
     * Static parameters used by this storage.
     * <p>
     * This field holds an instance of {@link Parameters}, which contains the configuration for the storage system. The parameters are
     * immutable and define the behavior of the storage, such as thresholds and exclusions.
     */
    public final Parameters parameters;

    /**
     * The parameters used for configuring the storage of variant data.
     * <p>
     * This record encapsulates various parameters that control the behavior of the storage system, including thresholds for coverage and
     * frequency, as well as exclusions for specific positions and variants. These parameters are immutable once set.
     *
     * @param minimalCoverage  The minimal coverage of a variant call to be accepted. Must be greater than or equal to 0. This ensures that
     *                         only variant calls with sufficient read depth are considered.
     * @param minimalFrequency The minimal frequency of a variant call to be accepted. Must be between 0.0 and 1.0, inclusive. This
     *                         parameter filters out low-frequency variants that may be due to sequencing errors.
     * @param maskFiltered     Whether to mask filtered calls as ambiguous bases (N) or ignore them.
     * @param skipAnnotation   Whether to skip the SnpEff annotation process. If true, the annotation step is bypassed, which can save time
     *                         if annotation is not required.
     * @param skipTyping       Whether to skip proteoform inference. If true, the inference of proteoforms (protein isoforms) is not
     *                         performed, which can be useful for non-coding regions.
     * @param masked           A map associating contig names with sets of positions to exclude from storage. Cannot be null but can be
     *                         empty. This allows specific genomic positions to be ignored during analysis.
     */
    public record Parameters(
            int minimalCoverage, // Minimum read depth required for a variant call to be accepted.
            double minimalFrequency, // Minimum allele frequency required for a variant call to be accepted.
            boolean maskFiltered, // Whether to mask filtered variants in the analysis.
            boolean skipAnnotation, // Flag to determine whether SnpEff annotation should be skipped.
            boolean skipTyping, // Flag to determine whether proteoform inference should be skipped.
            Map<String, Set<Integer>> masked // Map of contig names to sets of positions to exclude from analysis.
    ) {
        // Compact constructor with validation logic omitted for simplicity.

        /**
         * Whether {@code position} is excluded on {@code contig}.
         *
         * @param contig   Contig (id) to check for exclusion.
         * @param position Position to check for exclusion.
         * @return True if {@code position} on {@code contig} is excluded from analysis.
         */
        public boolean isPositionMasked(String contig, int position) {
            return masked.containsKey(contig) && masked.get(contig).contains(position);
        }

    }

    /**
     * Transient accessor to (indexed) reference sequences.
     */
    private transient IndexedFastaSequenceFile reference = null;

    /**
     * Contigs/chromosomes/plasmids of the reference sequence.
     */
    private final Map<String, Contig> contigs;

    /**
     * Genomic features.
     */
    private final Map<String, Feature> features;

    /**
     * Individual samples, i.e., variant calls from one distinct biological sample.
     */
    private final Map<String, Sample> samples;

    /**
     * Constructs a new {@link Storage} instance with the specified parameters.
     * <p>
     * This constructor initializes the {@link Storage} object with the provided configuration parameters. It also initializes empty
     * containers for contigs, features, and samples using {@link LinkedTreeMap}, ensuring that the data is stored in a sorted and efficient
     * manner.
     *
     * @param parameters The {@link Parameters} object containing the configuration for the storage. This includes thresholds, exclusions,
     *                   and other settings for managing genomic data.
     */
    public Storage(Parameters parameters) {
        this.parameters = parameters;
        // Todo: Initialize with capacity estimates.
        this.contigs = new LinkedHashMap<>(); // Initialize empty contig container.
        this.features = new LinkedHashMap<>(); // Initialize empty feature container.
        this.samples = new LinkedHashMap<>(); // Initialize empty sample container.
    }

    /**
     * Checks if the reference sequence is set for the storage.
     * <p>
     * This method determines whether the {@link #reference} field has been initialized with a non-null value. The {@link #reference} field
     * is a transient accessor to the indexed reference sequences used in the storage. If the reference is set, it indicates that the
     * storage has access to the reference sequence for further operations.
     *
     * @return {@code true} if the {@link #reference} field is non-null, indicating that the reference sequence is set; {@code false}
     * otherwise.
     */
    public boolean hasReference() {
        return Objects.nonNull(this.reference);
    }

    /**
     * Sets the reference sequence for the storage and populates contigs based on the reference.
     * <p>
     * This method assigns the provided {@link IndexedFastaSequenceFile} as the reference sequence for the storage. It clears any previously
     * stored contigs and repopulates them based on the sequences available in the reference.
     * <p>
     * The method iterates through all sequences in the reference file, adding each sequence as a contig to the storage if it does not
     * already exist. After processing all sequences, the reference file is reset to its initial state.
     * <p>
     * This method should only be called when creating new instances of {@link Storage} and not during deserialization!
     *
     * @param indexedFastaSequenceFile The {@link IndexedFastaSequenceFile} instance representing the reference sequence. Must not be null.
     * @throws IOException    If an error occurs while adding contigs to the storage.
     * @throws AssertionError If the provided {@code indexedFastaSequenceFile} is null.
     */
    public void setReference(IndexedFastaSequenceFile indexedFastaSequenceFile) throws IOException {
        // Ensure the provided reference file is not null.
        assert Objects.nonNull(indexedFastaSequenceFile) : "indexedFastaSequenceFile cannot be null.";

        // Assign the reference file to the storage and clear existing contigs.
        this.reference = indexedFastaSequenceFile;
        this.contigs.clear();

        // Iterate through all sequences in the reference file.
        ReferenceSequence referenceSequence = this.reference.nextSequence();
        while (Objects.nonNull(referenceSequence)) {
            // Add the sequence as a contig if it does not already exist in the storage.
            if (!this.hasContig(referenceSequence.getName())) {
                this.addContig(referenceSequence.getName(), referenceSequence.getBaseString());
            }
            referenceSequence = this.reference.nextSequence();
        }

        // Reset the reference file to its initial state.
        this.reference.reset();
    }

    /**
     * Checks if a contig is present in the storage by its unique identifier.
     * <p>
     * This method verifies whether a contig, identified by the given id, exists in the storage's contig map. It is useful for determining
     * the presence of a specific contig before performing operations on it.
     *
     * @param identifier The unique id of the contig to check. This id typically represents the name of the chromosome or plasmid.
     * @return {@code true} if the contig is present in the storage; {@code false} otherwise.
     */
    public boolean hasContig(String identifier) {
        return this.contigs.containsKey(identifier);
    }

    /**
     * Adds a contig to the storage under the specified id and sequence.
     * <p>
     * This method is responsible for adding a contig (chromosome or plasmid) to the storage. It ensures that the contig is only added if it
     * does not already exist in the storage. The contig is represented by its name (id) and its sequence.
     * <p>
     * The sequence is compressed using GZIP for efficient storage. If the provided sequence is null or empty, the method assigns an empty
     * string as the compressed sequence and sets the length of the sequence to 0. This ensures that the contig is still added to the
     * storage with its attributes, even if no sequence data is available.
     *
     * @param identifier The unique identifier (id) of the contig to add. This typically represents the name of the chromosome or plasmid.
     * @param sequence   The nucleotide sequence of the contig. This can be null or empty, in which case an empty sequence is stored.
     * @throws IOException If an error occurs during the compression of the sequence data.
     */
    public void addContig(String identifier, String sequence) throws IOException {
        this.contigs.putIfAbsent(identifier, new Contig(identifier, sequence));
    }

    /**
     * Retrieves a contig by its unique identifier.
     * <p>
     * This method looks up a contig in the storage using its unique id. If the contig is found, it returns the corresponding {@link Contig}
     * object. If no contig with the specified id exists, the method returns {@code null}.
     *
     * @param identifier The unique identifier of the contig to retrieve. This id typically represents the name of the chromosome or
     *                   plasmid.
     * @return The {@link Contig} object associated with the specified id, or {@code null} if no such contig exists in the storage.
     */
    public Contig getContig(String identifier) {
        return this.contigs.getOrDefault(identifier, null);
    }

    /**
     * Retrieves an unmodifiable collection view of the contigs stored in the storage.
     * <p>
     * This method provides a read-only view of the contigs stored in the storage. The returned collection reflects the current state of the
     * contigs map but cannot be modified directly. This ensures that the integrity of the underlying data structure is maintained.
     * <p>
     * This method is useful for accessing all contigs in the storage without allowing external modifications.
     *
     * @return An unmodifiable collection of {@link Contig} objects stored in the storage.
     */
    public Collection<Contig> getContigs() {
        return Collections.unmodifiableCollection(this.contigs.values());
    }

    /**
     * Query whether a feature is stored in this instance by its id.
     *
     * @param identifier The id of the feature.
     * @return True if a feature is stored for {@code id}.
     */
    public boolean hasFeature(String identifier) {
        return this.features.containsKey(identifier);
    }

    /**
     * Adds feature information from a {@link FeatureI} object to the storage.
     * <p>
     * This method extracts the necessary details from the provided {@link FeatureI} object, such as the parent contig, start and end
     * positions, strand, and type, and delegates the processing to the overloaded
     * {@link #addFeature(String, String, Number, Number, char, String, Map)} method.
     *
     * @param featureI   The {@link FeatureI} object containing the feature information to transfer.
     * @param name       The name of the feature.
     * @param attributes A map of attributes associated with the feature.
     * @throws MusialException If an error occurs while adding the feature to the storage.
     */
    public void addFeature(FeatureI featureI, String name, Map<String, String> attributes) throws MusialException {
        addFeature(name, featureI.seqname(), featureI.location().bioStart(), featureI.location().bioEnd(),
                featureI.location().bioStrand(),
                featureI.type(), attributes);
    }

    /**
     * Adds feature information to the storage.
     * <p>
     * This method processes and validates the provided feature information, including its type, location, and attributes. It checks if the
     * feature type is supported by the Sequence Ontology (SO) map and determines a unique identifier (UID) for the feature. If a feature
     * with the same UID already exists, it validates compatibility with the parent feature and updates the "children" attribute if
     * applicable. Otherwise, it creates a new feature and adds it to the storage. Processed attributes are removed from the attributes map,
     * and the remaining attributes are extended for the feature.
     * <p>
     * The method performs the following steps:
     * <ul>
     *   <li>Validates the feature type against the Sequence Ontology (SO) map.</li>
     *   <li>Validates the feature's location and compatibility with stored contigs.</li>
     *   <li>Determines a unique identifier (UID) for the feature based on its attributes.</li>
     *   <li>Checks if a feature with the same UID already exists and updates or creates the feature accordingly.</li>
     *   <li>Removes processed attributes and extends the feature's attributes with the remaining ones.</li>
     * </ul>
     *
     * @param name             The name of the feature.
     * @param contigIdentifier The chromosome where the feature is located matching the identifier of a stored contig.
     * @param start            The start position of the feature.
     * @param end              The end position of the feature.
     * @param strand           The strand of the feature ('+' or '-').
     * @param type             The type of the feature (e.g., "gene", "mRNA").
     * @param attributes       A map of attributes associated with the feature. @throws MusialException If an error occurs while adding the
     *                         feature to the storage.
     */
    public void addFeature(String name, String contigIdentifier, Number start, Number end, char strand, String type,
                           Map<String, String> attributes) throws MusialException {
        // Validate the feature type against the Sequence Ontology (SO) map.
        if (!SEQUENCE_ONTOLOGY_HIERARCHY.containsKey(type)) {
            Logging.logWarningOnce(
                    "INVALID_FEATURE_TYPE_%s".formatted(type),
                    "Features of type %s are currently not supported and will be ignored (%s).".formatted(type, name)
            );
            return;
        }

        // Validate the feature location.
        if ((int) start >= (int) end || (int) start < 1) {
            Logging.logWarning("Feature %s has an invalid location (%s:g.%d_%d) and will be ignored."
                    .formatted(name, contigIdentifier, (int) start, (int) end));
            return;
        }

        // Validate the compatibility with stored contigs.
        if (!this.hasContig(contigIdentifier)) {
            Logging.logWarning("Feature %s specifies an unknown parent locus (%s) and will be ignored."
                    .formatted(name, contigIdentifier));
            return;
        } else {
            int length = this.getContig(contigIdentifier).getSequenceLength();
            if (length != 0 && length < (int) end) {
                Logging.logWarning(("Feature %s (%s:g.%d_%d) exceeds the length of its parent (%d) and will be ignored.")
                        .formatted(name, contigIdentifier, (int) start, (int) end, length));
                return;
            }
        }

        // Determine the identifier for the feature.
        final String identifier;
        if (attributes.containsKey("Parent")) {
            identifier = attributes.get("Parent").matches("^.*-.*$") ? attributes.get("Parent").split("-")[1] : attributes.get("Parent");
        } else if (attributes.containsKey("ID")) {
            identifier = attributes.get("ID").matches("^.*-.*$") ? attributes.get("ID").split("-")[1] : attributes.get("ID");
        } else {
            identifier = attributes.containsKey("locus_tag") ? attributes.get("locus_tag") : "%s:g.%d_%d=".formatted(contigIdentifier,
                    (int) start,
                    (int) end);
        }

        // Check if a feature with the same identifier already exists.
        // Note: If a feature with the same identifier exists, it is ALWAYS assumed the parent feature (due to the sorting of GFF3 files).
        Optional<Feature> optional = this.features.values().stream().filter(feature -> feature._id.equals(identifier)).findFirst();
        Feature feature;

        if (optional.isPresent()) {
            feature = optional.get();

            // Validate compatibility with the parent feature if the "Parent" attribute is present.
            if (attributes.containsKey("Parent") && attributes.get("Parent").matches("^.*-%s$".formatted(identifier))) {
                if ((int) start < feature.start || (int) end > feature.end || !feature.contig.equals(contigIdentifier) || feature.strand != strand) {
                    Logging.logWarning(("Feature %s (identifier %s) has an incompatible location with its parent " +
                            "feature %s.")
                            .formatted(feature.name, feature._id, name));
                    return;
                }

                // Add or extend sub-features of the feature.
                feature.addSubFeature(type, (int) start, (int) end);
            } else {
                Logging.logWarning("Feature %s (identifier %s) already exists, but is no valid parent of feature %s."
                        .formatted(feature.name, feature._id, name));
                return;
            }
        } else {
            // Create a new feature if it does not already exist.
            feature = new Feature(name, contigIdentifier, start, end, strand, type, identifier);
            this.features.put(identifier, feature);
        }

        // Remove processed attributes from the attributes map.
        // Note: These attributes are already stored in the structure of the storage or not relevant. May be extended.
        attributes.remove("ID");
        attributes.remove("Parent");
        attributes.remove("Name");
        attributes.remove("old_locus_tag");
        attributes.remove("gbkey");
        // Extend the feature's attributes with the remaining attributes.
        feature.extendAttributes(attributes);
    }

    /**
     * Replaces an existing feature in the storage with a new feature.
     * <p>
     * This method checks if the feature to be replaced exists in the storage. If the feature exists, it is replaced with the provided
     * replacement feature. If the feature does not exist, a warning is logged, and no changes are made.
     * <p>
     * <b>This method is and should only be used for updating features during feature validation.</b>
     *
     * @param stored      The {@link Feature} object to be replaced. This feature must already exist in the storage.
     * @param replacement The {@link Feature} object to replace the existing feature with.
     * @return The replaced {@link Feature} object if the replacement was successful; {@code null} if the feature to be replaced does not
     * exist.
     */
    public Feature replaceFeature(Feature stored, Feature replacement) {
        if (!this.features.containsKey(stored._id)) {
            Logging.logWarning("Feature %s (identifier %s) cannot be replaced as it does not exist."
                    .formatted(stored.name, stored._id));
            return null;
        }
        return this.features.put(stored._id, replacement);
    }

    /**
     * Retrieves a feature from the storage by its unique identifier.
     * <p>
     * This method looks up a feature in the storage using its unique id. If the feature is found, it returns the corresponding
     * {@link Feature} object. If no feature with the specified id exists, the method returns {@code null}.
     * <p>
     * Todo: It may be relevant for some applications to retrieve features by other attributes, such as "locus_tag" or "ID".
     *
     * @param identifier The unique identifier of the feature to retrieve. This id typically represents the name or id of the genomic
     *                   feature.
     * @return The {@link Feature} object associated with the specified id, or {@code null} if no such feature exists in the storage.
     */
    public Feature getFeature(String identifier) {
        return this.features.getOrDefault(identifier, null);
    }

    /**
     * Retrieves an unmodifiable collection view of the features stored in the storage.
     * <p>
     * This method provides a read-only view of the features stored in the storage. The returned collection reflects the current state of
     * the features map but cannot be modified directly. This ensures that the integrity of the underlying data structure is maintained.
     *
     * @return An unmodifiable collection of {@link Feature} objects stored in the storage.
     */
    public Collection<Feature> getFeatures() {
        return Collections.unmodifiableCollection(this.features.values());
    }

    /**
     * Removes a feature from the storage by its unique identifier.
     * <p>
     * This method deletes the feature associated with the given identifier from the `features` map. It is useful for managing the storage
     * by allowing the removal of specific genomic features when they are no longer needed or relevant.
     *
     * @param identifier The unique identifier of the feature to remove. This id typically represents the name or id of the genomic
     *                   feature.
     */
    public void removeFeature(String identifier) {
        this.features.remove(identifier);
    }

    /**
     * Checks if a sample is present in the storage by its unique identifier.
     * <p>
     * This method verifies whether a sample, identified by the given id, exists in the storage's sample map. It is useful for determining
     * the presence of a specific sample before performing operations on it.
     *
     * @param identifier The unique id of the sample to check. This id typically represents the name or identifier of the biological
     *                   sample.
     * @return {@code true} if the sample is present in the storage; {@code false} otherwise.
     */
    public boolean hasSample(String identifier) {
        return this.samples.containsKey(identifier);
    }

    /**
     * Adds a sample to the storage if it does not already exist.
     * <p>
     * This method checks whether a sample with the given identifier is already present in the storage. If the sample does not exist, it
     * creates a new {@link Sample} object, associates it with the storage, and adds it to the `samples` map. The new {@link Sample} object
     * is initialized with the current number of contigs and features in the storage. If the sample already exists, the method does
     * nothing.
     *
     * @param sampleIdentifier The unique identifier of the sample to add. This typically represents the name or ID of the biological
     *                         sample.
     */
    public void addSample(String sampleIdentifier) {
        // Create a new Sample object with the given identifier, initialized with the current number of contigs and features.
        Sample sample = new Sample(sampleIdentifier, features.size());

        // Add the sample to the samples map if it does not already exist.
        samples.putIfAbsent(sample._id, sample);
    }

    /**
     * Retrieves a sample from the storage by its unique identifier.
     * <p>
     * This method looks up a sample in the storage using its unique identifier. If the sample is found, it returns the corresponding
     * {@link Sample} object. If no sample with the specified identifier exists, the method returns {@code null}.
     *
     * @param identifier The unique identifier of the sample to retrieve. This identifier typically represents the name or ID of the
     *                   biological sample.
     * @return The {@link Sample} object associated with the specified identifier, or {@code null} if no such sample exists in the storage.
     */
    public Sample getSample(String identifier) {
        return this.samples.get(identifier);
    }

    /**
     * Retrieves an unmodifiable collection view of the samples stored in the storage.
     * <p>
     * This method provides a read-only view of the samples stored in the storage. The returned collection reflects the current state of the
     * samples map but cannot be modified directly. This ensures that the integrity of the underlying data structure is maintained.
     * <p>
     * This method is useful for accessing all samples in the storage without allowing external modifications.
     *
     * @return An unmodifiable collection of {@link Sample} objects stored in the storage.
     */
    public Collection<Sample> getSamples() {
        return Collections.unmodifiableCollection(this.samples.values());
    }

    /**
     * Adds a variant to the specified contig in the storage.
     * <p>
     * This method ensures that the variant is in a canonical padded format. If the variant is not canonical, an
     * {@link IllegalArgumentException} is thrown. If the variant does not already exist in the contig, it is created and added. The method
     * also associates the variant with the specified sample and caches the variant for further processing.
     *
     * @param contig           The {@link Contig} object to which the variant belongs. Represents the genomic region where the variant is
     *                         located.
     * @param sampleIdentifier The unique identifier of the sample associated with the variant. Used to track the sample's variant calls.
     * @param position         The position of the variant within the contig. Represents the genomic coordinate of the variant.
     * @param reference        The reference allele of the variant. This is the expected sequence at the given position.
     * @param alternative      The alternative allele of the variant. This is the observed sequence differing from the reference.
     * @param variantCalls     A set of {@link VariantCall} objects representing the variant calls associated with the sample.
     * @throws IllegalArgumentException If the variant is not in a canonical padded format. Ensures data consistency and correctness.
     */
    public void addVariant(Contig contig, String sampleIdentifier, int position, String reference, String alternative,
                           Set<VariantCall> variantCalls) {
        // Validate that the variant is in a canonical padded format.
        if (!Bio.isPaddedCanonical(reference, alternative)) {
            throw new IllegalArgumentException("Failed to add non-canonical variant %s > %s at position %d to contig %s."
                    .formatted(reference, alternative, position, contig._id));
        }

        // Retrieve the variant from the contig, or create it if it does not exist.
        Variant variant = contig.getVariant(position, alternative);
        if (variant == null) {
            variant = new Variant(position, reference, alternative);
            contig.addVariant(variant);
        }

        // Associate the variant with the sample and cache it for further processing.
        variant.addRelation(sampleIdentifier, variantCalls);
    }

    /**
     * Calculates the total number of variants across all contigs in the storage.
     * <p>
     * This method iterates through all contigs stored in the {@code contigs} map and sums up the variant counts for each contig. The
     * variant count for each contig is retrieved using the {@link Contig#getVariantsCount()} method.
     *
     * @return The total number of variants across all contigs.
     */
    public int getVariantsCount() {
        return (int) contigs.values().stream().mapToLong(Contig::getVariantsCount).sum();
    }

    /**
     * Calculates the total number of active variants across all contigs in the storage.
     * <p>
     * This method iterates through all contigs stored in the {@code contigs} map and sums up the active variant counts for each contig. The
     * active variant count for each contig is retrieved using the {@link Contig#getActiveVariantsCount()} method.
     *
     * @return The total number of active variants across all contigs.
     */
    public int getActiveVariantsCount() {
        return (int) contigs.values().stream().mapToLong(Contig::getActiveVariantsCount).sum();
    }

    /**
     * Creates a custom {@link TypeAdapter} for the {@link Storage} class.
     * <p>
     * This method defines a custom {@link TypeAdapter} to handle the serialization and deserialization of {@link Storage} objects. The
     * adapter uses Gson's default adapter for most operations but adds custom behavior during deserialization to initialize the transient
     * {@link #reference} field.
     *
     * @return A {@link TypeAdapter} for the {@link Storage} class.
     */
    public static TypeAdapter<Storage> typeAdapter() {

        return new TypeAdapter<>() {

            // Default adapter for Contig objects provided by Gson
            final TypeAdapter<Storage> defaultAdapter = new Gson().getAdapter(Storage.class);

            /**
             * Serializes a {@link Storage} object to JSON.
             * <p>
             * This method delegates the serialization process to the default adapter.
             *
             * @param out   The {@link JsonWriter} to write the JSON output.
             * @param value The {@link Storage} object to serialize.
             * @throws IOException If an I/O error occurs during writing.
             */
            @Override
            public void write(JsonWriter out, Storage value) throws IOException {
                defaultAdapter.write(out, value);
            }

            /**
             * Deserializes a {@link Storage} object from JSON.
             * <p>
             * This method delegates the deserialization process to the default adapter and then initializes
             * the transient {@code cache} field to ensure the {@link Storage} object is fully functional.
             *
             * @param in The {@link JsonReader} to read the JSON input.
             * @return The deserialized {@link Storage} object.
             * @throws IOException If an I/O error occurs during reading.
             */
            @Override
            public Storage read(JsonReader in) throws IOException {
                Storage storage = defaultAdapter.read(in); // Deserialize using the default adapter

                try { // Build IndexedFastaSequenceFile from contigs if they have non-empty sequences
                    if (!storage.contigs.isEmpty()) {
                        List<String> fastaEntries = new ArrayList<>();
                        for (Contig contig : storage.contigs.values()) {
                            // Extract sequence of contig.
                            if (contig.hasSequence()) {
                                fastaEntries.add(">" + contig._id);
                                fastaEntries.add(contig.getSequence());
                            }
                            // Set transient sequence cache for each contig.
                            contig.sequenceCache = new HashMap<>();
                            // Set all variants to known (not active) after deserialization.
                            contig.getAllVariants().forEach(variant -> variant.active = false);
                        }
                        if (!fastaEntries.isEmpty()) {
                            Path tempFasta = Files.createTempFile(IO.md5Hash(Logging.getTimestamp()), ".fasta");
                            tempFasta.toFile().deleteOnExit();
                            Files.writeString(tempFasta, String.join("\n", fastaEntries));
                            storage.reference = new IndexedFastaSequenceFile(tempFasta,
                                    FastaSequenceIndexCreator.buildFromFasta(tempFasta));
                        } else {
                            storage.reference = null;
                        }
                    } else {
                        storage.reference = null;
                    }
                } catch (Exception e) {
                    throw new IOException("Failed to build IndexedFastaSequenceFile from contigs.", e);
                }

                // Initialize the transient sequence cache for each contig in the storage.
                storage.contigs.values().forEach(contig -> contig.sequenceCache = new HashMap<>());

                return storage;
            }
        };
    }

}