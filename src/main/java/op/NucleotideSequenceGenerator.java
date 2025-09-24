package op;

import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import model.Contig;
import model.Feature;
import model.Storage;
import model.Variant;
import uk.co.omegaprime.btreemap.BTreeMap;
import util.Bio;
import util.Constants;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * The {@code NucleotideSequenceGenerator} class is responsible for generating nucleotide sequences based on genomic data contained in a
 * {@link Storage} instance.
 * <p>
 * This class is initialized with a fixed {@link Contig} and optionally a {@link Feature}. It can be restricted to specific sample
 * identifiers and automatically handles the option specified by {@link model.Storage.Parameters#maskFiltered}.
 */
public class NucleotideSequenceGenerator {

    /**
     * The storage instance containing genomic data.
     */
    final Storage storage;

    /**
     * The contig associated with the sequence generation.
     */
    final Contig contig;

    /**
     * The feature associated with the sequence generation, if any.
     */
    final Feature feature;

    /**
     * A cache for storing previously generated sequences.
     */
    HashMap<String, String> cache;

    /**
     * Flag indicating whether generated sequences contain conserved sites.
     */
    boolean conserved;

    /**
     * Flag indicating whether generated sequences are aligned.
     */
    boolean aligned;

    /**
     * The interval (from, to) for sequence generation with respect to this {@link #contig}.
     */
    Tuple<Integer, Integer> interval;

    /**
     * The nucleotide context for sequence generation.
     */
    NavigableMap<Integer, Bio.ReferenceContext> context;

    /**
     * Sample identifiers to which the sequence generation is restricted.
     */
    final Set<String> sampleIdentifiers;

    /**
     * Constructs a new instance of the NucleotideSequenceGenerator class.
     * <p>
     * This constructor initializes the generator for a given contig and optional sample identifiers. It sets the interval to cover the
     * entire contig and determines whether the sequence generation should include conserved sites and/or align sequences.
     * <p>
     * If sample identifiers are provided, the sequence generation will be restricted to those samples.
     * <p>
     * The constructor validates the contig to ensure it exists in the storage and meets the requirements for conserved sequence generation,
     * if applicable.
     *
     * @param storage           The storage object containing genomic data.
     * @param contig            The contig associated with the sequence generation.
     * @param conserved         A flag indicating whether the sequence generation includes conserved sites.
     * @param aligned           A flag indicating whether the generated sequences are aligned.
     * @param sampleIdentifiers Optional sample identifiers to restrict the sequence generation.
     * @throws IOException     If an error occurs during sequence retrieval.
     * @throws MusialException If an error occurs during context generation.
     */
    public NucleotideSequenceGenerator(Storage storage, Contig contig, boolean conserved, boolean aligned, String... sampleIdentifiers)
            throws IOException, MusialException {
        // Validate the contig to ensure it exists in the storage and meets the requirements.
        validateContig(storage, contig, conserved);

        // Initialize instance variables.
        this.storage = storage;
        this.contig = contig;
        this.feature = null;
        this.cache = null;
        this.conserved = conserved;
        this.aligned = aligned;
        this.interval = new Tuple<>(1, contig.getSequenceLength());
        this.sampleIdentifiers = Arrays.stream(sampleIdentifiers).collect(Collectors.toSet());

        // Generate the nucleotide context based on the conserved flag.
        if (conserved) generateConservedContext();
        else generateContext();
    }

    /**
     * Constructs a new instance of the NucleotideSequenceGenerator class with a specified interval.
     * <p>
     * This constructor initializes the generator for a given contig and interval, along with optional sample identifiers. It sets the
     * interval to the specified range and determines whether the sequence generation should include conserved sites and/or align
     * sequences.
     * <p>
     * If sample identifiers are provided, the sequence generation will be restricted to those samples.
     * <p>
     * The constructor validates the contig to ensure it exists in the storage and meets the requirements for conserved sequence generation,
     * if applicable.
     *
     * @param storage           The storage object containing genomic data.
     * @param contig            The contig associated with the sequence generation.
     * @param from              The start position of the interval for sequence generation.
     * @param to                The end position of the interval for sequence generation.
     * @param conserved         A flag indicating whether the sequence generation includes conserved sites.
     * @param aligned           A flag indicating whether the generated sequences are aligned.
     * @param sampleIdentifiers Optional sample identifiers to restrict the sequence generation.
     * @throws IOException     If an error occurs during sequence retrieval.
     * @throws MusialException If an error occurs during context generation.
     */
    public NucleotideSequenceGenerator(Storage storage, Contig contig, int from, int to, boolean conserved, boolean aligned,
                                       String... sampleIdentifiers) throws IOException, MusialException {
        // Validate the contig to ensure it exists in the storage and meets the requirements.
        validateContig(storage, contig, conserved);

        // Initialize instance variables.
        this.storage = storage;
        this.contig = contig;
        this.feature = null;
        this.cache = null;
        this.conserved = conserved;
        this.aligned = aligned;
        this.interval = new Tuple<>(from, to);
        this.sampleIdentifiers = Arrays.stream(sampleIdentifiers).collect(Collectors.toSet());

        // Generate the nucleotide context based on the conserved flag.
        if (conserved) generateConservedContext();
        else generateContext();
    }

    /**
     * Constructs a new instance of the NucleotideSequenceGenerator class with a specified feature.
     * <p>
     * This constructor initializes the generator for a given contig and feature, along with optional sample identifiers. It sets the
     * interval to the range defined by the feature and determines whether the sequence generation should include conserved sites and/or
     * align sequences.
     * <p>
     * If sample identifiers are provided, the sequence generation will be restricted to those samples.
     * <p>
     * The constructor validates the contig and feature to ensure they exist in the storage and meet the requirements for conserved sequence
     * generation, if applicable.
     *
     * @param storage           The storage object containing genomic data.
     * @param contig            The contig associated with the sequence generation.
     * @param feature           The feature associated with the sequence generation.
     * @param conserved         A flag indicating whether the sequence generation includes conserved sites.
     * @param aligned           A flag indicating whether the generated sequences are aligned.
     * @param sampleIdentifiers Optional sample identifiers to restrict the sequence generation.
     * @throws IOException     If an error occurs during sequence retrieval.
     * @throws MusialException If an error occurs during context generation or feature validation.
     */
    public NucleotideSequenceGenerator(Storage storage, Contig contig, Feature feature, boolean conserved, boolean aligned,
                                       String... sampleIdentifiers) throws IOException, MusialException {
        // Validate the contig to ensure it exists in the storage and meets the requirements.
        validateContig(storage, contig, conserved);

        // Validate the feature to ensure it exists in the storage.
        validateFeature(storage, feature);

        // Initialize instance variables.
        this.storage = storage;
        this.contig = contig;
        this.feature = feature;
        this.cache = new HashMap<>(feature.getAlleleCount() + 1);
        this.conserved = conserved;
        this.aligned = aligned;
        this.interval = new Tuple<>(feature.start, feature.end);
        this.sampleIdentifiers = Arrays.stream(sampleIdentifiers).collect(Collectors.toSet());

        // Generate the nucleotide context based on the conserved flag.
        if (conserved) generateConservedContext();
        else generateContext();
    }

    /**
     * Retrieves the nucleotide sequence for a given sample identifier.
     * <p>
     * This method determines whether the sequence should be generated based on a feature or a contig. If a feature is defined, the sequence
     * is generated using the `fromFeature` method. Otherwise, the sequence is generated using the `fromContig` method.
     *
     * @param sampleIdentifier The identifier of the sample for which the sequence is generated.
     * @return A {@link String} representing the nucleotide sequence for the given sample.
     */
    public String getSequence(String sampleIdentifier) throws MusialException {
        return (feature != null) ? fromFeature(sampleIdentifier) : fromContig(sampleIdentifier);
    }

    /**
     * Generates a nucleotide sequence for a given sample identifier based on the associated feature.
     * <p>
     * This method retrieves the allele identifier related to the specified sample and feature. If the sequence for the allele is already
     * cached, it is returned directly. Otherwise, the sequence is generated by integrating variants associated with the allele into the
     * reference context. The generated sequence is then cached for future use.
     *
     * @param sampleIdentifier The identifier of the sample for which the sequence is generated.
     * @return A {@link String} representing the nucleotide sequence for the given sample and feature.
     * @throws MusialException          If {@link #feature} is {@code null} or from {@link Bio#integrateVariants}.
     * @throws IllegalArgumentException If the sample identifier is invalid.
     */
    private String fromFeature(String sampleIdentifier) throws MusialException {
        // Ensure the feature is defined.
        if (feature == null) {
            throw new MusialException("Cannot generate feature based sequence without a defined feature.");
        }

        // Validate the provided sample identifier.
        validateSample(sampleIdentifier);

        // Retrieve the allele identifier associated with the sample and feature.
        String alleleIdentifier = storage.getSample(sampleIdentifier).getRelatedAllele(feature._id);

        // Return the cached sequence if it is already available.
        if (cache.containsKey(alleleIdentifier)) {
            return cache.get(alleleIdentifier);
        }

        // Determine the sequence based on the allele identifier.
        String sequence = alleleIdentifier.equals(Constants.REFERENCE)
                ? Bio.integrateVariants(context, Collections.emptyMap(), !aligned)
                : Bio.integrateVariants(context, feature.getAllele(alleleIdentifier).getVariants(), !aligned);

        // Cache the generated sequence and return it.
        cache.put(alleleIdentifier, sequence);
        return sequence;
    }

    /**
     * Generates a nucleotide sequence for a given sample identifier by integrating variants.
     * <p>
     * This method retrieves the variants associated with the specified sample within the defined interval and integrates them into the
     * reference sequence. If no variants are found, the reference sequence is returned as is. Variants are filtered based on the sample
     * identifier and storage parameters.
     *
     * @param sampleIdentifier The identifier of the sample for which the sequence is generated.
     * @return A {@link String} representing the nucleotide sequence for the given sample.
     * @throws MusialException From {@link Bio#integrateVariants}.
     */
    private String fromContig(String sampleIdentifier) throws MusialException {
        // Validate the provided sample identifier.
        validateSample(sampleIdentifier);

        // Retrieve variants for the given sample within the interval.
        List<Variant> variants = contig.getVariantsOfSamplesWithin(interval.a, interval.b, Collections.singleton(sampleIdentifier));

        // If no variants are found, return the integrated reference sequence.
        if (variants.isEmpty()) {
            return Bio.integrateVariants(context, Collections.emptyMap(), !aligned);
        }

        // Map variants to their positions and alternative alleles.
        Map<Integer, String> variantsMap = new HashMap<>();
        for (Variant variant : variants) {
            boolean isFiltered = variant.isFiltered(sampleIdentifier);
            if (isFiltered && storage.parameters.maskFiltered()) {
                // Add a placeholder for filtered variants if masking is enabled.
                variantsMap.put(variant.position, Constants.ANY_NUCLEOTIDE);
            } else if (!isFiltered) {
                // Add the alternative allele for unfiltered variants.
                variantsMap.put(variant.position, variant.alternative);
            }
        }

        // Return the integrated sequence with the mapped variants.
        return Bio.integrateVariants(context, variantsMap, !aligned);
    }

    /**
     * Validates the provided contig to ensure it is present in the storage and meets the requirements for conserved sequence generation.
     *
     * @param storage   The storage object containing contigs and associated data.
     * @param contig    The contig to validate.
     * @param conserved A flag indicating whether the sequence generation is conserved.
     * @throws IllegalArgumentException If the contig is not found in the storage or if a conserved sequence is requested without a
     *                                  reference.
     */
    void validateContig(Storage storage, Contig contig, boolean conserved) {
        // Check if the contig exists in the storage.
        if (!storage.hasContig(contig._id)) {
            throw new IllegalArgumentException("Contig %s does not exist in storage.".formatted(contig._id));
        }
        // Ensure a reference sequence is available for conserved sequence generation.
        if (!contig.hasSequence() && conserved) {
            throw new IllegalArgumentException("Cannot initialize conserved sequence generation without a reference.");
        }
    }

    /**
     * Validates the provided feature to ensure it is present in the storage.
     *
     * @param storage The storage object containing features and associated data.
     * @param feature The feature to validate.
     * @throws IllegalArgumentException If the feature is not found in the storage.
     */
    void validateFeature(Storage storage, Feature feature) {
        // Check if the feature exists in the storage.
        if (!storage.hasFeature(feature._id)) {
            throw new IllegalArgumentException("Feature %s does not exist in storage.".formatted(feature._id));
        }
    }

    /**
     * Validates the provided sample identifier to ensure it exists in the storage and is part of the allowed sample identifiers.
     * <p>
     * This method checks if the sample identifier exists in the storage. If multiple sample identifiers are provided, it also verifies that
     * the given sample identifier is included in the list of allowed identifiers.
     *
     * @param sampleIdentifier The sample identifier to validate.
     * @throws IllegalArgumentException If the sample does not exist in the storage or is not in the allowed sample identifiers.
     */
    void validateSample(String sampleIdentifier) {
        // Check if the sample exists in the storage.
        if (!storage.hasSample(sampleIdentifier))
            throw new IllegalArgumentException("Sample %s does not exist in storage.".formatted(sampleIdentifier));

        // If multiple sample identifiers are provided, ensure the given identifier is in the list.
        if (!(sampleIdentifiers.isEmpty() || sampleIdentifiers.contains(sampleIdentifier)))
            throw new IllegalArgumentException("Sample " + sampleIdentifier + " not in sample identifiers.");
    }

    /**
     * Determines if a given variant is filtered. A variant is considered filtered if (1.) there are sample identifiers provided, and the
     * variant is filtered for those samples or (2.) the variant is filtered globally.
     *
     * @param variant The variant to check.
     * @return true if the variant is filtered, false otherwise.
     */
    private boolean filtered(Variant variant) {
        return (!sampleIdentifiers.isEmpty() && variant.isFiltered(sampleIdentifiers)) || variant.isFiltered();
    }

    /**
     * Checks if a given variant is unrelated to the current sample identifiers. A variant is considered unrelated if there are sample
     * identifiers provided and the variant is not associated with those sample identifiers.
     *
     * @param variant The variant to check.
     * @return true if the variant is unrelated to the sample identifiers, false otherwise.
     */
    private boolean unrelated(Variant variant) {
        return !sampleIdentifiers.isEmpty() && !variant.ofSamples(sampleIdentifiers);
    }

    /**
     * Generates the nucleotide context for the sequence generator.
     * <p>
     * This method initializes the `context` map by processing variants within the specified interval or the entire contig if no interval is
     * defined. It handles different types of variants (SNV, DELETION, INSERTION) and integrates them into the context map.
     * <p>
     * Variants that are unrelated or filtered (and not masked) are skipped. The context map is implemented using a BTreeMap for efficient
     * storage and retrieval.
     *
     * @throws IOException           See {@link AminoacidSequenceGenerator#generateContext()}.
     * @throws MusialException       See {@link AminoacidSequenceGenerator#generateContext()}.
     * @throws IllegalStateException If an unexpected variant type is encountered.
     */
    void generateContext() throws IOException, MusialException {
        List<Variant> variants;

        // Retrieve the relevant variants based on the interval or the entire contig.
        variants = (interval != null)
                ? contig.getVariantsWithin(interval.a, interval.b)
                : contig.getAllVariants();

        // Initialize the context map using a BTreeMap.
        context = BTreeMap.create();

        // Process each variant in the list.
        for (Variant variant : variants) {
            // Skip unrelated or filtered variants if masking is not enabled.
            if (unrelated(variant) || (filtered(variant) && !storage.parameters.maskFiltered())) continue;

            // Handle the variant based on its type.
            switch (variant.type) {
                case SNV ->
                    // Merge single nucleotide variants (SNVs) into the context map.
                        context.merge(variant.position, new Bio.ReferenceContext(variant.reference.charAt(0), 0),
                                (x, y) -> new Bio.ReferenceContext(variant.reference.charAt(0), Math.max(x.extension(), y.extension())));
                case DELETION -> {
                    // Handle deletions by iterating through the reference sequence.
                    for (int i = 0; i < variant.reference.length(); i++) {
                        char c = variant.reference.charAt(i);
                        context.merge(variant.position + i, new Bio.ReferenceContext(c, 0),
                                (x, y) -> new Bio.ReferenceContext(c, Math.max(x.extension(), y.extension())));
                    }
                }
                case INSERTION -> {
                    // Handle insertions by calculating the insertion length.
                    int insertionLength = filtered(variant) ? 0 : variant.alternative.length() - 1;
                    context.merge(variant.position, new Bio.ReferenceContext(variant.reference.charAt(0), insertionLength),
                            (x, y) -> new Bio.ReferenceContext(variant.reference.charAt(0), Math.max(x.extension(), y.extension())));
                }
                default ->
                    // Throw an exception for unexpected variant types.
                        throw new IllegalStateException("Unexpected variant type " + variant.type);
            }
        }
    }

    /**
     * Sets the nucleotide context for the sequence generator.
     * <p>
     * This method initializes the `context` map including conserved positions. It processes the reference sequence and calculates the
     * maximum insertion length for each position. Variants of type INSERTION are checked at each position, and their lengths are used to
     * update the context. The context map is implemented using a BTreeMap for efficient storage and retrieval.
     *
     * @throws IOException     See {@link AminoacidSequenceGenerator#generateConservedContext()}.
     * @throws MusialException See {@link AminoacidSequenceGenerator#generateConservedContext()}.
     */
    void generateConservedContext() throws IOException, MusialException {
        // Define the start and end positions of the sequence.
        int start, end;
        if (interval == null) {
            start = 1;
            end = contig.getSequenceLength();
        } else {
            start = interval.a;
            end = interval.b;
        }

        // Get the reference sequence as a character array.
        char[] referenceCharacters = (interval != null)
                ? contig.getSequence(interval.a, interval.b).toCharArray()
                : contig.getSequence().toCharArray();

        // Initialize the context map with the size of the reference sequence.
        context = BTreeMap.create();

        // Iterate through the sequence positions.
        for (int i = start, x = 0; i <= end; i++, x++) {
            int maxInsertionLength = 0;

            // Check for insertion variants at the current position.
            for (Variant variant : contig.getVariantsAt(i)) {
                // If the variant is an insertion and not filtered, update the maximum insertion length.
                if (variant.type.equals(Variant.Type.INSERTION) && !filtered(variant)) {
                    maxInsertionLength = Math.max(maxInsertionLength, variant.alternative.length() - 1);
                }
            }

            // Add the base and maximum insertion length to the context.
            context.put(i, new Bio.ReferenceContext(referenceCharacters[x], maxInsertionLength));
        }
    }

}
