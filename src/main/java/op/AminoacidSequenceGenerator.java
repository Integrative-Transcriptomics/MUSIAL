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

public class AminoacidSequenceGenerator extends NucleotideSequenceGenerator {

    /**
     * Reference amino acid sequence of the feature.
     */
    private char[] reference;

    /**
     * Allele identifiers related to the provided sample identifiers.
     */
    private String[] alleleIdentifiers;

    /**
     * Constructs a new instance of the AminoacidSequenceGenerator class.
     * <p>
     * This constructor initializes the generator for amino acid sequence generation by calling the superclass constructor. It passes the
     * provided storage, contig, feature, conserved flag, aligned flag, and sample identifiers to the superclass for initialization.
     * <p>
     * In contrast to the {@link NucleotideSequenceGenerator}, this class always defines a {@link Feature} and requires an associated parent
     * contig with a reference sequence. The feature must be a coding feature to ensure that amino acid sequences can be generated. In
     * addition, the {@link #interval} is always defined as the full length of the amino acid sequence derived from the feature's nucleotide
     * sequence.
     *
     * @param storage           The storage object containing genomic data.
     * @param contig            The contig associated with the sequence generation.
     * @param feature           The feature associated with the sequence generation.
     * @param conserved         A flag indicating whether the sequence generation includes conserved sites.
     * @param aligned           A flag indicating whether the generated sequences are aligned.
     * @param sampleIdentifiers Optional sample identifiers to restrict the sequence generation.
     * @throws IOException     If an error occurs during sequence retrieval.
     * @throws MusialException If an error occurs during initialization.
     */
    public AminoacidSequenceGenerator(Storage storage, Contig contig, Feature feature, boolean conserved, boolean aligned,
                                      String... sampleIdentifiers) throws IOException, MusialException {
        super(storage, contig, feature, conserved, aligned, sampleIdentifiers);
    }

    /**
     * Retrieves the amino acid sequence for a given sample identifier.
     * <p>
     * This method generates the sequence for the specified sample by integrating variants associated with the sample's allele and
     * proteoform. If the sequence for the proteoform is already cached, it is returned directly. Otherwise, the sequence is generated,
     * cached, and returned.
     *
     * @param sampleIdentifier The identifier of the sample for which the sequence is generated.
     * @return A {@link String} representing the nucleotide sequence for the given sample.
     * @throws MusialException If the sample identifier is invalid or an error occurs during sequence generation.
     */
    @Override
    public String getSequence(String sampleIdentifier) throws MusialException {
        // Validate the provided sample identifier.
        super.validateSample(sampleIdentifier);

        // Retrieve the allele identifier associated with the sample and feature.
        String alleleIdentifier = storage.getSample(sampleIdentifier).getRelatedAllele(feature._id);
        String proteoformIdentifier;

        // Determine the proteoform identifier based on the allele identifier.
        if (alleleIdentifier.equals(Constants.REFERENCE)) {
            proteoformIdentifier = Constants.SYNONYMOUS;
        } else {
            proteoformIdentifier = feature.getAllele(alleleIdentifier).getRelatedProteoform();
        }

        // Return the cached sequence if it is already available.
        if (cache.containsKey(proteoformIdentifier)) {
            return cache.get(proteoformIdentifier);
        }

        // Generate the sequence based on the proteoform identifier.
        String sequence = proteoformIdentifier.equals(Constants.SYNONYMOUS)
                ? Bio.integrateVariants(context, Collections.emptyMap(), !aligned)
                : Bio.integrateVariants(context, feature.getProteoform(proteoformIdentifier).getVariants(), !aligned);

        // Cache the generated sequence and return it.
        cache.put(proteoformIdentifier, sequence);
        return sequence;
    }

    /**
     * Validates the provided contig to ensure it is suitable for amino acid sequence generation.
     * <p>
     * This method first calls the superclass implementation to perform general contig validation. It then performs additional validation
     * specific to amino acid sequence generation, ensuring that the contig has an associated reference sequence.
     *
     * @param storage   The storage object containing genomic data.
     * @param contig    The contig to validate.
     * @param conserved A flag indicating whether the sequence generation includes conserved sites.
     * @throws IllegalArgumentException If the contig fails the superclass validation or does not have a reference sequence.
     */
    @Override
    protected void validateContig(Storage storage, Contig contig, boolean conserved) {
        // Validate contig as in super class.
        super.validateContig(storage, contig, conserved);

        // Additional validation for coding features.
        if (!contig.hasSequence()) {
            throw new IllegalArgumentException("Cannot initialize amino acid sequence generation without reference.");
        }
    }

    /**
     * Validates the provided feature to ensure it is suitable for amino acid sequence generation.
     * <p>
     * This method first calls the superclass implementation to perform general feature validation. It then performs additional validation
     * specific to amino acid sequence generation, ensuring that the feature is a coding feature.
     *
     * @param storage The storage object containing genomic data.
     * @param feature The feature to validate.
     * @throws IllegalArgumentException If the feature is not coding or fails the superclass validation.
     */
    @Override
    protected void validateFeature(Storage storage, Feature feature) {
        // Validate feature as in super class.
        super.validateFeature(storage, feature);

        // Additional validation for coding features.
        if (!feature.isCoding()) {
            throw new IllegalArgumentException("Cannot initialize amino acid sequence generation with non coding feature %s.".formatted(feature._id));
        }
    }

    /**
     * Generates the amino acid context for the sequence generator.
     * <p>
     * This method initializes the `context` map by processing variants related to the {@link #alleleIdentifiers} inferred from the
     * {@link #sampleIdentifiers} and the associated {@link #feature}. It collects variants from the proteoforms of the feature, filters
     * them based on their relation to the allele identifiers, and processes each variant to update the context map. The context map is
     * implemented using a BTreeMap for efficient storage and retrieval.
     *
     * @throws IOException     If an error occurs during initialization or sequence retrieval.
     * @throws MusialException If an error occurs during initialization or variant processing.
     */
    @Override
    protected void generateContext() throws IOException, MusialException {
        // Initialize the reference sequence and allele identifiers.
        initialize();

        // Collect variants related to the allele identifiers.
        Set<Variant.Stub> variants = this.feature.getProteoforms().stream()
                .filter(proteoform -> Arrays.stream(this.alleleIdentifiers).anyMatch(proteoform::hasRelation))
                .flatMap(proteoform -> proteoform.getVariants().entrySet().stream())
                .map(entry -> new Variant.Stub(entry.getKey(), entry.getValue()))
                .collect(Collectors.toSet());

        // Initialize the context map using a BTreeMap.
        this.context = BTreeMap.create();

        // Process each variant in the list.
        for (Variant.Stub variant : variants) {
            char r = reference[variant.position() - 1];
            if (Bio.isSubstitution(variant.alternative())) {
                // Handle substitution variants.
                this.context.merge(variant.position(), new Bio.ReferenceContext(r, 0),
                        (x, y) -> new Bio.ReferenceContext(r, Math.max(x.extension(), y.extension())));
            } else if (Bio.isDeletion(variant.alternative())) {
                // Handle deletion variants by iterating through the alternative sequence.
                for (int i = 0; i < variant.alternative().length(); i++) {
                    char refChar = reference[variant.position() - 1 + i];
                    this.context.merge(variant.position() + i, new Bio.ReferenceContext(refChar, 0),
                            (x, y) -> new Bio.ReferenceContext(refChar, Math.max(x.extension(), y.extension())));
                }
            } else if (Bio.isInsertion(variant.alternative())) {
                // Handle insertion variants by calculating the insertion length.
                int insertionLength = variant.alternative().length() - 1;
                this.context.merge(variant.position(), new Bio.ReferenceContext(r, insertionLength),
                        (x, y) -> new Bio.ReferenceContext(r, Math.max(x.extension(), y.extension())));
            } else {
                // Throw an exception for unsupported variant types.
                throw new IllegalStateException("Unable to determine type for variant %s:p.%d?>%s.".formatted(
                        feature._id, variant.position(), variant.alternative()));
            }
        }
    }

    /**
     * Generates the conserved nucleotide context for the sequence generator.
     * <p>
     * This method initializes the `context` map by processing the reference sequence and calculating the maximal insertion lengths for
     * positions related to the allele identifiers. The context map is populated with reference bases and their corresponding insertion
     * lengths.
     *
     * @throws IOException     If an error occurs during initialization or sequence retrieval.
     * @throws MusialException If an error occurs during initialization or variant processing.
     */
    @Override
    protected void generateConservedContext() throws IOException, MusialException {
        // Initialize the reference sequence and allele identifiers.
        initialize();

        // Collect maximal insertion lengths for positions related to allele identifiers.
        Map<Integer, Integer> maximalInsertionLengths = feature.getProteoforms().stream()
                .filter(proteoform -> Arrays.stream(alleleIdentifiers).anyMatch(proteoform::hasRelation))
                .flatMap(proteoform -> proteoform.getVariants().entrySet().stream())
                .filter(entry -> Bio.isInsertion(entry.getValue()))
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        entry -> entry.getValue().length() - 1,
                        Math::max
                ));

        // Initialize the context map and populate it with reference bases and insertion lengths.
        this.context = BTreeMap.create();
        for (int i = interval.a; i <= interval.b; i++) {
            int insertionLength = maximalInsertionLengths.getOrDefault(i, 0);
            this.context.put(i, new Bio.ReferenceContext(reference[i - 1], insertionLength));
        }
    }

    /**
     * Initializes the state of the AminoacidSequenceGenerator.
     * <p>
     * This method prepares the reference sequence, allele identifiers, interval, and cache required for amino acid sequence generation. It
     * translates the nucleotide sequence of the contig into an amino acid sequence, identifies relevant alleles, and sets up the interval
     * and cache for further processing.
     *
     * @throws IOException     If an error occurs during sequence retrieval.
     * @throws MusialException If an error occurs during initialization.
     */
    private void initialize() throws IOException, MusialException {
        // Translate the nucleotide sequence of the contig into an amino acid sequence.
        this.reference = Bio.translateSequence(this.contig.getSequence(this.feature.start, this.feature.end),
                this.feature.isReverse()).toCharArray();

        // Define the interval for sequence generation.
        this.interval = new Tuple<>(1, this.reference.length);

        // Identify allele identifiers related to the provided sample identifiers.
        this.alleleIdentifiers = Arrays.stream(this.sampleIdentifiers)
                .map(sampleIdentifier -> this.storage.getSample(sampleIdentifier).getRelatedAllele(this.feature._id))
                .filter(alleleIdentifiers -> !Objects.equals(alleleIdentifiers, Constants.REFERENCE))
                .toArray(String[]::new);

        // Initialize the cache for storing proteoform-related data.
        this.cache = new HashMap<>(feature.getProteoformCount() + 1);
    }

}
