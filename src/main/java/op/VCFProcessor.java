package op;

import htsjdk.samtools.util.Tuple;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import main.Musial;
import model.*;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.Strings;
import org.apache.commons.lang3.tuple.Triple;
import org.ehcache.Cache;
import org.ehcache.CacheManager;
import org.ehcache.config.builders.CacheConfigurationBuilder;
import org.ehcache.config.builders.CacheManagerBuilder;
import org.ehcache.config.builders.ExpiryPolicyBuilder;
import org.ehcache.config.builders.ResourcePoolsBuilder;
import org.ehcache.config.units.EntryUnit;
import org.ehcache.config.units.MemoryUnit;
import uk.co.omegaprime.btreemap.BTreeMap;
import util.Bio;
import util.Constants;
import util.IO;
import util.Logging;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Consumer;

/**
 * The {@code VCFProcessor} class is responsible for processing Variant Call Format (VCF) files.
 * <p>
 * This class handles the analysis of VCF files, including the extraction of variant data, imputation of contigs, and integration of the
 * processed data into the storage system. It provides methods to analyze VCF files, process variant contexts, and track statistics such as
 * the number of processed, ignored, and filtered variant calls.
 *
 * @noinspection DuplicatedCode
 */
public class VCFProcessor implements Closeable {

    /**
     * The storage object for managing processed genomic data.
     */
    private final Storage storage;

    /**
     * A list of file paths to the VCF files to be processed.
     * <p>
     * This field stores the paths to the Variant Call Format (VCF) files that will be read and processed.
     */
    private final List<Path> paths;

    /**
     * A flag indicating whether to impute contigs from the VCF files.
     * <p>
     * If set to {@code true}, the program will infer and add contigs to the storage based on the contigs present in the VCF files. This is
     * useful when the contigs are not explicitly defined in the input data.
     */
    private final boolean imputeContigs;

    /**
     * Tracks the total number of processed variant calls.
     * <p>
     * This field is used to count the number of variant calls that have been processed during the execution of the program. It is
     * incremented each time a variant call is processed, regardless of its outcome.
     */
    private long processedCallsCount = 0;

    /**
     * Tracks the total number of ignored variant calls.
     * <p>
     * This field is used to count the number of variant calls that were ignored during processing. A call may be ignored for various
     * reasons, such as missing data or being classified as a reference call.
     */
    private long ignoredCallsCount = 0;

    /**
     * Tracks the total number of filtered variant calls.
     * <p>
     * This field is used to count the number of variant calls that were filtered out during processing. Filtering may occur due to low
     * coverage, low frequency, or other criteria defined in the program.
     */
    private long filteredCallsCount = 0;

    /**
     * A map that tracks upstream deletions for each sample.
     * <p>
     * The key is the sample identifier, and the value is a {@link UpstreamDeletion} object that represents the details of the upstream
     * deletion, including its genomic range and filter status. This map is used to handle cases where a deletion affects downstream
     * positions in the genome.
     */
    private final Map<String, UpstreamDeletion> upstreamDeletions = new HashMap<>(1_000);

    /**
     * Represents an upstream deletion affecting a sample.
     * <p>
     * This record encapsulates the details of an upstream deletion, including:
     * <ul>
     *   <li>The contig identifier where the deletion occurs.</li>
     *   <li>The start position of the deletion.</li>
     *   <li>The end position of the deletion.</li>
     *   <li>A flag indicating whether the deletion is filtered.</li>
     * </ul>
     *
     * @param contigIdentifier The unique identifier of the contig where the deletion occurs.
     * @param start            The start position of the deletion.
     * @param end              The end position of the deletion.
     * @param isFiltered       A boolean flag indicating whether the deletion is filtered.
     */
    public record UpstreamDeletion(String contigIdentifier, int start, int end, boolean isFiltered) {
    }

    /**
     * A cache for storing and managing variant calls.
     * <p>
     * The {@link VariantCallCache} is responsible for indexing, storing, and retrieving variant calls based on sample identifiers, contig
     * identifiers, and genomic positions. It provides efficient access to variant data during processing and ensures that calls are merged
     * and updated as needed.
     */
    private final VariantCallCache vcCache;

    /**
     * A cache for storing and managing variant calls.
     * <p>
     * This class provides methods to index, store, and retrieve variant calls based on sample identifiers, contig identifiers, and genomic
     * positions. It uses a hierarchical map structure for efficient indexing and an Ehcache instance for caching the variant call objects.
     */
    private static class VariantCallCache {

        /**
         * A unique alias for the cache, generated as a random alphanumeric string.
         */
        final static String ALIAS = IO.randomAlphanumeric(4).toUpperCase();

        /**
         * Counter for indexing sample identifiers.
         */
        private int sampleIndex = 0;

        /**
         * Counter for indexing contig identifiers.
         */
        private int contigIndex = 0;

        /**
         * Counter for generating unique cache keys.
         */
        private long cacheKey = 0;

        /**
         * A map for storing sample identifiers and their corresponding indices.
         */
        final Map<String, Integer> sampleMap = new HashMap<>(1_000);

        /**
         * A map for storing contig identifiers and their corresponding indices.
         */
        final Map<String, Integer> contigMap = new HashMap<>(10);

        /**
         * A hierarchical map structure for indexing variant calls by sample, contig, and position.
         */
        final Map<Integer, Map<Integer, Map<Integer, Long>>> callMap = new HashMap<>(1_000);

        /**
         * The cache manager for managing the Ehcache instance.
         */
        final CacheManager manager = CacheManagerBuilder
                .newCacheManagerBuilder()
                .with(CacheManagerBuilder.persistence(Musial.tempDir))
                .withCache(ALIAS,
                        CacheConfigurationBuilder
                                .newCacheConfigurationBuilder(Long.class, VariantCall.class,
                                        ResourcePoolsBuilder.newResourcePoolsBuilder()
                                                .heap(100_000_000, EntryUnit.ENTRIES)
                                                .offheap(2, MemoryUnit.GB)
                                                .disk(20, MemoryUnit.GB, false))
                                .withValueSerializer(new VariantCall.VariantCallSerializer())
                                .withExpiry(ExpiryPolicyBuilder.noExpiration()))
                .build(true);

        /**
         * The Ehcache instance for storing variant calls.
         */
        private final Cache<Long, VariantCall> cache = manager.getCache(ALIAS, Long.class, VariantCall.class);

        /**
         * Private constructor to initialize the VariantCallCache.
         */
        private VariantCallCache() {
        }

        /**
         * Indexes a sample identifier by assigning it a unique index.
         *
         * @param identifier The sample identifier to index.
         */
        void indexSample(String identifier) {
            sampleMap.putIfAbsent(identifier, sampleIndex++);
        }

        /**
         * Indexes a contig identifier by assigning it a unique index.
         *
         * @param identifier The contig identifier to index.
         */
        void indexContig(String identifier) {
            contigMap.putIfAbsent(identifier, contigIndex++);
        }

        /**
         * Stores a variant call in the cache.
         *
         * @param sampleIdentifier The sample identifier associated with the variant call.
         * @param contigIdentifier The contig identifier associated with the variant call.
         * @param position         The genomic position of the variant call.
         * @param call             The {@link VariantCall} object to store.
         */
        void put(String sampleIdentifier, String contigIdentifier, int position, VariantCall call) {
            Long key = callMap.computeIfAbsent(sampleMap.get(sampleIdentifier), k -> new HashMap<>(10))
                    .computeIfAbsent(contigMap.get(contigIdentifier), k -> BTreeMap.create())
                    .putIfAbsent(position, cacheKey);
            if (Objects.nonNull(key))  // Replace existing call.
                cache.put(key, call);
            else  // New call.
                cache.put(cacheKey++, call);
        }

        /**
         * Retrieves a variant call from the cache based on sample, contig, and position.
         *
         * @param sampleIdentifier The sample identifier associated with the variant call.
         * @param contigIdentifier The contig identifier associated with the variant call.
         * @param position         The genomic position of the variant call.
         * @return The {@link VariantCall} object if found, or {@code null} if not found.
         */
        VariantCall get(String sampleIdentifier, String contigIdentifier, int position) {
            Map<Integer, Map<Integer, Long>> sampleMapEntry = callMap.get(sampleMap.get(sampleIdentifier));
            if (sampleMapEntry == null) return null;

            Map<Integer, Long> contigMapEntry = sampleMapEntry.get(contigMap.get(contigIdentifier));
            if (contigMapEntry == null) return null;

            Long positionKey = contigMapEntry.get(position);
            return positionKey != null ? cache.get(positionKey) : null;
        }

        /**
         * Retrieves all variant calls for a specific sample and contig.
         *
         * @param sampleIdentifier The sample identifier associated with the variant calls.
         * @param contigIdentifier The contig identifier associated with the variant calls.
         * @return A {@link List} of {@link Tuple} objects, where each tuple contains the position and the corresponding variant call.
         */
        List<Tuple<Integer, VariantCall>> get(String sampleIdentifier, String contigIdentifier) {
            if (!sampleMap.containsKey(sampleIdentifier) || !contigMap.containsKey(contigIdentifier))
                return Collections.emptyList();
            Set<Map.Entry<Integer, Long>> entries = callMap.getOrDefault(sampleMap.get(sampleIdentifier), Collections.emptyMap())
                    .getOrDefault(contigMap.get(contigIdentifier), Collections.emptyMap())
                    .entrySet();
            List<Tuple<Integer, VariantCall>> calls = new ArrayList<>(entries.size());
            for (var entry : entries) {
                calls.add(new Tuple<>(entry.getKey(), this.cache.get(entry.getValue())));
            }
            return calls;
        }
    }

    /**
     * Constructs a new {@link VCFProcessor} instance for processing VCF files.
     * <p>
     * This constructor initializes the processor with the specified list of VCF file paths, a storage object for managing genomic data, and
     * a flag indicating whether to impute contigs from the VCF files. The {@code imputeContigs} flag determines if contigs should be
     * inferred and added to the storage based on the VCF data.
     *
     * @param paths         A {@link List} of {@link Path} objects representing the file paths to the VCF files to be processed.
     * @param storage       The {@link Storage} object used to store and manage processed genomic data.
     * @param imputeContigs A {@code boolean} flag indicating whether to infer and add contigs from the VCF files to the storage.
     */
    public VCFProcessor(List<Path> paths, Storage storage, boolean imputeContigs) {
        this.paths = paths;
        this.storage = storage;
        this.vcCache = new VariantCallCache();
        this.imputeContigs = imputeContigs;
    }

    public void close() {
        this.vcCache.manager.close();
    }

    /**
     * Analyzes the provided VCF files and processes their variant data.
     * <p>
     * This method iterates through the list of VCF file paths, creates temporary indexed VCF files, and processes their variant contexts.
     * If the `imputeContigs` flag is set, it infers contigs from the VCF files and adds them to the storage. The method processes variant
     * contexts either for specific features in the storage or for all variants if no features are defined.
     *
     * @throws IOException If an I/O error occurs during file operations or VCF processing.
     */
    public void processFiles() throws IOException {
        // Iterate over each VCF file path in the list of paths.
        for (Path path : paths) {
            // Create a temporary index file for the VCF file if it does not already exist.
            File indexFile = new File(path + ".idx");
            if (!indexFile.exists()) {
                IndexFactory.createLinearIndex(path.toFile(), new VCFCodec()).write(indexFile);
            }

            // Open the VCF file for reading using a VCFFileReader.
            try (VCFFileReader vcfFileReader = new VCFFileReader(path)) {
                // Extract the header of the VCF file.
                VCFHeader vcfHeader = vcfFileReader.getHeader();

                // Check if the VCF file contains genotyping data.
                if (!vcfHeader.hasGenotypingData()) {
                    // Log a warning if no genotyping data is found.
                    Logging.logWarning("VCF file %s does not contain genotyping data.".formatted(path));
                } else {
                    // Index contigs from the VCF file and add them to the storage if the imputeContigs flag is set.
                    for (VCFContigHeaderLine contigLine : vcfHeader.getContigLines()) {
                        String contigId = contigLine.getID(); // Extract the contig ID.
                        this.vcCache.indexContig(contigId); // Index the contig in the store.
                        if (imputeContigs) {
                            storage.addContig(contigId, Constants.EMPTY); // Add the contig to storage.
                        }
                    }

                    // Index the sample identifiers from the VCF header.
                    vcfHeader.getGenotypeSamples().forEach(this.vcCache::indexSample);

                    // Check if specific features are defined in the storage.
                    if (!storage.getFeatures().isEmpty()) {
                        // Process variants for each feature in the storage.
                        for (Feature feature : storage.getFeatures()) {
                            processVariantContexts(vcfFileReader.query(feature.contig, feature.start, feature.end), path);
                        }
                    } else {
                        // Process all variants if no specific features are defined.
                        // Note: Optimization may be required for large files.
                        processVariantContexts(vcfFileReader.iterator(), path);
                    }
                }
            }
        }
    }

    /**
     * Retrieves the total number of processed variant calls.
     *
     * @return The total number of processed variant calls as a {@code long}.
     */
    public long getProcessedCallsCount() {
        return processedCallsCount;
    }

    /**
     * Retrieves the total number of ignored variant calls.
     * <p>
     * This method returns the count of variant calls that were ignored during processing. A variant call may be ignored for reasons such as
     * missing data, being classified as a reference call, or lacking sufficient information for analysis.
     *
     * @return The total number of ignored variant calls as a {@code long}.
     */
    public long getIgnoredCallsCount() {
        return ignoredCallsCount;
    }

    /**
     * Retrieves the total number of filtered variant calls.
     * <p>
     * This method returns the count of variant calls that were filtered out during processing. Filtering may occur due to criteria such as
     * low coverage, low frequency, or other conditions defined in the program.
     *
     * @return The total number of filtered variant calls as a {@code long}.
     */
    public long getFilteredCallsCount() {
        return filteredCallsCount;
    }

    /**
     * Retrieves an unmodifiable collection of sample identifiers.
     * <p>
     * This method returns a collection of all sample identifiers currently indexed in the variant call cache. The returned collection is
     * unmodifiable, ensuring that the underlying data cannot be altered.
     *
     * @return An unmodifiable {@link Collection} of {@link String} objects representing the sample identifiers.
     */
    public Collection<String> getSamples() {
        return Collections.unmodifiableCollection(this.vcCache.sampleMap.keySet());
    }

    /**
     * Updates the variants in the storage by processing variant calls from the cache.
     * <p>
     * This method iterates through all samples and contigs in the cache, processes their variant calls, and adds the resolved variants to
     * the storage. It handles complex cases such as deletions, insertions, and mixed InDels, ensuring that the variants are stored in a
     * canonical format.
     */
    public void updateVariants() {
        // Iterate through all samples in the cache's map.
        for (String sampleIdentifier : this.vcCache.sampleMap.keySet()) {
            // Add the sample to the storage.
            storage.addSample(sampleIdentifier);
            assert storage.hasSample(sampleIdentifier);
            Sample sample = storage.getSample(sampleIdentifier);

            // Iterate through all contigs in the storage.
            for (String contigIdentifier : this.vcCache.contigMap.keySet()) {
                assert storage.hasContig(contigIdentifier);
                Contig contig = storage.getContig(contigIdentifier);

                // Access all variant calls for the current sample and contig from the cache.
                List<Tuple<Integer, VariantCall>> entries = vcCache.get(sampleIdentifier, contigIdentifier);
                // Continue if no calls are present for the current sample and contig.
                if (entries.isEmpty()) continue;

                // Initialize builders and variables for processing variants, esp. handling deletions and mixed InDels.
                StringBuilder referenceBuilder = new StringBuilder();
                StringBuilder alternativeBuilder = new StringBuilder();
                Set<VariantCall> calls = new HashSet<>(); // Set of variant calls that constitute the current variant.
                int variantStartPosition = 0;
                int deletionExtension = 0;

                // Helper function to resolve variants and add them to the storage.
                Consumer<Integer> resolve = (position) -> {
                    String referenceContent = referenceBuilder.toString();
                    String alternativeContent = alternativeBuilder.toString();
                    if (!Bio.isPaddedCanonical(referenceContent, alternativeContent)) {
                        Tuple<String, String> alignment =
                                Bio.globalNucleotideSequenceAlignment(
                                        Bio.stripGaps(referenceBuilder.toString()),
                                        Bio.stripGaps(alternativeBuilder.toString()),
                                        5,
                                        2,
                                        true,
                                        false,
                                        0
                                );
                        ArrayList<Triple<Integer, String, String>> resolvedVariants =
                                Bio.getCanonicalVariants(alignment.a, alignment.b);
                        for (Triple<Integer, String, String> variant : resolvedVariants) {
                            storage.addVariant(contig, sample._id, position + variant.getLeft(), variant.getMiddle(), variant.getRight(),
                                    calls);

                        }
                    } else {
                        storage.addVariant(contig, sample._id, position, referenceContent, alternativeContent, calls);
                    }
                    // Reset outer state.
                    referenceBuilder.setLength(0);
                    alternativeBuilder.setLength(0);
                    calls.clear();
                };

                // Process variant calls for the current sample and contig.
                for (var entry : entries) {
                    int position = entry.a; // Extract the position of the variant.
                    VariantCall variantCall = entry.b; // Extract the variant call.

                    // Resolves the reference and alternative alleles for a variant call, handling special cases.
                    String reference = variantCall.isFiltered() && variantCall.flag() == VariantCall.Flag.MISSING_UPSTREAM_DELETION
                            ? variantCall.getReference(1) // Retrieve the reference content for the second alternative if flagged.
                            : variantCall.getReference(); // Retrieve the called reference content.
                    String alternative = variantCall.isFiltered() && variantCall.flag() == VariantCall.Flag.MISSING_UPSTREAM_DELETION
                            ? variantCall.getAlternative(1) // Retrieve the alternative content for the second position if flagged.
                            : variantCall.getAlternative(); // Retrieve the default alternative content.

                    // Replace the alternative content with an N, if it represents a filtered reference call.
                    if (alternative.equals(Constants.DOT)) {
                        alternative = Constants.ANY_NUCLEOTIDE;
                    }

                    // Resolve previous variant, if the current position is outside the stored upstream deletion range.
                    if (position > deletionExtension && referenceBuilder.length() > 0 && alternativeBuilder.length() > 0) {
                        resolve.accept(variantStartPosition);
                        deletionExtension = 0;
                    }

                    // Start processing a new variant if no ongoing deletion exists.
                    if (deletionExtension == 0 && referenceBuilder.length() == 0 && alternativeBuilder.length() == 0) {
                        if (Bio.isDeletion(reference, alternative, true)) {
                            referenceBuilder.append(reference);
                            alternativeBuilder.append(alternative);
                            calls.add(variantCall);
                            variantStartPosition = position;
                            deletionExtension = position + alternative.length() - 1;
                        } else {
                            referenceBuilder.append(reference);
                            alternativeBuilder.append(alternative);
                            calls.add(variantCall);
                            resolve.accept(position);
                        }
                        continue;
                    }

                    // Extend ongoing deletions or handle insertions within the deletion range.
                    if (position <= deletionExtension) {
                        if (Bio.isSubstitution(reference, alternative)) {
                            continue; // Skip substitutions within an upstream deletion.
                        }
                        if (Bio.isDeletion(reference, alternative, true)) {
                            int updatedDeletionExtension = position + alternative.length() - 1;
                            if (updatedDeletionExtension > deletionExtension) {
                                referenceBuilder.append(StringUtils.right(reference, updatedDeletionExtension - deletionExtension));
                                alternativeBuilder.append(StringUtils.right(alternative,
                                        updatedDeletionExtension - deletionExtension));
                                calls.add(variantCall);
                                deletionExtension = updatedDeletionExtension;
                            }
                            continue;
                        }
                        if (Bio.isInsertion(reference, alternative, true)) {
                            int offset = position - variantStartPosition;
                            alternativeBuilder.replace(offset, offset + 1,
                                    alternativeBuilder.charAt(offset) + alternative.substring(1));
                            referenceBuilder.replace(offset, offset + 1, referenceBuilder.charAt(offset) + reference.substring(1));
                            calls.add(variantCall);
                            continue;
                        }
                    }

                    // Log a warning if the variant cannot be handled.
                    Logging.logWarning("Failed to handle variant %s at site %s %d for sample %s."
                            .formatted(reference + " > " + alternative, contigIdentifier, position, sample._id));
                }

                // Resolve any remaining variants in the builders after processing all variants.
                if (referenceBuilder.length() > 0 && alternativeBuilder.length() > 0) {
                    resolve.accept(variantStartPosition);
                }
            }
        }
    }

    /**
     * Loads variant calls from the storage and processes them.
     * <p>
     * This method iterates through all contigs and samples in the variant call cache (`vcCache`), retrieves the associated variants, and
     * processes their variant calls. Each variant call string is parsed into a `VariantCall` object and passed to the `processVariantCall`
     * method for further processing.
     * <p>
     * The method ensures that only valid contigs and samples present in the storage are processed. It handles the relationship between
     * contigs, samples, and variants, and updates the storage with the processed variant calls.
     *
     * @return The total number of loaded variant calls as an {@code int}.
     */
    public int loadVariantCallsFromStorage() {
        // Counter for loaded variant calls.
        int c = 0;

        // Iterate through all contig identifiers in the variant call cache.
        for (String contigIdentifier : vcCache.contigMap.keySet()) {
            // Skip contigs that are not present in the storage.
            if (!storage.hasContig(contigIdentifier)) continue;

            // Retrieve the contig object from the storage.
            Contig contig = storage.getContig(contigIdentifier);

            // Iterate through all sample identifiers in the variant call cache.
            for (String sampleIdentifier : vcCache.sampleMap.keySet()) {
                // Skip samples that are not present in the storage.
                if (!storage.hasSample(sampleIdentifier)) continue;

                // Iterate through all variants associated with the current sample in the contig.
                for (Variant variant : contig.getVariantsOfSamples(sampleIdentifier)) {
                    // Split the variant call string into individual calls.
                    for (String variantCallString : variant.getSampleRelation(sampleIdentifier).split(Constants.PIPE)) {
                        // Parse the variant call string into a VariantCall object.
                        VariantCall variantCall = VariantCall.fromString(variantCallString);

                        // Process the variant call and update the storage.
                        processVariantCall(sampleIdentifier, contigIdentifier, variant.position,
                                variantCall.alternatives(), storage.parameters, Path.of("storage"));

                        c++;
                    }
                }
            }
        }

        return c;
    }

    /**
     * Processes a set of variant contexts from a VCF file and updates the storage with variant calls.
     * <p>
     * This method iterates through the provided {@link VariantContext} objects, processes each genotype, and extracts relevant information
     * such as reference and alternative alleles, allelic depths, and coverage data. The processed data is then added to the storage as
     * variant calls.
     *
     * @param variantContextIterator An {@link Iterator} of {@link VariantContext} objects representing the variants to process.
     * @param path                   The {@link Path} to the VCF file being processed, used for logging and error reporting.
     */
    private void processVariantContexts(Iterator<VariantContext> variantContextIterator, Path path) {
        while (variantContextIterator.hasNext()) {
            VariantContext variantContext = variantContextIterator.next();
            String contigIdentifier = variantContext.getContig();

            // Skip positions excluded by the configuration.
            if (storage.parameters.isPositionMasked(contigIdentifier, variantContext.getStart())) continue;

            // Process each genotype in the current VariantContext.
            for (Genotype genotype : variantContext.getGenotypes()) {
                processedCallsCount++; // Increment the count of processed genotype records.

                if (genotype.isNoCall()) {
                    ignoredCallsCount++; // Ignore no-call genotypes.
                    continue;
                }

                // Log a warning if AD and DP attributes are missing.
                if (!(genotype.hasDP() && (genotype.hasAD() || genotype.hasAnyAttribute("COV")))) {
                    if (!genotype.isHomRef()) {
                        Logging.logWarningOnce("MISSING_AD_DP_ATTRIBUTES",
                                String.format("Variants may be ignored as AD/COV/DP4 and DP attributes are unavailable (%s:g.%d %s in %s).",
                                        contigIdentifier, variantContext.getStart(), genotype.getSampleName(), path));
                    }
                    ignoredCallsCount++; // Ignore hom-ref genotypes without coverage information.
                    continue;
                }

                // Extract the sample identifier.
                String sampleIdentifier = genotype.getSampleName();

                // Extract allelic depth (AD) information for the genotype.
                int[] ADs;
                if (genotype.hasAD()) {
                    ADs = genotype.getAD();
                } else if (genotype.hasAnyAttribute("COV")) {
                    ADs = (int[]) genotype.getAnyAttribute("COV");
                } else if (genotype.getAlleles().size() == 2
                        && variantContext.getNSamples() == 1
                        && variantContext.hasAttribute("DP4")) {
                    // Use DP4 attribute as a fallback if AD is missing.
                    List<Integer> DP4 = variantContext.getAttributeAsIntList("DP4", 0);
                    ADs = new int[]{
                            DP4.get(0) + DP4.get(1), // Reference allele coverage.
                            DP4.get(2) + DP4.get(3)  // Alternative allele coverage.
                    };
                } else {
                    ignoredCallsCount++; // Skip genotypes with insufficient information.
                    continue;
                }

                // Compute the total depth of coverage (ADSum) from allelic depths.
                int ADSum = Arrays.stream(ADs).sum();

                // Log warnings for discrepancies between AD sum and DP.
                if (genotype.hasDP()) {
                    if (ADSum > genotype.getDP()) {
                        Logging.logWarningOnce("AD_SUM_GREATER_THAN_DP",
                                String.format("Erroneous data with summed allelic depth (%d) greater than total depth (%d) (%s:g.%d %s " +
                                                "in %s).", ADSum, genotype.getDP(), contigIdentifier, variantContext.getStart(),
                                        sampleIdentifier, path));
                    } else if (ADSum < 0.5 * genotype.getDP()) {
                        Logging.logWarningOnce("AD_SUM_LOWER_THAN_DP",
                                String.format("Low-quality data with summed allelic depth (%d) much lower than total depth (%d) (%s:g.%d " +
                                                "%s in %s).", ADSum, genotype.getDP(), contigIdentifier, variantContext.getStart(),
                                        sampleIdentifier, path));
                    }
                }

                // Process reference and alternative alleles for the genotype.
                String REF;
                String ALT;
                int AD;

                // Create a list to store alternatives for the current genotype.
                List<VariantCall.CallAlternative> alternatives = new ArrayList<>(ADs.length);
                VariantCall.CallAlternative referenceCall = null;

                for (int i = 0; i < ADs.length; i++) {
                    AD = ADs[i]; // Retrieve the allelic depth for the current allele.
                    if (i == 0) {
                        // For the reference allele, set REF to the first base and ALT to a dot.
                        REF = variantContext.getReference().getBaseString().substring(0, 1);
                        ALT = Constants.DOT;
                        referenceCall = new VariantCall.CallAlternative(REF, ALT, (short) AD);
                    } else {
                        // Skip alternative alleles with zero depth.
                        if (AD == 0) continue;

                        // For alternative alleles, retrieve the full reference and alternative sequences.
                        REF = variantContext.getReference().getBaseString();
                        ALT = variantContext.getAlleles().get(i).getBaseString();

                        // Handle upstream-deletion cases and ensure canonical formatting.
                        if (ALT.equals("*")) {
                            REF = REF.substring(0, 1);
                        } else if (Bio.isCanonical(REF, ALT)) {
                            REF = Bio.padGaps(REF, ALT.length());
                            ALT = Bio.padGaps(ALT, REF.length());
                        } else {
                            String commonSuffix = StringUtils.reverse(StringUtils.getCommonPrefix(StringUtils.reverse(REF),
                                    StringUtils.reverse(ALT)));
                            REF = Strings.CS.removeEnd(REF, commonSuffix);
                            ALT = Strings.CS.removeEnd(ALT, commonSuffix);

                            if (Bio.isCanonical(REF, ALT)) {
                                REF = Bio.padGaps(REF, ALT.length());
                                ALT = Bio.padGaps(ALT, REF.length());
                            } else {
                                Logging.logDebug("Realigning non-canonical variant call %s>%s at %s:g.%d".formatted(
                                        REF, ALT, contigIdentifier, variantContext.getStart()));
                                Tuple<String, String> alignment = Bio.globalNucleotideSequenceAlignment(
                                        REF, ALT, 5, 2, true, false, 0
                                );
                                REF = alignment.a;
                                ALT = alignment.b;
                            }
                        }

                        // Add the alternative allele to the list.
                        alternatives.add(new VariantCall.CallAlternative(REF, ALT, (short) AD));
                    }
                }

                // Add the reference allele if no alternatives are present or if it has read support.
                if (alternatives.isEmpty() || referenceCall.depth() > 0) {
                    alternatives.add(referenceCall);
                }

                // Add the variant call to the storage and handle the result.
                processVariantCall(sampleIdentifier, contigIdentifier, variantContext.getStart(), alternatives, storage.parameters, path);
            }
        }
    }

    /**
     * Processes a variant call for a specific sample, contig, and position.
     * <p>
     * This method handles the processing of variant calls by merging new alternatives with existing ones, calculating coverage depth,
     * entropy, and frequency, and determining the appropriate flag for the variant call. It also manages special cases such as missing
     * alleles and upstream deletions.
     *
     * @param sampleIdentifier The identifier of the sample associated with the variant call.
     * @param contigIdentifier The identifier of the contig where the variant call is located.
     * @param position         The genomic position of the variant call.
     * @param alternatives     A list of {@link VariantCall.CallAlternative} objects representing the alleles.
     * @param parameters       The {@link Storage.Parameters} object containing processing parameters.
     * @param origin           The {@link Path} to the source file for logging purposes.
     */
    private void processVariantCall(String sampleIdentifier, String contigIdentifier, int position,
                                    List<VariantCall.CallAlternative> alternatives, Storage.Parameters parameters,
                                    Path origin) {
        // Ensure that the list of alternatives is not empty.
        assert !alternatives.isEmpty();

        // Retrieve existing alternatives and merge with new ones.
        VariantCall existingCall = vcCache.get(sampleIdentifier, contigIdentifier, position);
        if (existingCall != null) {
            List<VariantCall.CallAlternative> existingAlternatives = existingCall.alternatives();

            for (VariantCall.CallAlternative alternative : alternatives) {
                int index = existingAlternatives.indexOf(alternative);

                if (index >= 0) {
                    // Update allelic depth for existing alternative.
                    VariantCall.CallAlternative existing = existingAlternatives.get(index);
                    existingAlternatives.set(index, new VariantCall.CallAlternative(
                            existing.reference(),
                            existing.alternative(),
                            (short) (existing.depth() + alternative.depth())
                    ));
                } else {
                    // Add new alternative to the list.
                    existingAlternatives.add(alternative);
                }
            }

            // Replace alternatives with the merged list.
            alternatives = existingAlternatives;
        }

        // Sort alleles in descending order by their allelic depth (AD).
        alternatives.sort((a, b) -> Short.compare(b.depth(), a.depth()));

        // Calculate the total observed depth of coverage.
        short depth = (short) alternatives.stream().mapToInt(VariantCall.CallAlternative::depth).sum();

        // Calculate the normalized entropy of the variant call.
        float entropy = 0;
        if (alternatives.size() > 1) {
            // Iterate through each alternative allele to compute entropy.
            for (VariantCall.CallAlternative alternative : alternatives) {
                float frequency = (float) alternative.depth() / depth;
                entropy += frequency * Math.log10(frequency);
            }
            // Normalize the entropy value based on the number of non-zero depth alleles.
            entropy = (float) (-entropy / Math.log10(alternatives.size())) + (float) 0.0;
        }

        // Access the allele with the highest depth of coverage.
        VariantCall.CallAlternative allele = alternatives.get(0);
        VariantCall.Flag flag = allele.alternative().equals(Constants.DOT) ? VariantCall.Flag.REFERENCE_CALL : VariantCall.Flag.PASS;

        // Compute the actual frequency of the selected allele.
        float frequency = allele.depth() / (float) depth;

        // Set call prefix for low frequency or coverage.
        if (frequency < parameters.minimalFrequency()) flag = VariantCall.Flag.LOW_FREQUENCY;
        if (depth < parameters.minimalCoverage()) flag = VariantCall.Flag.LOW_COVERAGE;
        boolean isFiltered = (flag.equals(VariantCall.Flag.LOW_FREQUENCY) || flag.equals(VariantCall.Flag.LOW_COVERAGE));

        // Handle special cases.
        if (allele.alternative().equals("*")) {
            // Validate missing allele against upstream deletion.
            UpstreamDeletion upstreamDeletion = this.upstreamDeletions.get(sampleIdentifier);
            if (upstreamDeletion == null ||
                    (upstreamDeletion.contigIdentifier().equals(contigIdentifier) && position > upstreamDeletion.end())) {
                Logging.logWarningOnce("MISSING_DELETION",
                        String.format("Erroneous data with missing allele (*) not explained by an upstream deletion (%s:g.%d %s in %s).",
                                contigIdentifier, position, sampleIdentifier, origin));
                flag = VariantCall.Flag.MISSING_UPSTREAM_DELETION;
            } else if (upstreamDeletion.contigIdentifier().equals(contigIdentifier) && upstreamDeletion.start() <= position) {
                if (upstreamDeletion.isFiltered()) {
                    Logging.logWarningOnce("UNEXPLAINED_DELETION",
                            String.format("Ambiguous data with missing allele (*) and filtered upstream deletion (%s:g.%d %s in %s).",
                                    contigIdentifier, position, sampleIdentifier, origin));
                    flag = VariantCall.Flag.MISSING_UPSTREAM_DELETION;
                } else {
                    flag = VariantCall.Flag.UPSTREAM_DELETION;
                }
            }
        } else if (Bio.isDeletion(allele.alternative())) {
            // Set upstream deletion, if the current accepted call is a deletion.
            this.upstreamDeletions.put(sampleIdentifier, new UpstreamDeletion(contigIdentifier,
                    position + StringUtils.indexOf(allele.alternative(), Constants.GAP_CHAR),
                    position + StringUtils.lastIndexOf(allele.alternative(), Constants.GAP_CHAR), isFiltered)
            );
        }

        // Store the variant call in the cache if it is not a reference call or missing due to an upstream deletion.
        if (flag != VariantCall.Flag.REFERENCE_CALL && flag != VariantCall.Flag.UPSTREAM_DELETION) {
            vcCache.put(sampleIdentifier, contigIdentifier, position, new VariantCall(flag, depth, entropy, alternatives));
        }

        // Update statistics based on the flag.
        switch (flag) {
            case PASS, UPSTREAM_DELETION -> {
                // No action needed.
            }
            case REFERENCE_CALL -> ignoredCallsCount++;
            case LOW_COVERAGE, LOW_FREQUENCY, MISSING_UPSTREAM_DELETION -> filteredCallsCount++;
            default -> throw new IllegalStateException("Unexpected state " + flag);
        }
    }

}
