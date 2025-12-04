package op;

import main.Musial;
import model.*;
import tech.tablesaw.api.*;
import tech.tablesaw.io.csv.CsvWriteOptions;
import util.Constants;
import util.Logging;

import java.util.*;
import java.util.stream.Collectors;

/**
 * The {@code StorageTable} class provides functionality to create and manage tabular representations of genomic data stored in a
 * {@link Storage} instance using the <a href="https://github.com/jtablesaw/tablesaw">tablesaw</a> library.
 */
public class StorageTable {

    /**
     * The storage object for managing genomic data.
     */
    private final Storage storage;

    /**
     * The {@link Table} object representing the data table.
     */
    private Table table;

    /**
     * Filters (inclusive {@link model.Contig#_id}s and associated {@link model.Variant#position}s) for building the table.
     */
    private final Map<String, int[]> contigFilter;

    /**
     * Filters (inclusive {@link model.Feature#_id}s or {@link model.Feature#name}s) for building the table.
     */
    private final Set<String> featureFilter;

    /**
     * Filters (inclusive {@link model.Sample#_id}s) for building the table.
     */
    private final Set<String> sampleFilter;

    /**
     * Constructs an instance of the {@link StorageTable} class with the specified storage.
     *
     * @param storage The {@link Storage} instance containing genomic data.
     */
    public StorageTable(Storage storage) {
        this.storage = storage;
        this.contigFilter = new HashMap<>(storage.getContigs().size());
        this.featureFilter = new HashSet<>(storage.getFeatures().size());
        this.sampleFilter = new HashSet<>(storage.getSamples().size());
    }

    /**
     * Clears all filters used for building the table.
     * <p>
     * This method removes all entries from the contig, feature, and sample filters. After calling this method, the filters will be empty,
     * and the user must reapply filters if needed before populating the table.
     */
    public void clearFilters() {
        contigFilter.clear();
        featureFilter.clear();
        sampleFilter.clear();
    }

    /**
     * Updates the contig filter with the specified contig identifier and positions.
     * <p>
     * If the contig identifier does not exist in the storage, no action is taken and {@code false} is returned.
     * <p>
     * If the contig identifier already exists in the filter, its associated positions are replaced. An empty array of positions indicates
     * that all positions within the contig should be included. If the filter is successfully updated, {@code true} is returned.
     * <p>
     * The table is not automatically updated after modifying the filter; the user must call one of the populate methods to refresh the
     * table.
     *
     * @param contigIdentifier The identifier of the contig to be added or updated in the filter.
     * @param positions        An array of positions (as integers) associated with the contig. If the array is empty, it indicates that all
     *                         positions within the contig should be included.
     * @return {@code true} if the filter was successfully updated; {@code false} if the contig identifier does not exist in the storage.
     */
    public boolean addContigFilter(String contigIdentifier, int[] positions) {
        if (storage.hasContig(contigIdentifier)) {
            contigFilter.put(contigIdentifier, positions);
            return true;
        } else {
            return false;
        }
    }

    /**
     * Adds a feature filter using the provided identifier or name.
     * <p>
     * This method checks if the given feature identifier or name exists in the storage. If it exists, the feature is added to the
     * {@link #featureFilter}. If the feature does not exist, no changes are made, and {@code false} is returned.
     * <p>
     * The table is not automatically updated after modifying the filter; the user must call one of the populate methods to refresh the
     * table.
     *
     * @param featureIdentifierOrName The identifier or name of the feature to include in the filter.
     * @return {@code true} if the feature was successfully added to the filter; {@code false} if the feature does not exist in the storage.
     */
    public boolean addFeatureFilter(String featureIdentifierOrName) {
        if (storage.hasFeature(featureIdentifierOrName)) {
            featureFilter.add(featureIdentifierOrName);
            return true;
        } else {
            Optional<Feature> optionalFeature = storage.getFeatures().stream()
                    .filter(feature -> feature.name.equals(featureIdentifierOrName))
                    .findFirst();
            if (optionalFeature.isPresent()) {
                featureFilter.add(optionalFeature.get()._id);
                return true;
            } else {
                return false;
            }
        }
    }

    /**
     * Adds a sample filter using the provided identifier.
     * <p>
     * This method checks if the given sample identifier exists in the storage. If it exists, the sample is added to the
     * {@link #sampleFilter}. If the sample does not exist, no changes are made, and {@code false} is returned.
     * <p>
     * The table is not automatically updated after modifying the filter; the user must call one of the populate methods to refresh the
     * table.
     *
     * @param sampleIdentifier The identifier of the sample to include in the filter.
     * @return {@code true} if the sample was successfully added to the filter; {@code false} if the sample does not exist in the storage.
     */
    public boolean addSampleFilter(String sampleIdentifier) {
        if (storage.hasSample(sampleIdentifier)) {
            sampleFilter.add(sampleIdentifier);
            return true;
        } else {
            return false;
        }
    }

    /**
     * Sets filters for building the table based on the provided query strings.
     * <p>
     * This method clears all existing filters and processes each query string to determine the type of filter to apply. Queries can specify
     * sample filters, feature filters, or contig filters with optional position ranges.
     * <p>
     * The method performs the following steps for each query:
     * <ul>
     *     <li>Trims the query string and skips blank queries.</li>
     *     <li>Attempts to add the query as a sample or feature filter.</li>
     *     <li>If not a sample or feature, parses the query as a contig filter with optional position ranges.</li>
     *     <li>Handles position ranges in the format "start-end" or "start" and generates an array of positions.</li>
     * </ul>
     * <p>
     * If the query cannot be parsed (e.g., due to invalid formatting), it is ignored.
     *
     * @param query A set of query strings specifying the filters to apply.
     */
    public void setFilters(Set<String> query) {
        clearFilters();
        for (String q : query) {
            if (q.isBlank()) continue;

            q = q.trim();
            if (addSampleFilter(q) || addFeatureFilter(q)) continue;

            try {
                String[] parts = q.split(":", 2);
                String contig = parts[0];
                int[] positions = {};

                if (parts.length == 2) {
                    String[] range = parts[1].split("-", 2);
                    int start = Integer.parseInt(range[0]);
                    int end = (range.length == 2) ? Integer.parseInt(range[1]) : start;

                    if (start > 0 && end >= start) {
                        positions = java.util.stream.IntStream.rangeClosed(start, end).toArray();
                    }
                }

                if (addContigFilter(contig, positions)) continue;

                Logging.logWarningOnce("UNKNOWN_QUERY", "Unable to match query %s to storage entries.".formatted(q));
            } catch (NumberFormatException | ArrayIndexOutOfBoundsException ignored) {
                // Ignored.
            }
        }
    }

    /**
     * Populates the {@link #table} with data of {@link Variant}s from the storage.
     * <p>
     * The table is created with predefined columns and filled with variant data filtered by the current contig, position, sample, and
     * feature filters. The table is then sorted in ascending order based on chromosome and position.
     * <p>
     * The columns of the table include:
     * <ul>
     *     <li>chrom: Chromosome/contig identifier.</li>
     *     <li>position: Variant position.</li>
     *     <li>reference: Reference allele/base content.</li>
     *     <li>alternative: Alternative allele/base content.</li>
     *     <li>type: Variant type, see {@link Variant.Type}.</li>
     *     <li>feature_id: Affected feature name or {@code null}.</li>
     *     <li>impact: Impact of the variant or {@code null}.</li>
     *     <li>effect: Effect of the variant or {@code null}.</li>
     *     <li>{@link Constants.AttributesKeys#VARIANT_FREQUENCY}: Variant frequency wrt. samples.</li>
     *     <li>{@link Constants.AttributesKeys#MEAN_COVERAGE}: Mean coverage wrt. samples.</li>
     *     <li>{@link Constants.AttributesKeys#MEAN_ENTROPY}: Mean entropy wrt. samples.</li>
     * </ul>
     */
    public void populateFromVariants() {
        // Initialize the table with predefined columns
        table = Table.create("Variants").addColumns(
                StringColumn.create("chrom"),
                IntColumn.create("position"),
                StringColumn.create("reference"),
                StringColumn.create("alternative"),
                StringColumn.create("type"),
                StringColumn.create(Constants.SNP_EFF_KEYS.get(6)),
                StringColumn.create(Constants.SNP_EFF_KEYS.get(2)),
                StringColumn.create(Constants.SNP_EFF_KEYS.get(1)),
                DoubleColumn.create(Constants.AttributesKeys.VARIANT_FREQUENCY),
                DoubleColumn.create(Constants.AttributesKeys.MEAN_COVERAGE),
                DoubleColumn.create(Constants.AttributesKeys.MEAN_ENTROPY)
        );

        // Filter contigs based on the contig filter
        Collection<Contig> filteredContigs = contigFilter.isEmpty()
                ? storage.getContigs()
                : storage.getContigs().stream()
                .filter(c -> contigFilter.containsKey(c._id))
                .toList();

        // Process each contig and its variants
        for (Contig contig : filteredContigs) {
            int[] positions = contigFilter.getOrDefault(contig._id, new int[]{});
            List<Variant> variants = positions.length == 0
                    ? contig.getAllVariants()
                    : contig.getVariantsAt(positions);

            for (Variant variant : variants) {
                // Skip variants that do not match the filters
                if (!sampleFilter.isEmpty() && !variant.ofSamples(sampleFilter)) continue;
                if (!featureFilter.isEmpty() && !variant.ofFeatures(featureFilter)) continue;

                // Add variant data to the table
                Row row = table.appendRow();
                row.setString("chrom", contig._id);
                row.setInt("position", variant.position);
                row.setString("reference", variant.reference);
                row.setString("alternative", variant.alternative);
                row.setString("type", variant.type.toString().toUpperCase());

                // Set feature-related attributes
                String featureId = variant.getAttribute(Constants.SNP_EFF_PREFIX + Constants.SNP_EFF_KEYS.get(6));
                if (storage.hasFeature(featureId)) {
                    row.setString(Constants.SNP_EFF_KEYS.get(6), storage.getFeature(featureId).name);
                    row.setString(Constants.SNP_EFF_KEYS.get(2),
                            variant.getAttribute(Constants.SNP_EFF_PREFIX + Constants.SNP_EFF_KEYS.get(2)));
                    row.setString(Constants.SNP_EFF_KEYS.get(1),
                            variant.getAttribute(Constants.SNP_EFF_PREFIX + Constants.SNP_EFF_KEYS.get(1)));
                } else {
                    row.setString(Constants.SNP_EFF_KEYS.get(6), null);
                    row.setString(Constants.SNP_EFF_KEYS.get(2), null);
                    row.setString(Constants.SNP_EFF_KEYS.get(1), null);
                }

                // Set additional attributes
                row.setDouble(Constants.AttributesKeys.VARIANT_FREQUENCY,
                        Double.parseDouble(variant.getAttributeOrDefault(Constants.AttributesKeys.VARIANT_FREQUENCY, "0")));
                row.setDouble(Constants.AttributesKeys.MEAN_COVERAGE,
                        Double.parseDouble(variant.getAttributeOrDefault(Constants.AttributesKeys.MEAN_COVERAGE, "0")));
                row.setDouble(Constants.AttributesKeys.MEAN_ENTROPY,
                        Double.parseDouble(variant.getAttributeOrDefault(Constants.AttributesKeys.MEAN_ENTROPY, "0")));
            }
        }

        // Sort the table by chromosome and position
        table.sortAscendingOn("chrom", "position");
    }

    /**
     * Populates the {@link #table} with data of {@link model.Feature}s from the storage.
     * <p>
     * The table is created with predefined columns and filled with feature data filtered by the current contig and feature filters.
     * Additional attributes are dynamically added as columns if they are not part of the standard attributes. The table is then sorted in
     * ascending order based on chromosome and start position.
     * <p>
     * The (default) columns of the table include:
     * <ul>
     *     <li>chrom: Chromosome/contig identifier.</li>
     *     <li>start: Start position of the feature.</li>
     *     <li>end: End position of the feature.</li>
     *     <li>strand: Strand information (+/-).</li>
     *     <li>id: Feature identifier.</li>
     *     <li>name: Feature name.</li>
     *     <li>type: Feature type.</li>
     *     <li>{@link Constants.AttributesKeys#NUMBER_OF_ALLELES}: Number of non-reference alleles associated with the feature.</li>
     *     <li>{@link Constants.AttributesKeys#DIVERSITY_ALLELE}: Diversity of alleles associated with the feature.</li>
     *     <li>{@link Constants.AttributesKeys#FREQUENCY_REFERENCE}: Fraction of samples with the reference allele wrt. the feature.</li>
     *     <li>{@link Constants.AttributesKeys#NUMBER_OF_PROTEOFORMS}: Number of non-reference proteoforms associated with the feature (0
     *     for non-coding features).</li>
     *     <li>{@link Constants.AttributesKeys#DIVERSITY_PROTEOFORM}: Diversity of proteoforms associated with the feature (0.0 for
     *     non-coding features).</li>
     *     <li>{@link Constants.AttributesKeys#FREQUENCY_DISRUPTED}: Fraction of samples with disrupted proteoforms associated with the
     *     feature (0.0 for non-coding features).</li>
     * </ul>
     */
    public void populateFromFeatures() {
        // Create a new table named "Features" with predefined columns.
        table = Table.create("Features").addColumns(
                StringColumn.create("chrom"),
                IntColumn.create("start"),
                IntColumn.create("end"),
                StringColumn.create("strand"),
                StringColumn.create("id"),
                StringColumn.create("name"),
                StringColumn.create("type"),
                IntColumn.create(Constants.AttributesKeys.NUMBER_OF_ALLELES),
                DoubleColumn.create(Constants.AttributesKeys.DIVERSITY_ALLELE),
                DoubleColumn.create(Constants.AttributesKeys.FREQUENCY_REFERENCE),
                IntColumn.create(Constants.AttributesKeys.NUMBER_OF_PROTEOFORMS),
                DoubleColumn.create(Constants.AttributesKeys.DIVERSITY_PROTEOFORM),
                DoubleColumn.create(Constants.AttributesKeys.FREQUENCY_DISRUPTED)
        );

        // Add custom attributes as columns if they are not standard attributes.
        Set<String> standardAttributes = Set.of(
                Constants.AttributesKeys.NUMBER_OF_ALLELES, Constants.AttributesKeys.DIVERSITY_ALLELE,
                Constants.AttributesKeys.FREQUENCY_REFERENCE,
                Constants.AttributesKeys.NUMBER_OF_PROTEOFORMS, Constants.AttributesKeys.DIVERSITY_PROTEOFORM,
                Constants.AttributesKeys.FREQUENCY_DISRUPTED
        );
        storage.getFeatureAttributeKeys()
                .stream()
                .filter(key -> !standardAttributes.contains(key))
                .forEach(key -> table.addColumns(StringColumn.create(key)));

        // Iterate over features and populate the table.
        for (var feature : storage.getFeatures()) {
            if ((contigFilter.isEmpty() || contigFilter.containsKey(feature.contig))
                    && (featureFilter.isEmpty() || featureFilter.contains(feature._id) || featureFilter.contains(feature.name))) {
                Row row = table.appendRow();
                row.setString("chrom", feature.contig);
                row.setInt("start", feature.start);
                row.setInt("end", feature.end);
                row.setString("strand", String.valueOf(feature.strand));
                row.setString("id", feature._id);
                row.setString("name", feature.name);
                row.setString("type", feature.type.toLowerCase());

                // Set standard attributes.
                row.setInt(Constants.AttributesKeys.NUMBER_OF_ALLELES, feature.getAlleleCount());
                row.setDouble(Constants.AttributesKeys.DIVERSITY_ALLELE,
                        Double.parseDouble(feature.getAttribute(Constants.AttributesKeys.DIVERSITY_ALLELE)));
                row.setDouble(Constants.AttributesKeys.FREQUENCY_REFERENCE,
                        Double.parseDouble(feature.getAttribute(Constants.AttributesKeys.FREQUENCY_REFERENCE)));

                // Set coding-specific attributes or defaults for non-coding features.
                if (feature.isCoding()) {
                    row.setInt(Constants.AttributesKeys.NUMBER_OF_PROTEOFORMS, feature.getProteoformCount());
                    row.setDouble(Constants.AttributesKeys.DIVERSITY_PROTEOFORM,
                            Double.parseDouble(feature.getAttribute(Constants.AttributesKeys.DIVERSITY_PROTEOFORM)));
                    row.setDouble(Constants.AttributesKeys.FREQUENCY_DISRUPTED,
                            Double.parseDouble(feature.getAttribute(Constants.AttributesKeys.FREQUENCY_DISRUPTED)));
                } else {
                    row.setInt(Constants.AttributesKeys.NUMBER_OF_PROTEOFORMS, 0);
                    row.setDouble(Constants.AttributesKeys.DIVERSITY_PROTEOFORM, 0.0);
                    row.setDouble(Constants.AttributesKeys.FREQUENCY_DISRUPTED, 0.0);
                }

                // Add additional attributes.
                feature.getAttributes().forEach((key, value) -> {
                    if (!standardAttributes.contains(key)) {
                        if (!table.containsColumn(key)) {
                            table.addColumns(StringColumn.create(key));
                        }
                        row.setString(key, value);
                    }
                });
            }
        }

        // Sort the table by chromosome and start position.
        table.sortAscendingOn("chrom", "start");
    }

    /**
     * Populates the {@link #table} with data of {@link model.Feature}s from the storage.
     * <p>
     * The table is created with predefined columns and filled with sample data filtered by the current sample filter. Additional attributes
     * are dynamically added as columns if they are not part of the standard attributes. The table is then sorted in ascending order based
     * on diversity allele and diversity proteoform.
     * <p>
     * The (default) columns of the table include:
     */
    public void populateFromSamples() {
        // Create the "Samples" table with predefined columns.
        table = Table.create("Samples").addColumns(
                StringColumn.create("id"),
                IntColumn.create(Constants.AttributesKeys.NUMBER_OF_SNVS),
                IntColumn.create(Constants.AttributesKeys.NUMBER_OF_INDELS),
                DoubleColumn.create(Constants.AttributesKeys.FREQUENCY_FILTERED_CALLS),
                DoubleColumn.create(Constants.AttributesKeys.MEAN_COVERAGE),
                DoubleColumn.create(Constants.AttributesKeys.FREQUENCY_REFERENCE),
                DoubleColumn.create(Constants.AttributesKeys.FREQUENCY_DISRUPTED)
        );

        // Add custom attributes as columns if they are not standard attributes.
        Set<String> standardAttributes = Set.of(
                Constants.AttributesKeys.NUMBER_OF_SNVS, Constants.AttributesKeys.NUMBER_OF_INDELS,
                Constants.AttributesKeys.FREQUENCY_FILTERED_CALLS,
                Constants.AttributesKeys.MEAN_COVERAGE, Constants.AttributesKeys.FREQUENCY_REFERENCE,
                Constants.AttributesKeys.FREQUENCY_DISRUPTED
        );
        storage.getSampleAttributeKeys()
                .stream()
                .filter(key -> !standardAttributes.contains(key))
                .forEach(key -> table.addColumns(StringColumn.create(key)));

        // Populate the table with sample data.
        storage.getSamples().stream()
                .filter(sample -> sampleFilter.isEmpty() || sampleFilter.contains(sample._id))
                .forEach(sample -> {
                    Row row = table.appendRow();
                    row.setString("id", sample._id);
                    row.setInt(Constants.AttributesKeys.NUMBER_OF_SNVS,
                            Integer.parseInt(sample.getAttribute(Constants.AttributesKeys.NUMBER_OF_SNVS)));
                    row.setInt(Constants.AttributesKeys.NUMBER_OF_INDELS,
                            Integer.parseInt(sample.getAttribute(Constants.AttributesKeys.NUMBER_OF_INDELS)));
                    row.setDouble(Constants.AttributesKeys.FREQUENCY_FILTERED_CALLS,
                            Double.parseDouble(sample.getAttribute(Constants.AttributesKeys.FREQUENCY_FILTERED_CALLS)));
                    row.setDouble(Constants.AttributesKeys.MEAN_COVERAGE,
                            Double.parseDouble(sample.getAttribute(Constants.AttributesKeys.MEAN_COVERAGE)));
                    row.setDouble(Constants.AttributesKeys.FREQUENCY_REFERENCE,
                            Double.parseDouble(sample.getAttribute(Constants.AttributesKeys.FREQUENCY_REFERENCE)));
                    row.setDouble(Constants.AttributesKeys.FREQUENCY_DISRUPTED,
                            Double.parseDouble(sample.getAttribute(Constants.AttributesKeys.FREQUENCY_DISRUPTED)));

                    // Add additional attributes.
                    sample.getAttributes().entrySet().stream()
                            .filter(attr -> !standardAttributes.contains(attr.getKey()))
                            .forEach(attr -> {
                                if (!table.containsColumn(attr.getKey())) {
                                    table.addColumns(StringColumn.create(attr.getKey()));
                                }
                                row.setString(attr.getKey(), attr.getValue());
                            });
                });

        // Sort the table by the number of SNVs and INDELs.
        table.sortAscendingOn(Constants.AttributesKeys.NUMBER_OF_SNVS, Constants.AttributesKeys.NUMBER_OF_INDELS);
    }

    /**
     * Populates the {@link #table} with allele profile data from the storage.
     * <p>
     * Each row in the table represents a sample, with the first column being the sample identifier ("id") and subsequent columns
     * representing features. The values in the feature columns are allele identifiers or numeric indices (0 represents the reference, 1 the
     * first non-reference allele, and so on), depending on the `enumerate` flag.
     * <p>
     * The data is filtered by the current contig and feature filters.
     *
     * @param enumerate A boolean flag indicating whether to replace allele identifiers with numeric indices.
     */
    public void populateWithAlleleProfile(boolean enumerate) {
        // Create the "AlleleProfile" table with an "id" column.
        table = Table.create("AlleleProfile").addColumns(StringColumn.create("id"));

        // Initialize feature and identifier maps.
        Map<String, Feature> featureMap = new HashMap<>();
        Map<String, Map<String, Integer>> identifierMap = new HashMap<>();

        // Filter and process features.
        storage.getFeatures().stream()
                .filter(feature -> (contigFilter.isEmpty() || contigFilter.containsKey(feature.contig))
                        && (featureFilter.isEmpty() || featureFilter.contains(feature.name) || featureFilter.contains(feature._id)))
                .forEach(feature -> {
                    // Add the feature to the feature map.
                    featureMap.put(feature.name, feature);
                    if (enumerate) {
                        // Initialize the identifier map for the feature if enumeration is enabled.
                        identifierMap.put(feature.name, new HashMap<>());
                        identifierMap.get(feature.name).put(Constants.REFERENCE, 0);
                    }
                });

        // Add feature names as columns to the table.
        featureMap.keySet().forEach(name -> table.addColumns(StringColumn.create(name)));

        // Populate the table with sample data.
        storage.getSamples().stream()
                .filter(sample -> sampleFilter.isEmpty() || sampleFilter.contains(sample._id))
                .forEach(sample -> {
                    // Append a new row for the sample.
                    Row row = table.appendRow();
                    row.setString("id", sample._id);

                    // Populate the row with allele data for each feature.
                    featureMap.forEach((name, feature) -> {
                        String alleleIdentifier = sample.getRelatedAllele(feature._id);
                        if (enumerate) {
                            // Replace allele identifiers with numeric indices if enumeration is enabled.
                            identifierMap.get(name).putIfAbsent(alleleIdentifier, identifierMap.get(name).size());
                            row.setString(name, String.valueOf(identifierMap.get(name).get(alleleIdentifier)));
                        } else {
                            // Use the allele identifier directly.
                            row.setString(name, alleleIdentifier);
                        }
                    });
                });
    }

    /**
     * Populates the {@link #table} with proteoform profile data from the storage.
     * <p>
     * Each row in the table represents a sample, with the first column being the sample identifier ("id") and subsequent columns
     * representing coding features. The values in the feature columns are proteoform identifiers or numeric indices (0 represents
     * synonymous, 1 the first non-synonymous proteoform, and so on), depending on the `enumerate` flag.
     * <p>
     * The data is filtered by the current contig and feature filters, and only coding features are included.
     *
     * @param enumerate A boolean flag indicating whether to replace proteoform identifiers with numeric indices.
     */
    public void populateWithProteoformProfile(boolean enumerate) {
        // Initialize the "ProteoformProfile" table with an "id" column.
        table = Table.create("ProteoformProfile").addColumns(StringColumn.create("id"));

        // Filter and map coding features.
        Map<String, Feature> featureMap = storage.getFeatures().stream()
                .filter(feature -> (contigFilter.isEmpty() || contigFilter.containsKey(feature.contig))
                        && (featureFilter.isEmpty() || featureFilter.contains(feature.name) || featureFilter.contains(feature._id))
                        && feature.isCoding())
                .collect(Collectors.toMap(feature -> feature.name, feature -> feature));

        // Initialize identifier map if enumeration is enabled.
        Map<String, Map<String, Integer>> identifierMap = new HashMap<>();
        if (enumerate) {
            featureMap.forEach((name, feature) -> {
                Map<String, Integer> map = new HashMap<>();
                map.put(Constants.SYNONYMOUS, 0);
                identifierMap.put(name, map);
            });
        }

        // Add feature names as columns to the table.
        featureMap.keySet().forEach(name -> table.addColumns(StringColumn.create(name)));

        // Populate the table with sample data.
        storage.getSamples().stream()
                .filter(sample -> sampleFilter.isEmpty() || sampleFilter.contains(sample._id))
                .forEach(sample -> {
                    Row row = table.appendRow();
                    row.setString("id", sample._id);

                    featureMap.forEach((name, feature) -> {
                        String alleleIdentifier = sample.getRelatedAllele(feature._id);
                        String proteoformIdentifier = alleleIdentifier.equals(Constants.REFERENCE)
                                ? Constants.SYNONYMOUS
                                : feature.getAllele(alleleIdentifier).getRelatedProteoform();

                        if (enumerate) {
                            identifierMap.get(name).putIfAbsent(proteoformIdentifier, identifierMap.get(name).size());
                            row.setString(name, String.valueOf(identifierMap.get(name).get(proteoformIdentifier)));
                        } else {
                            row.setString(name, proteoformIdentifier);
                        }
                    });
                });
    }

    /**
     * Populates the {@link #table} with variant profile data from the storage.
     * <p>
     * Each row in the table represents a unique variant position, with columns for chromosome, position, reference allele, and samples. The
     * values in the sample columns are either the alternative allele or a placeholder (e.g., ".") if the sample does not have the variant.
     * The `simple` flag determines whether to use a simplified representation of the variant calls.
     * <p>
     * The data is filtered by the current contig, position, and sample filters.
     *
     * @param simple A boolean flag indicating whether to use a simplified representation of variant calls.
     */
    public void populateWithVariantsProfile(boolean simple) {
        // Create the "VariantsProfile" table with predefined columns.
        table = Table.create("VariantsProfile").addColumns(
                StringColumn.create("chrom"),
                IntColumn.create("position"),
                StringColumn.create("reference")
        );

        // Add sample columns to the table.
        Collection<Sample> samples = sampleFilter.isEmpty()
                ? storage.getSamples()
                : storage.getSamples().stream()
                .filter(s -> sampleFilter.contains(s._id))
                .toList();
        samples.forEach(sample -> table.addColumns(StringColumn.create(sample._id)));

        // Filter contigs based on the contig filter.
        Collection<Contig> filteredContigs = contigFilter.isEmpty()
                ? storage.getContigs()
                : storage.getContigs().stream()
                .filter(contig -> contigFilter.containsKey(contig._id))
                .toList();

        int currentPosition = 0;
        Row currentRow = null;

        // Iterate over filtered contigs and their variants.
        for (Contig contig : filteredContigs) {
            int[] positions = contigFilter.getOrDefault(contig._id, new int[]{});
            List<Variant> variants = positions.length == 0
                    ? contig.getAllVariants()
                    : contig.getVariantsAt(positions);

            for (Variant variant : variants) {

                // Skip variants that do not match the sample or feature filters.
                if ((!sampleFilter.isEmpty() && sampleFilter.stream().noneMatch(variant::ofSample))
                        || (!featureFilter.isEmpty() && featureFilter.stream().noneMatch(variant::ofFeature))) continue;

                // Create a new row if the position changes.
                if (variant.position != currentPosition) {
                    currentRow = table.appendRow();
                    currentRow.setString("chrom", contig._id);
                    currentRow.setInt("position", variant.position);
                    currentRow.setString("reference", String.valueOf(variant.reference.charAt(0)));
                    for (Sample sample : samples) {
                        currentRow.setString(sample._id, Constants.DOT);
                    }
                }
                currentPosition = variant.position;

                // Update the row with sample-specific data.
                assert currentRow != null;
                for (var relation : variant.getRelatedSamples()) {
                    String value = simple
                            ? (VariantCall.isFiltered(relation.b) ? Constants.ANY_NUCLEOTIDE : variant.alternative)
                            : relation.b;
                    currentRow.setString(relation.a, value);
                }
            }
        }
    }

    /**
     * Writes the current {@link #table} to a tabular file with the specified separator.
     * <p>
     * If the table is not populated (i.e., it is {@code null} or empty), an {@link IllegalStateException} is thrown.
     *
     * @param file      The path to the output file.
     * @param separator The character to use as the separator in the file.
     */
    public void write(String file, char separator) {
        if (Objects.isNull(table))
            throw new IllegalStateException("The table is not populated. Please call one of the populate methods before writing.");
        if (table.isEmpty())
            Logging.logWarning("The table is empty.");
        table.write().csv(CsvWriteOptions.builder(file).separator(separator).build());
    }

    /**
     * Writes the current {@link #table} to a TSV file in the current user directory at {@code System.getProperty("user.dir")}.
     * <p>
     * If the table is not populated (i.e., it is {@code null} or empty), an {@link IllegalStateException} is thrown.
     */
    public void write() {
        if (Objects.isNull(table))
            throw new IllegalStateException("The table is not populated. Please call one of the populate methods before writing.");
        if (table.isEmpty())
            Logging.logWarning("The table is empty.");
        write(Musial.tempDir.getParentFile().getAbsolutePath() + "/" + this.table.name() + ".tsv", '\t');
    }

    /**
     * Prints the contents of the {@link #table} to the standard output.
     * <p>
     * This method checks if the {@link #table} is populated and not empty before printing. If the table is not populated, an
     * {@link IllegalStateException} is thrown. Depending on the value of the {@code all} parameter, it either prints the entire table or a
     * truncated version.
     *
     * @param all A boolean flag indicating whether to print the entire table ({@code true}) or a truncated version ({@code false}).
     * @throws IllegalStateException If the table is not populated or is empty.
     */
    public void print(boolean all) {
        if (Objects.isNull(table))
            throw new IllegalStateException("The table is not populated. Please call one of the populate methods before printing.");
        if (table.isEmpty())
            Logging.logWarning("The table is empty.");
        if (all) {
            table.write().csv(System.out);
        } else {
            System.out.println(table.print());
        }
    }

}
