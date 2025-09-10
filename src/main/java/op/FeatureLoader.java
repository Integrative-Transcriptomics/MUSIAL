package op;

import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import model.Feature;
import model.Storage;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import util.Constants;
import util.Logging;

import java.util.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;

/**
 * The FeatureLoader class is responsible for loading and validating genomic features.
 * <p>
 * This class provides methods to load features from a GFF3 annotation file into the storage system, validate the features against Sequence
 * Ontology (SO) hierarchy rules, and adjust or impute features as necessary. It interacts with the {@link Storage} object to manage
 * contigs, features, samples, and variant calls.
 * <p>
 * The FeatureLoader is initialized with a {@link Storage} object, a {@link FeatureList} containing parsed features, and a map of
 * user-defined feature specifications. It ensures that the genomic data is processed and stored in a consistent and hierarchical manner.
 */
public class FeatureLoader {

    /**
     * The storage object for managing genomic data.
     */
    private final Storage storage;

    /**
     * The list of features parsed from the GFF3 annotation file.
     * <p>
     * This field contains a collection of features that are extracted from the GFF3 file. Each feature represents a genomic element with
     * associated attributes such as type, location, and hierarchy.
     */
    private final FeatureList featureList;

    /**
     * A map of feature specifications provided by the user.
     * <p>
     * This field stores a mapping of feature identifiers to their corresponding specifications. Each specification is represented as a map
     * of key-value pairs that define the attributes and matching criteria for the feature.
     */
    private final Map<String, Map<String, String>> features;

    /**
     * Constructs a new instance of the {@link FeatureLoader} class.
     * <p>
     * This constructor initializes the {@link FeatureLoader} with the provided storage, feature list, and feature specifications. The
     * {@link FeatureLoader} is responsible for loading and validating genomic features based on the given data.
     *
     * @param storage     The {@link Storage} object used to manage genomic data, including contigs, features, samples, and variant calls.
     * @param featureList The {@link FeatureList} containing features parsed from the GFF3 annotation file.
     * @param features    A map of feature specifications provided by the user, where each key is a feature identifier and the value is a
     *                    map of attributes defining the feature's properties and matching criteria.
     */
    public FeatureLoader(Storage storage, FeatureList featureList, Map<String, Map<String, String>> features) {
        this.storage = storage;
        this.featureList = featureList;
        this.features = features;
    }

    /**
     * Loads features into the storage based on the provided annotations and specifications.
     * <p>
     * This method validates the input conditions to ensure that the reference sequence and annotations are properly specified. It processes
     * the features specified in the CLI or loads all annotated features from the provided GFF3 file.
     * <p>
     * The method performs the following steps:
     * <ul>
     *   <li>Validates the presence of a reference sequence and annotations.</li>
     *   <li>Processes attributes to ensure they are in the correct format.</li>
     *   <li>Matches and adds specified features to the storage.</li>
     *   <li>Loads all annotated features if no specific features are provided.</li>
     * </ul>
     *
     * @throws MusialException If the input conditions are invalid or if a feature is missing required specifications.
     */
    public void loadFeatures() throws MusialException {
        // Validate input conditions
        if (!storage.hasReference() && !featureList.isEmpty()) {
            throw new MusialException("Annotation (GFF3) specified without reference sequence (FASTA).");
        }
        if (featureList.isEmpty() && !features.isEmpty()) {
            throw new MusialException("Features specified without annotation (GFF3).");
        }

        // Helper to process attributes
        Consumer<Map<String, String>> reprocessAttributes = attributes ->
                attributes.replaceAll((k, v) -> Arrays.stream(v.split(Constants.COMMA))
                        .map(s -> s.split("\\|")[0])
                        .filter(s -> !s.isEmpty())
                        .collect(Collectors.joining(Constants.COMMA)));

        if (!features.isEmpty()) { // Process specified features.
            Logging.logConfig("Match %d specified features from annotation (GFF3).".formatted(features.size()));
            for (var feature : features.entrySet()) {
                String featureName = feature.getKey();
                Map<String, String> featureSpecification = feature.getValue();
                String matchKey = featureSpecification.remove("key");
                String matchValue = featureSpecification.remove("value");

                if (Objects.isNull(matchKey) || Objects.isNull(matchValue)) {
                    // Note: This should not happen, as the CLI parser already validated this.
                    throw new MusialException("Feature %s is missing key or value specification.".formatted(featureName));
                }

                FeatureList matchedFeatures = featureList.selectByAttribute(matchKey, matchValue);
                if (matchedFeatures.isEmpty()) {
                    Logging.logWarning("Failed to match feature by key `%s` and value `%s`.".formatted(matchKey, matchValue));
                    continue;
                }
                for (FeatureI matchedFeature : matchedFeatures) {
                    Map<String, String> attributes = matchedFeature.getAttributes();
                    attributes.putAll(featureSpecification);
                    reprocessAttributes.accept(attributes);
                    storage.addFeature(matchedFeature, featureName, attributes);
                }
            }
        } else if (!featureList.isEmpty()) { // Process all annotated features.
            Logging.logConfig("Load all %d annotated features from annotation (GFF3).".formatted(featureList.size()));
            for (FeatureI featureI : featureList) {
                if (!"region".equals(featureI.type())) {
                    Map<String, String> attributes = featureI.getAttributes();
                    reprocessAttributes.accept(attributes);
                    String name = attributes.getOrDefault("Name", "%s:g.%d_%d=".formatted(
                            featureI.seqname(), featureI.location().bioStart(), featureI.location().bioEnd()));
                    storage.addFeature(featureI, name, attributes);
                }
            }
        }
    }

    /**
     * Validates the features stored in the storage to ensure compliance with Sequence Ontology (SO) hierarchy rules.
     * <p>
     * This method iterates through all features in the storage and performs the following validations and adjustments:
     * <ul>
     *   <li>Removes children for features of level 0 SO term types, as they are not allowed to have children.</li>
     *   <li>Ensures that only one level 1 SO term exists for a feature or its children.</li>
     *   <li>Ensures that only one level 2 SO term exists for a feature or its children.</li>
     *   <li>Adjusts features with SO levels greater than 1 to type "gene" to maintain a consistent hierarchy.</li>
     *   <li>Imputes missing children based on the location ranges of existing sub-features.</li>
     * </ul>
     * <p>
     * Features that violate the rules are either adjusted or removed from the storage, and appropriate warnings are logged.
     *
     * @throws MusialException If an error occurs during the adjustment of features.
     */
    public void validateFeatures() throws MusialException {
        for (Feature feature : storage.getFeatures()) {
            List<Feature.SubFeature> subFeatures = feature.getSubFeatures();

            // Remove children for level 0 SO term types.
            if (Storage.SEQUENCE_ONTOLOGY_HIERARCHY.get(feature.type) == 0 && !subFeatures.isEmpty()) {
                feature.clearSubFeatures();
                Logging.logWarningOnce("REMOVE_SO0_CHILDREN",
                        "Features of type(s) %s are not supported to have children and associated children of %s will be removed."
                                .formatted(String.join(", ", getSOTerms(0)), feature.name));
                continue;
            }

            // Ensure only one level 1 SO term exists.
            // Note: This means that a feature of type "gene" cannot have a child of type "gene".
            if (countSOTerms(feature, 1) > 1) {
                Logging.logWarning("Only one of %s is allowed as the type of the feature or its children; %s is removed."
                        .formatted(String.join(", ", getSOTerms(1)), feature.name));
                storage.removeFeature(feature._id);
                continue;
            }

            // Ensure only one level 2 SO term exists.
            // Note: This means that a feature of type "mRNA" cannot have a child of type "mRNA" and features can not have multiple
            //  transcripts as children.
            if (countSOTerms(feature, 2) > 1) {
                Logging.logWarning("Only one of %s is allowed as the type of the feature or its children; %s is removed."
                        .formatted(String.join(", ", getSOTerms(2)), feature.name));
                storage.removeFeature(feature._id);
                continue;
            }

            // Adjust feature type and location for lower-level SO terms.
            // Note: To ensure a consistent hierarchy, features with a type of level > 1 - i.e., transcripts, CDS, and exon - are adjusted
            //  to type "gene".
            if (Storage.SEQUENCE_ONTOLOGY_HIERARCHY.get(feature.type) > 1) {
                feature = storage.replaceFeature(feature, adjustFeatureToGene(feature));
            }

            // Impute missing children.
            adjustSubFeature(feature, "CDS", "mRNA");
            adjustSubFeature(feature, "mRNA", "CDS");
            for (String soTerm : getSOTerms(2)) {
                if (!"mRNA".equals(soTerm)) {
                    adjustSubFeature(feature, soTerm, "exon");
                }
            }
        }
    }

    /**
     * Counts the number of features and sub-features of a specific Sequence Ontology (SO) hierarchy level.
     * <p>
     * This method calculates the total number of features and sub-features within a given {@link Feature} object that belong to the
     * specified SO hierarchy level. The hierarchy level is determined using the {@link Storage#SEQUENCE_ONTOLOGY_HIERARCHY} map.
     * <p>
     * The method includes both the main feature and its sub-features in the count if their types match the specified level.
     *
     * @param feature The {@link Feature} object whose features and sub-features are to be counted.
     * @param level   The SO hierarchy level to count (e.g., 0 for "region", 1 for "gene").
     * @return The total count of features and sub-features at the specified SO hierarchy level.
     */
    private int countSOTerms(Feature feature, int level) {
        // Retrieve the list of sub-features associated with the feature.
        List<Feature.SubFeature> subFeatures = feature.getSubFeatures();

        // Count the sub-features that match the specified SO hierarchy level and add 1 if the main feature matches the level.
        return (int) subFeatures.stream()
                .filter(sf -> Storage.SEQUENCE_ONTOLOGY_HIERARCHY.get(sf.type()) == level)
                .count() + (Storage.SEQUENCE_ONTOLOGY_HIERARCHY.get(feature.type) == level ? 1 : 0);
    }

    /**
     * Retrieve a collection of sequence ontology terms for a given level.
     *
     * @param level The level to retrieve the sequence ontology terms for.
     * @return A collection of sequence ontology terms for the specified level.
     */
    private Collection<String> getSOTerms(int level) {
        return Storage.SEQUENCE_ONTOLOGY_HIERARCHY.entrySet().stream().filter(sot -> sot.getValue().equals(level)).map(Map.Entry::getKey).collect(Collectors.toSet());
    }

    /**
     * Adjusts a feature to represent a "gene" type by updating its type, location, and children.
     * <p>
     * This method modifies the given {@link Feature} object to ensure it adheres to the "gene" type structure. It recalculates the start
     * and end positions of the feature based on its children's ranges, removes any children with level 1 Sequence Ontology (SO) terms, and
     * updates the feature's type to "gene".
     * <p>
     * The method also adds the current feature's type and location as a child and creates a new {@link Feature} object with the updated
     * attributes, type, and children.
     *
     * @param feature The {@link Feature} object to adjust. This object represents a genomic feature that needs to be converted to the
     *                "gene" type.
     * @return A new {@link Feature} object with the updated type, location, and children.
     * @throws MusialException If an error occurs during the adjustment process.
     */
    private Feature adjustFeatureToGene(Feature feature) throws MusialException {
        List<Feature.SubFeature> subFeatures = feature.getSubFeatures();

        // Calculate the new start position as the minimum of the children's start positions and the feature's current start.
        int start = Math.min(subFeatures.stream().mapToInt(Feature.SubFeature::start).min().orElse(Integer.MAX_VALUE), feature.start);

        // Calculate the new end position as the maximum of the children's end positions and the feature's current end.
        int end = Math.max(subFeatures.stream().mapToInt(Feature.SubFeature::end).max().orElse(Integer.MIN_VALUE), feature.end);

        // Remove children with level 1 SO terms.
        // Note: This may occur if a feature of type "CDS" has a child of type "gene", indicating a wrong structure in the input annotation.
        feature.clearSubFeatures(1);

        // Add the current feature's type and location as a child.
        feature.addSubFeature(feature.type, feature.start, feature.end);

        // Create a new feature with the updated type, location, and attributes.
        Feature adjustedFeature = new Feature(feature.name, feature.contig, start, end, feature.strand, "gene", feature._id);
        adjustedFeature.addAttributes(feature.getAttributes());
        for (Feature.SubFeature subFeature : subFeatures) {
            adjustedFeature.addSubFeature(subFeature.type(), subFeature.start(), subFeature.end());
        }

        return adjustedFeature;
    }

    /**
     * Imputes a feature's children based on the location ranges of an existing source children type.
     * <p>
     * This method ensures that if a source type (e.g., "CDS") exists among the sub-features of the given feature and the target type (e.g.,
     * "mRNA") does not, a new child of the target type is created. The new child is derived from the minimum start and maximum end
     * positions of the source type's ranges.
     * <p>
     * This procedure is useful for correcting missing children in genomic features, such as imputing an "mRNA" child for a "gene" when only
     * "CDS" children are present.
     *
     * @param feature    The {@link Feature} object whose children are being processed.
     * @param sourceType The type of the source sub-feature (e.g., "CDS") to derive ranges from.
     * @param targetType The type of the target sub-feature (e.g., "mRNA") to impute.
     * @throws MusialException If an error occurs during the adjustment process.
     */
    private void adjustSubFeature(Feature feature, String sourceType, String targetType) throws MusialException {
        // Retrieve the list of sub-features associated with the feature.
        List<Feature.SubFeature> subFeatures = feature.getSubFeatures();

        // Check if the source type exists and the target type does not.
        boolean hasSource = subFeatures.stream().anyMatch(sf -> sf.type().equals(sourceType));
        boolean hasTarget = subFeatures.stream().anyMatch(sf -> sf.type().equals(targetType));

        // Proceed only if the source type exists and the target type is missing.
        if (hasSource && !hasTarget) {
            // Group sub-features by their type and collect their position ranges.
            SortedMap<String, List<Tuple<Integer, Integer>>> children = subFeatures.stream()
                    .collect(Collectors.groupingBy(Feature.SubFeature::type,
                            TreeMap::new,
                            Collectors.mapping(sf -> new Tuple<>(sf.start(), sf.end()), Collectors.toList())));

            // Retrieve the position ranges for the source type.
            List<Tuple<Integer, Integer>> sourceRanges = children.get(sourceType);

            // Calculate the minimum start and maximum end positions from the source ranges.
            int start = sourceRanges.stream().mapToInt(t -> t.a).min().orElse(0);
            int end = sourceRanges.stream().mapToInt(t -> t.b).max().orElse(0);

            // If valid start and end positions are found, create a new sub-feature of the target type.
            if (start > 0 && end > 0) {
                feature.addSubFeature(targetType, start, end);
            }
        }
    }

}
