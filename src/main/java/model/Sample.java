package model;

import com.google.gson.Gson;
import com.google.gson.TypeAdapter;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonWriter;
import htsjdk.samtools.util.Tuple;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.MutableTriple;
import util.Bio;
import util.Constants;
import util.Logging;

import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents a single biological sample.
 * <p>
 * This class provides a structure to store and manage variant calls associated with a biological sample. Furthermore, it allows the
 * association with an allele per feature. It extends the {@link Attributes} class to inherit functionality for managing attributes
 * associated with the sample.
 * <p>
 * Samples are stored in the {@link Storage#samples} property of the model.
 */
public class Sample extends Attributes {

    /**
     * Unique identifier of this sample.
     * <p>
     * This field serves as the unique identifier for the sample and is used to reference it in the model. This should be at best a database
     * identifier, such as a NCBI accession number, or strain identifier.
     */
    public final String _id;

    /**
     * A map that assigns features to their corresponding alleles.
     * <p>
     * This {@link Map} stores the relationship between feature names and their associated allele identifiers. The keys represent the names
     * of the features, and the values represent the unique identifiers of the alleles. This structure is used to track which allele is
     * associated with each feature in the sample.
     */
    private final Map<String, String> alleles;

    /**
     * A map that organizes variant calls by contig identifiers and positions.
     * <p>
     * This {@link Map} stores contig identifiers as keys, each associated with a {@link Map} that maps positions (as {@link Integer}) to
     * {@link VariantCall} objects. This structure allows efficient storage and retrieval of variant calls for specific contigs and
     * positions.
     */
    private final Map<String, Map<Integer, VariantCall>> variantCalls;

    /**
     * A transient cache for storing novel variant calls.
     * <p>
     * This {@link HashMap} is used to temporarily store novel variant calls during runtime. The keys represent contig identifiers, and the
     * values are {@link HashSet} objects containing the positions of the novel variant calls within the contig.
     * <p>
     * The cache is initialized with a default capacity of 2 to optimize memory allocation for typical use cases. Being marked as
     * {@code transient}, this field is excluded from serialization, as it is only relevant during the execution of the program.
     */
    transient HashMap<String, HashSet<Integer>> novelCalls = new HashMap<>(2);

    /**
     * Stores information about an upstream deletion affecting this sample.
     * <p>
     * This transient field holds an {@link UpstreamDeletion} object that represents details about an upstream deletion, including the
     * contig identifier, start and end positions, and whether the deletion is filtered. The field is marked as {@code transient} to exclude
     * it from serialization.
     */
    private transient UpstreamDeletion upstreamDeletion = null;

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
     * @param filtered         A boolean flag indicating whether the deletion is filtered.
     */
    private record UpstreamDeletion(String contigIdentifier, int start, int end, boolean filtered) {
    }

    /**
     * Constructs a new {@link Sample} instance with the specified id and initial capacity for the alleles map.
     * <p>
     * This constructor initializes a {@link Sample} object with the given id and allocates a {@link HashMap} for the {@link #alleles} field
     * with the specified initial capacity. The {@link #_id} field is set to the provided id, and the superclass constructor is invoked to
     * initialize inherited properties.
     *
     * @param identifier The id of the sample, used as its unique identifier.
     * @param capacity   The expected initial capacity of the {@link #alleles} map.
     */
    Sample(String identifier, int capacity) {
        super();
        this._id = identifier;
        this.variantCalls = new HashMap<>(2);
        this.alleles = new HashMap<>(capacity);
    }

    /**
     * Associates a specific allele of a feature with this sample.
     * <p>
     * This method updates the {@link #alleles} map by setting the sequence type (allele) for the specified feature. The feature and allele
     * are identified by their id.
     *
     * @param featureIdentifier The id of the feature ({@link Feature#_id}) to associate with the allele.
     * @param alleleIdentifier  The unique identifier of the allele ({@link SequenceType#_id}) to set for the feature.
     */
    public void addRelation(String featureIdentifier, String alleleIdentifier) {
        this.alleles.put(featureIdentifier, alleleIdentifier);
    }

    /**
     * Retrieves the allele associated with a specific feature in this sample.
     * <p>
     * This method looks up the {@link #alleles} map to find the allele associated with the given feature identifier. If no association
     * exists, it returns the default reference allele defined in {@link Constants#REFERENCE}.
     *
     * @param featureIdentifier The unique identifier of the feature ({@link Feature#name}) to retrieve the associated allele for.
     * @return The unique identifier of the allele ({@link SequenceType#_id}) associated with the feature, or the default reference allele
     * if no association exists.
     */
    public String getRelatedAllele(String featureIdentifier) {
        return this.alleles.getOrDefault(featureIdentifier, Constants.REFERENCE);
    }

    /**
     * Retrieves all feature-allele associations in this sample.
     * <p>
     * This method converts the {@link #alleles} map, which stores feature names as keys and their associated allele identifiers as values,
     * into a collection of {@link Tuple} objects. Each tuple contains a feature name and its corresponding allele identifier.
     *
     * @return A {@link Set} of {@link Tuple} objects, where each tuple represents a feature-allele association.
     */
    public Set<Tuple<String, String>> getRelatedAlleles() {
        return this.alleles.entrySet().stream()
                .map(entry -> new Tuple<>(entry.getKey(), entry.getValue()))
                .collect(Collectors.toUnmodifiableSet());
    }

    /**
     * Retrieves the number of alternate alleles in this sample.
     * <p>
     * This method returns the size of the {@link #alleles} map, which represents the number of unique alleles associated with features in
     * this sample. This corresponds to the number of non-reference alleles present in the sample.
     *
     * @return The number of alleles in this sample.
     */
    public int getRelatedAllelesCount() {
        return this.alleles.size();
    }

    /**
     * Checks if there are any variant calls associated with this sample.
     * <p>
     * This method verifies whether the {@link #variantCalls} map contains any entries, indicating the presence of variant calls.
     *
     * @return {@code true} if there are variant calls, {@code false} otherwise.
     */
    public boolean hasVariantCalls() {
        return !this.variantCalls.isEmpty();
    }

    /**
     * Checks if there are any variant calls for a specific contig.
     * <p>
     * This method checks whether the {@link #variantCalls} map contains the specified contig identifier as a key and if the associated map
     * of positions is not empty.
     *
     * @param contigIdentifier The unique identifier of the contig to check for variant calls.
     * @return {@code true} if there are variant calls for the specified contig, {@code false} otherwise.
     */
    public boolean hasVariantCalls(String contigIdentifier) {
        return this.variantCalls.containsKey(contigIdentifier) && !this.variantCalls.get(contigIdentifier).isEmpty();
    }

    /**
     * Adds a variant call to the sample for a specific contig and position.
     * <p>
     * This method processes a list of alternative alleles and calculates various properties of the variant call, such as total depth,
     * normalized entropy, and the most likely allele. It also handles cases such as low frequency, low coverage, and upstream deletions.
     * The variant call is then stored in the {@link #variantCalls} map.
     *
     * @param contigIdentifier The unique identifier of the contig where the variant call occurs.
     * @param position         The position within the contig where the variant call occurs.
     * @param alternatives     A list of {@link VariantCall.CallAlternative} objects representing the alternative alleles for the variant
     *                         call.
     * @param parameters       The {@link Storage.Parameters} object containing thresholds for frequency and coverage.
     * @param origin           The origin of the variant call (e.g., the file path from which the call was derived).
     * @return The {@link VariantCall.Flag} indicating the status of the variant call (e.g., PASS, LOW_FREQUENCY, LOW_COVERAGE).
     */
    public VariantCall.Flag addVariantCall(String contigIdentifier, int position, List<VariantCall.CallAlternative> alternatives,
                                           Storage.Parameters parameters, Path origin) {

        // Retrieve the current upstream deletion affecting this sample.
        UpstreamDeletion upstreamDeletion = this.upstreamDeletion;

        // Check if a variant call already exists for the given contig and position.
        if (variantCalls.containsKey(contigIdentifier) && variantCalls.get(contigIdentifier).containsKey(position)) {
            List<VariantCall.CallAlternative> _alternatives = variantCalls.get(contigIdentifier).get(position).alternatives();
            int i;
            for (VariantCall.CallAlternative alternative : alternatives) {
                i = _alternatives.indexOf(alternative);
                if (i >= 0) {
                    // Update the allelic depth for an existing alternative allele.
                    VariantCall.CallAlternative existing = _alternatives.get(i);
                    _alternatives.set(i, new VariantCall.CallAlternative(existing.reference(), existing.alternative(),
                            existing.allelicDepth() + alternative.allelicDepth()));
                } else {
                    // Add a new alternative allele to the list.
                    _alternatives.add(alternative);
                }
            }
            alternatives = _alternatives;
        }

        // Sort alleles in descending order by their allelic depth (AD).
        alternatives.sort((a, b) -> Integer.compare(b.allelicDepth(), a.allelicDepth()));

        // Calculate the total observed depth of coverage.
        int totalDepth = alternatives.stream().mapToInt(VariantCall.CallAlternative::allelicDepth).sum();

        // Calculate normalized entropy for the call context.
        double callEntropy = alternatives.size() == 1 ? 0.0 : -1 * (alternatives.stream().mapToDouble(alternative -> {
            float frequency = alternative.allelicDepth() / (float) totalDepth;
            return frequency == 0 ? 0 : frequency * (Math.log(frequency) / Constants.LOG2);
        }).sum()) / (Math.log(alternatives.size()) / Constants.LOG2);

        // Access the allele with the highest depth of coverage.
        VariantCall.CallAlternative allele = alternatives.get(0);
        VariantCall.Flag flag = allele.alternative().equals(Constants.DOT) ? VariantCall.Flag.REFERENCE_CALL : VariantCall.Flag.PASS;

        // Compute the actual frequency of the selected allele.
        float frequency = allele.allelicDepth() / (float) totalDepth;

        // Set call prefix for low frequency or coverage.
        if (frequency < parameters.minimalFrequency()) flag = VariantCall.Flag.LOW_FREQUENCY;
        if (totalDepth < parameters.minimalCoverage()) flag = VariantCall.Flag.LOW_COVERAGE;
        boolean isFiltered = (flag.equals(VariantCall.Flag.LOW_FREQUENCY) || flag.equals(VariantCall.Flag.LOW_COVERAGE));

        // Handle missing allele due to an upstream deletion.
        if (!isFiltered && allele.alternative().equals("*")) {
            if (Objects.isNull(upstreamDeletion)
                    || (upstreamDeletion.contigIdentifier.equals(contigIdentifier) && upstreamDeletion.start <= position && position <= upstreamDeletion.end && upstreamDeletion.filtered)
                    || (upstreamDeletion.contigIdentifier.equals(contigIdentifier) && position > upstreamDeletion.end)) {
                flag = VariantCall.Flag.MISSING_UPSTREAM_DELETION;
                isFiltered = true;
                Logging.logWarningOnce("UNEXPLAINED_DELETION",
                        String.format("Possible error in genotype data. Called deleted allele (*) is not explained by an " +
                                        "upstream deletion at site %s %d for sample %s in file %s.",
                                contigIdentifier, position, _id, origin));
            }
        }

        // Set deleted downstream positions if the current accepted call is a deletion.
        if (Bio.isDeletion(allele.reference(), allele.alternative(), true)) {
            this.upstreamDeletion = new UpstreamDeletion(
                    contigIdentifier, position + StringUtils.indexOf(allele.alternative(), Constants.GAP_CHAR),
                    position + StringUtils.lastIndexOf(allele.alternative(), Constants.GAP_CHAR), isFiltered);
        }

        // Store the variant call in the calls map if it is not a reference call.
        if (!flag.equals(VariantCall.Flag.REFERENCE_CALL)) {
            variantCalls.computeIfAbsent(contigIdentifier, k -> new HashMap<>(128));
            variantCalls.get(contigIdentifier).put(position, new VariantCall(flag, totalDepth, callEntropy, alternatives));
            // Todo: There may be a more efficient way to track novel variant calls.
            novelCalls.computeIfAbsent(contigIdentifier, k -> new HashSet<>(128));
            novelCalls.get(contigIdentifier).add(position);
        }

        return flag;
    }

    /**
     * Retrieves all variant calls in this sample.
     * <p>
     * This method processes the {@link #variantCalls} map to extract all variant calls across all contigs and positions. It converts the
     * entries into a list of {@link MutableTriple} objects, where each triple contains:
     * <ul>
     *   <li>The contig identifier as a {@link String}.</li>
     *   <li>The position of the variant as an {@link Integer}.</li>
     *   <li>The {@link VariantCall} object representing the variant call.</li>
     * </ul>
     *
     * @return A {@link List} of {@link ImmutableTriple} objects representing all variant calls in this sample.
     */
    public List<ImmutableTriple<String, Integer, VariantCall>> getVariantCalls() {
        return this.variantCalls.entrySet().stream()
                .flatMap(entry -> entry.getValue().entrySet().stream()
                        .map(callEntry -> new ImmutableTriple<>(entry.getKey(), callEntry.getKey(), callEntry.getValue())))
                .toList();
    }

    /**
     * Retrieves a list of variant calls for a specific contig.
     * <p>
     * This method retrieves variant calls associated with the specified contig identifier. It can return either novel variant calls or all
     * variant calls based on the value of the {@code novel} parameter.
     * <p>
     * If {@code novel} is {@code true}, the method retrieves only the novel variant calls for the specified contig. These are fetched from
     * the {@link #novelCalls} map. If {@code novel} is {@code false}, the method retrieves all variant calls for the contig from the
     * {@link #variantCalls} map.
     * <p>
     * The result is a list of {@link Tuple} objects, where each tuple contains:
     * <ul>
     *   <li>The position of the variant as an {@link Integer}.</li>
     *   <li>The {@link VariantCall} object representing the variant call.</li>
     * </ul>
     * If the contig identifier does not exist in the respective map, an empty list is returned.
     *
     * @param contigIdentifier The unique identifier of the contig to retrieve variant calls for.
     * @param novel            A boolean flag indicating whether to retrieve only novel variant calls ({@code true}) or all variant calls
     *                         ({@code false}).
     * @return A {@link List} of {@link Tuple} objects representing the variant calls for the specified contig.
     */
    public List<Tuple<Integer, VariantCall>> getVariantCalls(String contigIdentifier, boolean novel) {
        if (novel) {
            return this.novelCalls.getOrDefault(contigIdentifier, new HashSet<>()).stream()
                    .map(position -> new Tuple<>(position, this.variantCalls.get(contigIdentifier).get(position))).toList();
        } else {
            return this.variantCalls.getOrDefault(contigIdentifier, new HashMap<>()).entrySet().stream()
                    .map(callEntry -> new Tuple<>(callEntry.getKey(), callEntry.getValue())).toList();
        }
    }

    /**
     * Converts the sample to its string representation.
     * <p>
     * This method returns the unique identifier of the sample as its string representation.
     *
     * @return A {@link String} representing the unique identifier of the sample.
     */
    public String toString() {
        return this._id;
    }

    /**
     * Computes the hash code for this sample.
     * <p>
     * This method calculates the hash code of the sample based on its unique identifier.
     *
     * @return The hash code of the sample.
     */
    public int hashCode() {
        return this._id.hashCode();
    }

    /**
     * Compares this sample to another object for equality.
     * <p>
     * This method checks if the provided object is the same instance as this sample. If not, it verifies that the object is of the same
     * class and compares their unique identifiers for equality.
     *
     * @param obj The object to compare with this {@link Sample} instance.
     * @return {@code true} if the objects are the same instance or if their unique identifiers are equal; {@code false} otherwise.
     */
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        Sample that = (Sample) obj;
        return this._id.equals(that._id);
    }

    /**
     * Creates a custom {@link TypeAdapter} for the {@link Sample} class.
     * <p>
     * This method defines a custom {@link TypeAdapter} to handle the serialization and deserialization of {@link Sample} objects. The
     * adapter uses Gson's default adapter for most operations but adds custom behavior during deserialization to initialize the transient
     * {@link #novelCalls} field.
     *
     * @return A {@link TypeAdapter} for the {@link Sample} class.
     */
    public static TypeAdapter<Sample> typeAdapter() {

        return new TypeAdapter<>() {

            // Default adapter for Sample objects provided by Gson
            final TypeAdapter<Sample> defaultAdapter = new Gson().getAdapter(Sample.class);

            /**
             * Serializes a {@link Sample} object to JSON.
             * <p>
             * This method delegates the serialization process to the default adapter.
             *
             * @param out   The {@link JsonWriter} to write the JSON output.
             * @param value The {@link Sample} object to serialize.
             * @throws IOException If an I/O error occurs during writing.
             */
            @Override
            public void write(JsonWriter out, Sample value) throws IOException {
                defaultAdapter.write(out, value);
            }

            /**
             * Deserializes a {@link Sample} object from JSON.
             * <p>
             * This method delegates the deserialization process to the default adapter and then initializes
             * the transient {@code cache} field to ensure the {@link Sample} object is fully functional.
             *
             * @param in The {@link JsonReader} to read the JSON input.
             * @return The deserialized {@link Sample} object.
             * @throws IOException If an I/O error occurs during reading.
             */
            @Override
            public Sample read(JsonReader in) throws IOException {
                Sample sample = defaultAdapter.read(in); // Deserialize using the default adapter
                sample.novelCalls = new HashMap<>();
                return sample;
            }
        };
    }

}
