package datastructure;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import exceptions.MusialException;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import main.Constants;
import utility.Compression;
import utility.IO;
import utility.Logger;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.nio.file.Files;
import java.util.*;

/**
 * Main data structure to store information about variants wrt. features and samples.
 *
 * @author Simon Hackl
 */
@SuppressWarnings("unused")
public class MusialStorage extends InfoContainer {

    /**
     * Set of parameters used for filtering variant calls.
     */
    public final MusialStorageParameters parameters;
    /**
     * Container to store maintained genomic features.
     */
    private final HashMap<String, Feature> features;
    /**
     * Container to store samples, i.e., variant call sets from one biological sample.
     */
    private final HashMap<String, Sample> samples;
    /**
     * Specification of positions to exclude from the analysis per reference contig.
     */
    private final HashMap<String, TreeSet<Integer>> excludedPositions;
    /**
     * Specification of positions and explicit variants to exclude from the analysis per reference contig.
     */
    public final HashMap<String, HashMap<Integer, HashSet<String>>> excludedVariants;
    /**
     * Stores already seen variants to avoid redundant form names. No permanent storage in MUSIAL session.
     */
    private final transient HashMap<Integer, String> formNameCodes = new HashMap<>();
    /**
     * Indexed reference sequences. No permanent storage in MUSIAL session.
     */
    private transient IndexedFastaSequenceFile reference;

    /**
     * Constructor of {@link MusialStorage}.
     *
     * @param minimalCoverage  {@link Double}; Minimal read coverage to use for allele filtering.
     * @param minimalFrequency {@link Double}; Minimal allele frequency to use for allele filtering.
     */
    public MusialStorage(IndexedFastaSequenceFile reference, Double minimalCoverage, Double minimalFrequency) {
        super();
        this.parameters = new MusialStorageParameters(minimalCoverage, minimalFrequency);
        this.reference = reference;
        addInfo("reference_length", String.valueOf(reference.nextSequence().length()));
        this.features = new HashMap<>();
        this.samples = new HashMap<>();
        this.excludedPositions = new HashMap<>();
        this.excludedVariants = new HashMap<>();
    }

    /**
     * Constructor of {@link MusialStorage}.
     *
     * @param minimalCoverage  {@link Double}; Minimal read coverage to use for allele filtering.
     * @param minimalFrequency {@link Double}; Minimal allele frequency to use for allele filtering.
     */
    public MusialStorage(IndexedFastaSequenceFile reference, Double minimalCoverage, Double minimalFrequency, Collection<Feature> features, Collection<Sample> samples, HashMap<String, TreeSet<Integer>> excludedPositions, HashMap<String, HashMap<Integer, HashSet<String>>> excludedVariants) throws MusialException {
        this.parameters = new MusialStorageParameters(minimalCoverage, minimalFrequency);
        this.reference = reference;
        this.features = new HashMap<>(features.size());
        for (Feature feature : features) {
            addFeature(feature);
        }
        this.samples = new HashMap<>(samples.size());
        for (Sample sample : samples) {
            sample.imputeVcfFileReader();
            addSample(sample);
        }
        this.excludedPositions = excludedPositions;
        this.excludedVariants = excludedVariants;
        addInfo("reference_length", String.valueOf(reference.nextSequence().length()));
    }

    /**
     * Sets the reference sequence associated with this instance to the passed indexed reference feature file.
     *
     * @param reference Instance of {@link IndexedFastaSequenceFile}.
     */
    public void setReference(IndexedFastaSequenceFile reference) {
        this.reference = reference;
    }

    /**
     * Checks whether {@code position} is excluded on {@code contig} of {@link #reference}.
     *
     * @param contig   Contig (name) to check fo exclusion.
     * @param position Position to check for exclusion.
     * @return True if {@code position} on feature {@code featureName} is excluded from the analysis.
     */
    public boolean isPositionExcluded(String contig, int position) {
        return excludedPositions.containsKey(contig) && excludedPositions.get(contig).contains(position);
    }

    /**
     * Checks whether {@code variant} at {@code position} is excluded on {@code contig} of {@link #reference}.
     *
     * @param contig   Contig (name) to check fo exclusion.
     * @param position Position to check for exclusion.
     * @param variant  Variant to check for exclusion.
     * @return True if {@code variant} at {@code position} on feature {@code featureName} is excluded from the analysis.
     */
    public boolean isVariantExcluded(String contig, int position, String variant) {
        return excludedVariants.containsKey(contig) && excludedVariants.get(contig).containsKey(position) && excludedVariants.get(contig).get(position).contains(variant);
    }

    /**
     * Add an entry to this {@link #features}. Existing entries will be overwritten. The {@link Feature#name} value will be used as key.
     *
     * @param feature Entry to add.
     */
    public void addFeature(Feature feature) {
        this.features.put(feature.name, feature);
    }

    /**
     * Getter method for {@link #features}.
     *
     * @param name Name of the feature to query.
     * @return The queried feature or null, if no feature is stored at key name.
     */
    public Feature getFeature(String name) {
        return this.features.getOrDefault(name, null);
    }

    /**
     * Query whether a feature is stored in this instance by its name.
     *
     * @param name Name of the feature to query.
     * @return Whether the queried feature is present in this instance.
     */
    @SuppressWarnings("BooleanMethodIsAlwaysInverted")
    public boolean hasFeature(String name) {
        return this.features.containsKey(name);
    }

    /**
     * @return String iterator over the key set of {@link #features}.
     */
    public Iterator<String> getFeatureNameIterator() {
        return this.features.keySet().iterator();
    }

    /**
     * @param position The position to restrict features to.
     * @return String iterator over the key set of {@link #features} covering the specified position.
     */
    public Iterator<String> getFeatureNameIterator(int position) {
        return this.features.keySet().stream().filter(featureName -> this.getFeature(featureName).start <= position && this.getFeature(featureName).end >= position).iterator();
    }

    /**
     * Constructs the base sequence of the specified feature wrt. {@link #reference}.
     *
     * @param featureName The name of the feature to query.
     * @return The base sequence of the specified feature or null if the feature does not exist.
     */
    public String getReferenceSequenceOfFeature(String featureName) {
        Feature feature = getFeature(featureName);
        if (feature == null || reference == null) {
            return null;
        } else {
            return reference.getSubsequenceAt(feature.contig, feature.start, feature.end).getBaseString();
        }
    }

    /**
     * Removes a feature from the currently stored features.
     *
     * @param featureName The name of the feature to remove.
     */
    public void removeFeature(String featureName) {
        this.features.remove(featureName);
    }

    /**
     * Add an entry to this {@link #samples}. Existing entries will be overwritten. The {@link Sample#name} value will be used as key.
     *
     * @param sample Entry to add.
     */
    public void addSample(Sample sample) {
        this.samples.put(sample.name, sample);
    }

    /**
     * Getter method for {@link #samples}.
     *
     * @param name Name of the sample to query.
     * @return The queried sample or null, if no feature is stored at key name.
     */
    public Sample getSample(String name) {
        return this.samples.getOrDefault(name, null);
    }

    /**
     * Query whether a sample is stored in this instance by its name.
     *
     * @param name Name of the sample to query.
     * @return Whether the queried sample is present in this instance.
     */
    @SuppressWarnings("BooleanMethodIsAlwaysInverted")
    public boolean hasSample(String name) {
        return this.samples.containsKey(name);
    }

    /**
     * @return String iterator over the key set of {@link #samples}.
     */
    public Iterator<String> getSampleNameIterator() {
        return this.samples.keySet().iterator();
    }


    /**
     * Getter method for {@link #samples} size.
     *
     * @return The number of entries stored in {@link #samples}.
     */
    public int getNumberOfSamples() {
        return this.samples.size();
    }

    /**
     * Dumps this {@link MusialStorage} instance into a .json(.gz) format file. The output is optionally compressed with gzip.
     *
     * @param path     {@link String} specifying file path to which {@link MusialStorage} should be written to.
     * @param compress {@link Boolean} indicating whether the output file should be compressed.
     * @throws MusialException If the output generation fails.
     */
    public void dump(String path, boolean compress) throws MusialException {
        if (compress && !path.endsWith(".gz"))
            path += ".gz";
        else if (!compress && !path.endsWith(".json"))
            path += ".json";
        File outFile = new File(path);
        Logger.logStatus("Dump to file `" + path + "`");
        try {
            String dump;
            if (compress) {
                // Write to compressed (gzip) JSON.
                dump = new Gson().toJson(this);
                Files.write(
                        outFile.toPath(),
                        Compression.gzip(dump));
            } else {
                // Write to pretty JSON.
                GsonBuilder gsonBuilder = new GsonBuilder().setPrettyPrinting();
                dump = gsonBuilder.create().toJson(this, MusialStorage.class);
                Writer writer = new FileWriter(outFile);
                writer.write(dump);
                writer.close();
            }
        } catch (IOException e) {
            throw new MusialException("(I/O) Failed to write to output file (" + path + "):\t" + e.getMessage());
        }
    }

    /**
     * Constructs a dummy vcf file from all variants stored in this instance.
     *
     * @return The dummy vcf file content.
     */
    public String toVcf() {
        StringBuilder vcfContent = new StringBuilder();
        vcfContent.append("##fileformat=VCFv4.2").append(IO.LINE_SEPARATOR).append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").append(IO.LINE_SEPARATOR);
        for (Feature feature : this.features.values()) {
            for (Integer variantPosition : feature.getNucleotideVariantPositions()) {
                for (Map.Entry<String, VariantInformation> variantEntry : feature.getNucleotideVariantsAt(variantPosition).entrySet()) {
                    String variantContent = variantEntry.getKey();
                    if (variantContent.contains(Constants.ANY_NUCLEOTIDE_STRING))
                        continue;
                    VariantInformation variantInformation = variantEntry.getValue();
                    vcfContent
                            .append(feature.contig).append("\t")
                            .append(variantPosition).append("\t")
                            .append(".\t").append(variantInformation.referenceContent.replace("-", "")).append("\t").append(variantContent.replace("-", "")).append("\t")
                            .append("1000\t")
                            .append(".\t")
                            .append("\t").append(IO.LINE_SEPARATOR);

                }
            }
        }
        return vcfContent.toString();
    }

    /**
     * Used to collide form names based on mapping stored in {@link MusialStorage#formNameCodes}.
     * <p>
     * Specifically, this is used to map hash keys of variant fingerprints (see {@link Form#variants}) to human-readable names.
     *
     * @param hashKey Key to store.
     * @param value   Value to store.
     * @return Returns {@code value} or the value stored at {@code hashKey}, if present.
     */
    public String getFormName(int hashKey, String value) {
        if (formNameCodes.containsKey(hashKey)) {
            return formNameCodes.get(hashKey);
        } else {
            this.formNameCodes.put(hashKey, value);
            return value;
        }
    }

}
