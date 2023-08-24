package datastructure;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import exceptions.MusialException;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import main.MusialConstants;
import utility.Compression;
import utility.IO;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

/**
 * Main data structure to store information about variants wrt. features and samples.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.1
 */
@SuppressWarnings("unused")
public class MusialStorage {

    /**
     * Set of parameters used for filtering variant calls.
     */
    public final MusialStorageParameters parameters;
    /**
     * Container to store maintained genomic features.
     */
    private final HashMap<String, Feature> features = new HashMap<>();
    /**
     * Container to store samples, i.e., variant call sets from one biological sample.
     */
    private final HashMap<String, Sample> samples = new HashMap<>();
    /**
     * Specification of positions to exclude from the analysis per reference contig.
     */
    public final HashMap<String, TreeSet<Integer>> excludedPositions = new HashMap<>();
    /**
     * Indexed reference genome.
     */
    public transient IndexedFastaSequenceFile reference;

    /**
     * Constructor of {@link MusialStorage}.
     *
     * @param minimalCoverage              {@link Double}; Minimal read coverage to use for variant filtering.
     * @param minimalFrequency             {@link Double}; Minimal allele frequency to use for hom. variant filtering.
     * @param minimalHeterozygousFrequency {@link Double}; Minimal allele frequency to use for het. variant filtering.
     * @param maximalHeterozygousFrequency {@link Double}; Maximal allele frequency to use for het. variant filtering.
     * @param minimalQuality               {@link Double}; Minimal Phred scaled call quality to use for variant filtering.
     */
    public MusialStorage(IndexedFastaSequenceFile reference, Double minimalCoverage, Double minimalFrequency, Double minimalHeterozygousFrequency, Double maximalHeterozygousFrequency, Double minimalQuality) {
        this.parameters = new MusialStorageParameters(minimalCoverage, minimalFrequency, minimalHeterozygousFrequency, maximalHeterozygousFrequency, minimalQuality);
        this.reference = reference;
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
        return this.features.get(name);
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
        if (feature == null) {
            return null;
        } else {
            return reference.getSubsequenceAt(feature.chromosome, feature.start, feature.end).getBaseString();
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
        return this.samples.get(name);
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
     * Dumps this {@link MusialStorage} instance into a .json(.br) format file. The output is optionally compressed with brotli.
     *
     * @param outfile  {@link File} object pointing to the file this {@link MusialStorage} should be written to.
     * @param compress {@link Boolean} indicating whether the output file should be compressed.
     * @throws MusialException If the output generation fails.
     */
    public void dump(File outfile, boolean compress) throws MusialException {
        if (!outfile.exists()) {
            throw new MusialException("(I/O) The specified output file does not exist:\t" + outfile.getAbsolutePath());
        }
        try {
            String dumpString;
            String dumpFilePath = outfile.getAbsolutePath();
            if (compress) {
                // Write to compressed (brotli) JSON.
                dumpString = new Gson().toJson(this);
                Files.write(
                        Paths.get(dumpFilePath.endsWith(".br") ? dumpFilePath : dumpFilePath + ".br"),
                        Compression.brotliEncodeStringToBytes(dumpString));
            } else {
                // Write to pretty JSON.
                dumpString = new GsonBuilder().setPrettyPrinting().create().toJson(this);
                Writer dumpWriter = new FileWriter(dumpFilePath.endsWith(".json") ? dumpFilePath : dumpFilePath + ".json");
                dumpWriter.write(dumpString);
                dumpWriter.close();
            }
        } catch (IOException e) {
            throw new MusialException("(I/O) Failed to write to output file:\t" + outfile.getAbsolutePath());
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
                for (Map.Entry<String, VariantAnnotation> variantEntry : feature.getNucleotideVariants(variantPosition).entrySet()) {
                    String variantContent = variantEntry.getKey();
                    VariantAnnotation variantAnnotation = variantEntry.getValue();
                    vcfContent
                            .append(feature.chromosome).append("\t")
                            .append(variantPosition).append("\t")
                            .append(".\t")
                            .append(variantAnnotation.getProperty(MusialConstants.REFERENCE_CONTENT)).append("\t")
                            .append(variantContent.replace("-", "")).append("\t")
                            .append("1000\t")
                            .append(".\t")
                            .append("\t").append(IO.LINE_SEPARATOR);

                }
            }
        }
        return vcfContent.toString();
    }

}
