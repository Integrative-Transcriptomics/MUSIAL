package datastructure;

import com.aayushatharva.brotli4j.Brotli4jLoader;
import com.aayushatharva.brotli4j.encoder.Encoder;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import components.Bio;
import components.Logging;
import exceptions.MusialException;
import main.Musial;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

/**
 * Main data structure to store information about variants.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class VariantsDictionary {

    /**
     * Set of parameters, a {@link VariantsDictionaryParameters} instance, used for variant filtering.
     */
    public final VariantsDictionaryParameters parameters;
    /**
     * Map of {@link String} (key) / {@link FeatureEntry} (value) pairs; All (gene) features that are maintained.
     */
    public final HashMap<String, FeatureEntry> features = new HashMap<>();
    /**
     * Map of {@link String} (key) / {@link SampleEntry} (value) pairs; All samples that are maintained.
     */
    public final HashMap<String, SampleEntry> samples = new HashMap<>();
    /**
     * Hierarchical map structure to store variants wrt. the nucleotide sequence. The first layer represents the
     * position on the chromosome. The second layer represents the variant content.
     */
    public final ConcurrentSkipListMap<Integer, ConcurrentSkipListMap<String, NucleotideVariantEntry>>
            nucleotideVariants =
            new ConcurrentSkipListMap<>(Integer::compare);
    /**
     * Software name used as meta-information.
     */
    @SuppressWarnings("unused")
    public final String software = Musial.NAME + Musial.VERSION;
    /**
     * Date tag used as meta-information.
     */
    @SuppressWarnings("unused")
    public final String date = new SimpleDateFormat("dd/MM/yyyy").format(new Date());
    /**
     * TODO
     */
    public final static String FIELD_SEPARATOR_1 = "_";
    /**
     * TODO
     */
    public final static String FIELD_SEPARATOR_2 = ";";
    /**
     * TODO
     */
    public final static String NULL_VALUE = "N/A";

    /**
     * Constructor of {@link VariantsDictionary}.
     *
     * @param minCoverage  {@link Double}; Minimal read coverage to use for variant filtering.
     * @param minFrequency {@link Double}; Minimal allele frequency to use for hom. variant filtering.
     * @param minHet       {@link Double}; Minimal allele frequency to use for het. variant filtering.
     * @param maxHet       {@link Double}; Maximal allele frequency to use for het. variant filtering.
     * @param minQuality   {@link Double}; Minimal Phred scaled genotyping call quality to use for variant filtering.
     */
    public VariantsDictionary(Double minCoverage, Double minFrequency, Double minHet, Double maxHet, Double minQuality) {
        this.parameters = new VariantsDictionaryParameters(minCoverage, minFrequency, minHet, maxHet, minQuality);
    }

    /**
     * Adds information about a nucleotide variant to {@link VariantsDictionary#nucleotideVariants}.
     *
     * @param featureId  {@link String}; The internal name of the feature for which this variant was detected.
     * @param position   {@link Integer}; The position on the reference chromosome (1-based indexing).
     * @param altContent {@link String}; The alternate content of the variant; Single letter contents reflect
     *                   SNVs, multi letter contents reflect SVs.
     * @param refContent {@link String}; The reference content of the variant; The length in relation to the
     *                   variantContent parameter indicates whether a SV is a insertion or deletion.
     * @param sampleId   {@link String}; The internal name of the sample of which the variant was called.
     * @param isPrimary  {@link Boolean}; Whether the variant call is a primary call, i.e. if it is the variant
     *                   with the highest frequency.
     */
    public void addNucleotideVariant(String featureId, int position, String altContent, String refContent, String sampleId,
                                     boolean isPrimary) {
        // Update variant information.
        if (!this.nucleotideVariants.containsKey(position)) {
            this.nucleotideVariants.put(position, new ConcurrentSkipListMap<>());
        }
        if (!this.nucleotideVariants.get(position).containsKey(altContent)) {
            this.nucleotideVariants.get(position).put(altContent, new NucleotideVariantEntry());
            this.nucleotideVariants.get(position).get(altContent).annotations.put(NucleotideVariantEntry.PROPERTY_NAME_PRIMARY, String.valueOf(isPrimary));
            this.nucleotideVariants.get(position).get(altContent).annotations.put(NucleotideVariantEntry.PROPERTY_NAME_REFERENCE_CONTENT, refContent);
            this.nucleotideVariants.get(position).get(altContent).annotations.put(NucleotideVariantEntry.PROPERTY_NAME_LOCATION, this.features.get(featureId).chromosome);
        }
        // Add information about variant to temp. annotation of sample. The information is used to infer alleles later
        // and is then discarded.
        if (isPrimary) {
            assert this.samples.containsKey(sampleId);
            assert this.features.containsKey(featureId);
            assert this.nucleotideVariants.containsKey(position);
            assert this.nucleotideVariants.get(position).containsKey(altContent);
            if (!this.features.get(featureId).novelNucleotideVariants.containsKey(sampleId)) {
                this.features.get(featureId).novelNucleotideVariants.put(sampleId, position + FIELD_SEPARATOR_1 + altContent);
            } else {
                this.features.get(featureId).novelNucleotideVariants.put(
                        sampleId,
                        this.features.get(featureId).novelNucleotideVariants.get(sampleId) + FIELD_SEPARATOR_2 + position + FIELD_SEPARATOR_1 + altContent
                );
            }
        }
    }

    /**
     * Adds annotations to an existing variant.
     *
     * @param position    {@link Integer}; The position at which the variant is stored.
     * @param content     {@link String}; The content of the variant.
     * @param annotations {@link HashMap} of {@link String} key/value pairs which reflect the annotations that shall be
     *                    added.
     */
    @SuppressWarnings("unused")
    public void addVariantAnnotation(int position, String content, Map<String, String> annotations) {
        if (this.nucleotideVariants.containsKey(position) && this.nucleotideVariants.get(position).containsKey(content)) {
            this.nucleotideVariants.get(position).get(content).annotations.putAll(annotations);
        }
    }

    /**
     * Dumps this {@link VariantsDictionary} instance into a JSON format file. The output is compressed with brotli
     * dependent on the value of {@link Musial#COMPRESS}.
     *
     * @param outfile {@link File} object pointing to the file this {@link VariantsDictionary} should be written to.
     * @throws MusialException If the output generation fails.
     */
    public void dump(File outfile) throws MusialException {
        if (!outfile.exists()) {
            throw new MusialException("(I/O) The specified output file does not exist:\t" + outfile.getAbsolutePath());
        }
        try {
            Gson gson;
            String dumpString;
            String dumpFilePath = outfile.getAbsolutePath();
            if (Musial.COMPRESS) {
                // Write compressed (gzip) JSON.
                gson = new GsonBuilder().create();
                dumpString = gson.toJson(this);
                try (FileOutputStream output =
                             new FileOutputStream(dumpFilePath);
                     Writer writer = new OutputStreamWriter(new GZIPOutputStream(output), StandardCharsets.UTF_8)) {
                    writer.write(dumpString);
                }
                // Write compressed (brotli) JSON.
                gson = new GsonBuilder().create();
                dumpString = gson.toJson(this);
                Brotli4jLoader.ensureAvailability();
                byte[] compressed = Encoder.compress(dumpString.getBytes());
                Files.write(Paths.get(dumpFilePath.endsWith(".br") ? dumpFilePath : dumpFilePath + ".br"), compressed);
            } else {
                // Write to pretty JSON.
                gson = new GsonBuilder().setPrettyPrinting().create();
                dumpString = gson.toJson(this);
                Writer variantDBWriter = new FileWriter(dumpFilePath);
                variantDBWriter.write(dumpString);
                variantDBWriter.close();
            }
        } catch (IOException e) {
            throw new MusialException("(I/O) Failed to write to output file:\t" + outfile.getAbsolutePath());
        }
    }

    /**
     * Extracts a {@link HashMap} of {@link Integer}/{@link String} key/value pairs reflecting all variants wrt. one
     * sample and feature.
     * <p>
     * If specified, variants are returned wrt. the features start and end coordinates and its strand orientation.
     *
     * @param featureId         {@link String}; The name/id of the feature for which variants shall be extracted.
     * @param sampleId          {@link String}; The name/id of the sample for which variants shall be extracted.
     * @param relativeToFeature {@link Boolean}; Whether to extract variants position and content wrt. feature (true).
     * @return {@link HashMap} of variants extracted for the specified feature and sample.
     */
    public HashMap<Integer, String> getNucleotideVariants(String featureId, String sampleId, boolean relativeToFeature) {
        String sampleAllele = this.samples.get(sampleId).annotations.get("AL" + FIELD_SEPARATOR_1 + featureId);
        String concatVariants = this.features.get(featureId).alleles.get(sampleAllele).annotations.get(AlleleEntry.PROPERTY_NAME_VARIANTS);
        HashMap<Integer, String> variants = new HashMap<>();
        if (!sampleAllele.equals(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)) {
            int variantPosition;
            int relativeVariantPosition;
            String variantContent;
            String relativeVariantContent;
            for (String variant : concatVariants.split(FIELD_SEPARATOR_2)) {
                variantPosition = Integer.parseInt(variant.split(FIELD_SEPARATOR_1)[0]);
                variantContent = variant.split(FIELD_SEPARATOR_1)[1];
                if (relativeToFeature) {
                    relativeVariantPosition = features.get(featureId).isSense ? variantPosition - features.get(featureId).start + 1 :
                            features.get(featureId).end - variantPosition + 1;
                    relativeVariantContent = features.get(featureId).isSense ? variantContent : Bio.reverseComplement(variantContent);
                    variants.put(relativeVariantPosition, relativeVariantContent);
                } else {
                    variants.put(variantPosition, variantContent);
                }
            }
        }
        return variants;
    }

    /**
     * Extracts a {@link String} yielding the nucleotide sequence with all variants of one sample being incorporated for one feature.
     * <p>
     * All deletions are removed from the resulting sequence.
     *
     * @param featureId {@link String}; The name/id of the feature for which the sequence shall be extracted.
     * @param sampleId  {@link String}; The name/id of the sample for which the sequence shall be extracted.
     * @return {@link String}; Nucleotide sequence of one sample wrt. one feature.
     */
    public String getSampleNucleotideSequence(String featureId, String sampleId) {
        if (samples.get(sampleId).annotations.get("AL" + FIELD_SEPARATOR_1 + featureId).equals(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)) {
            return features.get(featureId).nucleotideSequence;
        } else {
            char[] referenceSequence = features.get(featureId).isSense ? features.get(featureId).nucleotideSequence.toCharArray() : Bio.reverseComplement(features.get(featureId).nucleotideSequence).toCharArray();
            StringBuilder sequenceBuilder = new StringBuilder();
            HashMap<Integer, String> variants = getNucleotideVariants(featureId, sampleId, true);
            String variant;
            long skipPositions = 0;
            for (int i = 1; i <= referenceSequence.length; i++) {
                if (skipPositions > 0) {
                    skipPositions -= 1;
                    continue;
                }
                if (variants.containsKey(i)) {
                    variant = variants.get(i);
                    if (variant.contains("-")) {
                        skipPositions = variant.chars().filter(c -> c == '-').count();
                        variant = variant.replace("-", "");
                    }
                    sequenceBuilder.append(variant);
                } else {
                    sequenceBuilder.append(referenceSequence[i - 1]);
                }
            }
            String sampleSequence = sequenceBuilder.toString();
            return features.get(featureId).isSense ? sampleSequence : Bio.reverseComplement(sampleSequence);
        }
    }

    /**
     * Extracts a {@link HashMap} of {@link String}/{@link String} key/value pairs reflecting all (filtered) variants
     * for one proteoform and feature.
     * <p>
     * - Returns null if the feature is none-coding.
     * - Keys of the returned map are formatted as X+Y; X is the position wrt. the reference protein and Y the number of
     * inserted positions.
     *
     * @param featureId    {@link String}; The name/id of the feature for which variants shall be extracted.
     * @param proteoformId {@link String}; The name/id of the proteoform for which variants shall be extracted.
     * @return {@link HashMap} of variants extracted for the specified feature and proteoform.
     */
    public HashMap<String, String> getProteoformAminoacidVariants(String featureId, String proteoformId) {
        if (!features.get(featureId).isCodingSequence) {
            return null;
        }
        String concatVariants = features.get(featureId).proteoforms.get(proteoformId).annotations.get(ProteoformEntry.PROPERTY_NAME_VARIANTS);
        HashMap<String, String> variants = new HashMap<>();
        if (!concatVariants.equals("")) {
            String variantPosition;
            String variantContent;
            for (String variant : concatVariants.split(FIELD_SEPARATOR_2)) {
                variantPosition = variant.split(FIELD_SEPARATOR_1)[0];
                variantContent = variant.split(FIELD_SEPARATOR_1)[1];
                variants.put(variantPosition, variantContent);
            }
        }
        return variants;
    }

    /**
     * Extracts a {@link String} yielding the aminoacid sequence with all variants of one proteoform being incorporated
     * for one feature.
     * <p>
     * The translated reference feature sequence will be used as the reference aminoacid sequence.
     *
     * @param featureId    {@link String}; The name/id of the feature for which the sequence shall be extracted.
     * @param proteoformId {@link String}; The name/id of the proteoform for which the sequence shall be extracted.
     * @return {@link String}; Amino-acid sequence of one proteoform wrt. one feature.
     */
    public String getProteoformSequence(String featureId, String proteoformId) {
        if (!features.get(featureId).isCodingSequence) {
            return null;
        }
        //noinspection DuplicatedCode
        TreeMap<String, String> proteoformContent = new TreeMap<>((s1, s2) -> {
            int p1 = Integer.parseInt(s1.split("\\+")[0]);
            int p2 = Integer.parseInt(s2.split("\\+")[0]);
            if (p1 != p2) {
                return Integer.compare(p1, p2);
            } else {
                String i1 = s1.split("\\+")[1];
                String i2 = s2.split("\\+")[1];
                return i1.compareTo(i2);
            }
        });
        char[] referenceProteoformSequenceContent = features.get(featureId).translatedNucleotideSequence.toCharArray();
        for (int i = 0; i < referenceProteoformSequenceContent.length; i++) {
            proteoformContent.put((i + 1) + "+0", String.valueOf(referenceProteoformSequenceContent[i]));
        }
        proteoformContent.putAll(getProteoformAminoacidVariants(featureId, proteoformId));
        StringBuilder proteoformSequenceBuilder = new StringBuilder();
        proteoformContent.values().forEach(proteoformSequenceBuilder::append);
        return proteoformSequenceBuilder.toString().replace("-", "");
    }

    /**
     * Removes any sample ids from the passed list that are not associated with an actual {@link SampleEntry} in
     * {@link VariantsDictionary#samples}.
     * <p>
     * If an empty list ist provided, a list of the keys stored in {@link VariantsDictionary#samples} is returned.
     *
     * @param sampleIds {@link ArrayList} of sample ids ({@link String}s) to filter.
     * @return {@link ArrayList} of sample ids that are associated with an entry in this instance.
     * @throws MusialException If none of the passed sample ids is associated with an entry in this instance.
     */
    public ArrayList<String> removeInvalidSampleIds(ArrayList<String> sampleIds) throws MusialException {
        List<String> filteredSampleIds;
        if (sampleIds.size() == 0) {
            return new ArrayList<>(this.samples.keySet());
        } else {
            filteredSampleIds = sampleIds.stream().filter(this.samples::containsKey).collect(Collectors.toList());
            if (sampleIds.size() != filteredSampleIds.size()) {
                Logging.logWarning(
                        "Removed " + (sampleIds.size() - filteredSampleIds.size()) + " sample ids for which no entry is stored."
                );
            }
        }
        if (filteredSampleIds.size() == 0) {
            throw new MusialException("The number of specified valid samples is zero.");
        } else {
            return (ArrayList<String>) filteredSampleIds;
        }
    }

    /**
     * Removes any feature ids from the passed list that are not associated with an actual {@link FeatureEntry} in
     * {@link VariantsDictionary#features} or do not match with the specified criteria.
     * <p>
     * If an empty list is provided, a list of the keys stored in {@link VariantsDictionary#features} that pass the
     * specified criteria is returned.
     *
     * @param featureIds {@link ArrayList} of feature ids ({@link String}s) to filter.
     * @param criteria   A {@link Function} that is additionally used to filter the input list.
     * @return {@link ArrayList} of feature ids that are associated with an entry in this instance and pass the specified criteria.
     * @throws MusialException If none of the passed/stored feature ids passes the filtering.
     */
    public ArrayList<String> removeInvalidFeatureIds(ArrayList<String> featureIds, Function<String, Boolean> criteria) throws MusialException {
        List<String> filteredFeatureIds;
        if (featureIds.size() == 0) {
            featureIds = new ArrayList<>(this.features.keySet());
        }
        filteredFeatureIds = featureIds.stream().filter(
                fId -> this.features.containsKey(fId) && criteria.apply(fId)
        ).collect(Collectors.toList());
        if (featureIds.size() != filteredFeatureIds.size()) {
            Logging.logWarning(
                    "Removed " + (featureIds.size() - filteredFeatureIds.size()) + " feature ids."
            );
        }
        if (filteredFeatureIds.size() == 0) {
            throw new MusialException("The number of specified valid features is zero.");
        } else {
            return (ArrayList<String>) filteredFeatureIds;
        }
    }

    /**
     * Removes any feature ids from the passed list that are not associated with an actual {@link FeatureEntry} in
     * {@link VariantsDictionary#features}.
     * <p>
     * If an empty list is provided, a list of the keys stored in {@link VariantsDictionary#features} is returned.
     *
     * @param featureIds {@link ArrayList} of feature ids ({@link String}s) to filter.
     * @return {@link ArrayList} of feature ids that are associated with an entry in this instance.
     * @throws MusialException If none of the passed feature ids passes the filtering.
     */
    public ArrayList<String> removeInvalidFeatureIds(ArrayList<String> featureIds) throws MusialException {
        return removeInvalidFeatureIds(featureIds, s -> true);
    }


}
