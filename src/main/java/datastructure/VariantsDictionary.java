package datastructure;

import com.aayushatharva.brotli4j.Brotli4jLoader;
import com.aayushatharva.brotli4j.encoder.Encoder;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import components.Bio;
import exceptions.MusialBioException;
import exceptions.MusialIOException;
import main.Musial;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentSkipListMap;
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
     * Software name for software version meta-information generation.
     */
    @SuppressWarnings("unused")
    public final String software = Musial.NAME + Musial.VERSION;
    /**
     * Date tag for time stamp generation.
     */
    @SuppressWarnings("unused")
    public final String date = new SimpleDateFormat("dd/MM/yyyy").format(new Date());
    /**
     * The chromosome name on which all maintained {@link VariantsDictionary#features} are located.
     */
    public final String chromosome;
    /**
     * Hierarchical map structure to store variants wrt. maintained features and samples. The first layer represents the
     * position on the chromosome. The second layer represents the variant content.
     */
    public final ConcurrentSkipListMap<Integer, ConcurrentSkipListMap<String, NucleotideVariantAnnotationEntry>>
            variants =
            new ConcurrentSkipListMap<>(Integer::compare);
    /**
     * Static property to use as the sample name/id for the wild type.
     * TODO: Move to separate class with other attribute names.
     */
    public final static String WILD_TYPE_SAMPLE_ID = "WildType";
    /**
     * Static property to use as suffix for the 'VSWAB' annotation key in {@link FeatureEntry}s and {@link SampleEntry}s.
     */
    public final static String ATTRIBUTE_VARIANT_SWAB_NAME = "_VSWAB";
    /**
     * Static property to use as key for the 'REFCONTENT' annotation key for {@link VariantsDictionary#variants}.
     */
    public final static String ATTRIBUTE_VARIANT_REFERENCE_CONTENT = "REFCONTENT";

    /**
     * Constructor of {@link VariantsDictionary}.
     *
     * @param minCoverage  {@link Double}; Minimal read coverage to use for variant filtering.
     * @param minFrequency {@link Double}; Minimal allele frequency to use for hom. variant filtering.
     * @param minHet       {@link Double}; Minimal allele frequency to use for het. variant filtering.
     * @param maxHet       {@link Double}; Maximal allele frequency to use for het. variant filtering.
     * @param minQuality   {@link Double}; Minimal Phred scaled genotyping call quality to use for variant filtering.
     * @param chromosome   {@link String}; Name of the chromosome on which maintained {@link FeatureEntry}s are located;
     *                     Has to reflect the value in the used input .vcf, .fasta and .gff files.
     */
    public VariantsDictionary(Double minCoverage, Double minFrequency, Double minHet, Double maxHet, Double minQuality,
                              String chromosome) {
        this.parameters = new VariantsDictionaryParameters(minCoverage, minFrequency, minHet, maxHet, minQuality);
        this.chromosome = chromosome;
    }

    /**
     * Adds information about a variant to {@link VariantsDictionary#variants}.
     *
     * @param featureId         {@link String}; The internal name of the feature for which this variant was detected.
     * @param referencePosition {@link Integer}; The position on the reference chromosome (1-based indexing).
     * @param variantContent    {@link String}; The alternate content of the variant; Single letter contents reflect
     *                          SNVs, multi letter contents reflect SVs.
     * @param referenceContent  {@link String}; The reference content of the variant; The length in relation to the
     *                          variantContent parameter indicates whether a SV is a insertion or deletion.
     * @param sampleId          {@link String}; The internal name of the sample of which the variant was called.
     * @param isPrimary         {@link Boolean}; Whether the variant call is a primary call, i.e. if it is the variant
     *                          with the highest frequency.
     * @param isRejected        {@link Boolean}; Whether the variant was rejected, i.e. if it failed any filtering
     *                          criteria.
     * @param quality           {@link Double}; The Phred scaled GT quality of the call.
     * @param coverage          {@link Double}; The depth of coverage of the call.
     * @param frequency         {@link Double}; The allelic frequency of the call.
     */
    public void addVariant(String featureId, int referencePosition, String variantContent,
                           @SuppressWarnings("unused") String referenceContent, String sampleId,
                           boolean isPrimary, boolean isRejected, double quality, double coverage, double frequency) {
        // Update variant information in `this.variants`.
        if (!this.variants.containsKey(referencePosition)) {
            this.variants.put(referencePosition, new ConcurrentSkipListMap<>());
        }
        if (!this.variants.get(referencePosition).containsKey(variantContent)) {
            this.variants.get(referencePosition).put(variantContent, new NucleotideVariantAnnotationEntry());
        }
        if (!this.variants.get(referencePosition).get(variantContent).occurrence.containsKey(sampleId)) {
            this.variants.get(referencePosition).get(variantContent).occurrence.put(sampleId,
                    NucleotideVariantAnnotationEntry
                            .constructSampleSpecificAnnotation(isRejected, isPrimary, quality, frequency, coverage)
            );
            this.variants.get(referencePosition).get(variantContent).annotations.put(VariantsDictionary.ATTRIBUTE_VARIANT_REFERENCE_CONTENT, referenceContent);
            if (isPrimary && !isRejected) {
                if (!samples.get(sampleId).annotations.containsKey(featureId + VariantsDictionary.ATTRIBUTE_VARIANT_SWAB_NAME)) {
                    this.samples.get(sampleId).annotations.put(featureId + VariantsDictionary.ATTRIBUTE_VARIANT_SWAB_NAME, variantContent + "@" + referencePosition);
                } else {
                    this.samples.get(sampleId).annotations
                            .put(featureId + VariantsDictionary.ATTRIBUTE_VARIANT_SWAB_NAME, this.samples.get(sampleId).annotations.get(featureId + VariantsDictionary.ATTRIBUTE_VARIANT_SWAB_NAME)
                                    .concat("|" + variantContent + "@" + referencePosition)
                            );
                }
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
        if (this.variants.containsKey(position) && this.variants.get(position).containsKey(content)) {
            this.variants.get(position).get(content).annotations.putAll(annotations);
        }
    }

    /**
     * Dumps this {@link VariantsDictionary} instance into a JSON format file. The output is compressed with gzip
     * dependent on the value of {@link Musial#COMPRESS}.
     *
     * @param outfile {@link File} object pointing to the file the VDict should be written to.
     * @throws MusialIOException If the output generation fails.
     */
    public void dump(File outfile) throws MusialIOException {
        if (!outfile.exists()) {
            throw new MusialIOException("The specified output file does not exist:\t" + outfile.getAbsolutePath());
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
                             new FileOutputStream(dumpFilePath.endsWith(".gz") ? dumpFilePath : dumpFilePath + ".gz");
                     Writer writer = new OutputStreamWriter(new GZIPOutputStream(output), StandardCharsets.UTF_8)) {
                    writer.write(dumpString);
                }

                // Write compressed (brotli) JSON.
                gson = new GsonBuilder().create();
                dumpString = gson.toJson(this);
                Brotli4jLoader.ensureAvailability();
                byte[] compressed = Encoder.compress(dumpString.getBytes());
                Files.write( Paths.get( dumpFilePath.endsWith(".br") ? dumpFilePath : dumpFilePath + ".br" ), compressed );
            }
            // Write to pretty JSON.
            gson = new GsonBuilder().setPrettyPrinting().create();
            dumpString = gson.toJson(this);
            FileWriter variantDBWriter = new FileWriter(dumpFilePath.endsWith(".gz") ? dumpFilePath.replace(".gz", "") :
                    dumpFilePath);
            variantDBWriter.write(dumpString);
            variantDBWriter.close();
        } catch (IOException e) {
            throw new MusialIOException("Failed to write to output file:\t" + outfile.getAbsolutePath());
        }
    }

    /**
     * Extracts a {@link HashMap} of {@link Integer}/{@link String} key/value pairs reflecting all (filtered) variants
     * for one sample and feature.
     *
     * @param fId {@link String}; The name/id of the feature for which variants shall be extracted.
     * @param sId {@link String}; The name/id of the sample for which variants shall be extracted.
     * @return {@link HashMap} of variants extracted for the specified feature and sample.
     */
    public HashMap<Integer, String> getSampleVariants(String fId, String sId) {
        String sampleVSwab = samples.get(sId).annotations.get(fId + VariantsDictionary.ATTRIBUTE_VARIANT_SWAB_NAME);
        if (sampleVSwab == null) {
            return null;
        } else {
            HashMap<Integer, String> sampleVariants = new HashMap<>();
            int variantPosition;
            int relativeVariantPosition;
            String variantContent;
            String relativeVariantContent;
            for (String sV : sampleVSwab.split("\\|")) {
                variantContent = sV.split("@")[0];
                variantPosition = Integer.parseInt(sV.split("@")[1]);
                relativeVariantPosition = features.get(fId).isSense ? variantPosition - features.get(fId).start + 1 :
                        features.get(fId).end - variantPosition + 1;
                relativeVariantContent = features.get(fId).isSense ? variantContent : Bio.reverseComplement(variantContent);
                sampleVariants.put(relativeVariantPosition, relativeVariantContent);
            }
            return sampleVariants;
        }
    }

    /**
     * Extracts a {@link String} yielding the nucleotide sequence with all variants of one sample being incorporated for one feature.
     * <p>
     * All deletions are removed from the resulting sequence.
     *
     * @param fId {@link String}; The name/id of the feature for which the sequence shall be extracted.
     * @param sId {@link String}; The name/id of the sample for which the sequence shall be extracted.
     * @return {@link String}; Nucleotide sequence of one sample wrt. one feature.
     */
    public String getNucleotideSequence(String fId, String sId) {
        if (sId.equals(VariantsDictionary.WILD_TYPE_SAMPLE_ID)
                || samples.get(sId).annotations.get(fId + VariantsDictionary.ATTRIBUTE_VARIANT_SWAB_NAME) == null) {
            return features.get(fId).nucleotideSequence;
        } else {
            char[] referenceSequence = features.get(fId).isSense ? features.get(fId).nucleotideSequence.toCharArray() : Bio.reverseComplement(features.get(fId).nucleotideSequence).toCharArray();
            StringBuilder sampleSequenceBuilder = new StringBuilder();
            HashMap<Integer, String> sampleVariants = getSampleVariants(fId, sId);
            String sampleVariant;
            long skipPositions = 0;
            for (int i = 1; i <= referenceSequence.length; i++) {
                if (skipPositions > 0) {
                    skipPositions -= 1;
                    continue;
                }
                if (sampleVariants.containsKey(i)) {
                    sampleVariant = sampleVariants.get(i);
                    if (sampleVariant.contains("-")) {
                        skipPositions = sampleVariant.chars().filter(c -> c == '-').count();
                        sampleVariant = sampleVariant.replace("-", "");
                    }
                    sampleSequenceBuilder.append(sampleVariant);
                } else {
                    sampleSequenceBuilder.append(referenceSequence[i - 1]);
                }
            }
            String sampleSequence = sampleSequenceBuilder.toString();
            return features.get(fId).isSense ? sampleSequence : Bio.reverseComplement(sampleSequence);
        }
    }

    /**
     * Extracts a {@link HashMap} of {@link String}/{@link String} key/value pairs reflecting all (filtered) variants
     * for one proteoform and feature.
     * <p>
     * - Returns null if no protein is allocated to the feature.
     * - Keys of the returned map are formatted as X+Y; X is the position wrt. the reference protein and Y the number of
     * inserted positions.
     *
     * @param fId  {@link String}; The name/id of the feature for which variants shall be extracted.
     * @param pfId {@link String}; The name/id of the proteoform for which variants shall be extracted.
     * @return {@link HashMap} of variants extracted for the specified feature and proteoform.
     */
    public HashMap<String, String> getProteoformVariants(String fId, String pfId) {
        // TODO: Store proteoform annotation attribute name as constant.
        if (features.get(fId).allocatedProtein == null) {
            return null;
        }
        String proteoformVSwab = features.get(fId).allocatedProtein.proteoforms.get(pfId).annotations.get("VSWAB");
        HashMap<String, String> proteoformVariants = new HashMap<>();
        if (!proteoformVSwab.equals("")) {
            String variantPosition;
            String variantContent;
            for (String pV : proteoformVSwab.split("\\|")) {
                variantContent = pV.split("@")[0];
                variantPosition = pV.split("@")[1];
                proteoformVariants.put(variantPosition, variantContent);
            }
        }
        return proteoformVariants;
    }

    /**
     * Extracts a {@link String} yielding the amino-acid sequence with all variants of one proteoform being incorporated
     * for one feature.
     * <p>
     * The translated reference feature sequence will be used as the reference amino-acid sequence.
     *
     * @param fId  {@link String}; The name/id of the feature for which the sequence shall be extracted.
     * @param pfId {@link String}; The name/id of the proteoform for which the sequence shall be extracted.
     * @return {@link String}; Amino-acid sequence of one proteoform wrt. one feature.
     * @throws MusialBioException If the translation of the reference sequence fails.
     */
    public String getProteoformSequence(String fId, String pfId) throws MusialBioException {
        if (features.get(fId).allocatedProtein == null) {
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
        char[] referenceProteoformSequenceContent = Bio.translateNucSequence(features.get(fId).nucleotideSequence, true, true, features.get(fId).isSense).toCharArray();
        for (int i = 0; i < referenceProteoformSequenceContent.length; i++) {
            proteoformContent.put((i + 1) + "+0", String.valueOf(referenceProteoformSequenceContent[i]));
        }
        proteoformContent.putAll(getProteoformVariants(fId, pfId));
        StringBuilder proteoformSequenceBuilder = new StringBuilder();
        proteoformContent.values().forEach(proteoformSequenceBuilder::append);
        return proteoformSequenceBuilder.toString().replace("-", "");
    }
}
