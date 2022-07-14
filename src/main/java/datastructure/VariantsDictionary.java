package datastructure;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import exceptions.MusialBioException;
import exceptions.MusialIOException;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;

import main.Musial;
import components.Bio;

/**
 * TODO
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class VariantsDictionary {

    public final VariantsDictionaryParameters parameters;

    public final HashMap<String, FeatureEntry> features = new HashMap<>();

    public final HashMap<String, SampleEntry> samples = new HashMap<>();

    @SuppressWarnings("unused")
    public final String software = Musial.NAME + Musial.VERSION;

    @SuppressWarnings("unused")
    public final String date = new SimpleDateFormat("dd/MM/yyyy").format(new Date());

    public final String chromosome;

    public final ConcurrentSkipListMap<Integer, ConcurrentSkipListMap<String, NucleotideVariantAnnotationEntry>>
            variants =
            new ConcurrentSkipListMap<>(Integer::compare);

    public final static String WILD_TYPE_SAMPLE_ID = "WildType";

    public final static String ATTRIBUTE_VARIANT_SWAB_NAME = "_VSWAB";

    public transient ConcurrentSkipListSet<String> novelVariants = new ConcurrentSkipListSet<>((s1, s2) -> {
        int p1 = Integer.parseInt(s1.split("@")[0]);
        int p2 = Integer.parseInt(s2.split("@")[0]);
        if (p1 != p2) {
            return Integer.compare(p1, p2);
        } else {
            String c1 = s1.split("@")[2];
            String c2 = s2.split("@")[2];
            return c1.compareTo(c2);
        }
    });

    public VariantsDictionary(Double minCoverage, Double minFrequency, Double minHet, Double maxHet, Double minQuality,
                              String chromosome) {
        this.parameters = new VariantsDictionaryParameters(minCoverage, minFrequency, minHet, maxHet, minQuality);
        this.chromosome = chromosome;
    }

    public void addVariant(String featureId, int referencePosition, String variantContent,
                           String referenceContent, String sampleId,
                           boolean isPrimary, boolean isRejected, double quality, double coverage, double frequency) {
        // Update variant information in `this.variants`.
        if (!this.variants.containsKey(referencePosition)) {
            this.variants.put(referencePosition, new ConcurrentSkipListMap<>());
        }
        if (!this.variants.get(referencePosition).containsKey(variantContent)) {
            this.variants.get(referencePosition).put(variantContent, new NucleotideVariantAnnotationEntry());
            this.novelVariants.add(referencePosition + "@" + referenceContent + "@" + variantContent);
        }
        if (!this.variants.get(referencePosition).get(variantContent).occurrence.containsKey(sampleId)) {
            this.variants.get(referencePosition).get(variantContent).occurrence.put(sampleId,
                    NucleotideVariantAnnotationEntry
                            .constructSampleSpecificAnnotation(isRejected, isPrimary, quality, frequency, coverage)
            );
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

    public void addVariantAnnotation(int position, String content, Map<String, String> annotations) {
        if (this.variants.containsKey(position) && this.variants.get(position).containsKey(content)) {
            this.variants.get(position).get(content).annotations.putAll(annotations);
        }
    }

    public void dump(File outfile) throws MusialIOException {
        if (!outfile.exists()) {
            throw new MusialIOException("The specified output file does not exist:\t" + outfile.getAbsolutePath());
        }
        try {
            Gson gson;
            String dumpString;
            String dumpFilePath = outfile.getAbsolutePath();
            /* Write compressed JSON.
              gson = new GsonBuilder().create();
              dumpString = gson.toJson(this);
              try (FileOutputStream output =
                       new FileOutputStream(dumpFilePath.endsWith(".gz") ? dumpFilePath : dumpFilePath + ".gz");
                   Writer writer = new OutputStreamWriter(new GZIPOutputStream(output), StandardCharsets.UTF_8)) {
                writer.write(dumpString);
              }
            */
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

    public HashMap<Integer, String> getSampleVariants(String fId, String sId) throws MusialBioException {
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

    public String getNucleotideSequence(String fId, String sId) throws MusialBioException {
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
}
