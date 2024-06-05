package datastructure;

import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import main.Constants;

import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Representation of a sample, i.e., variant calls from a single biological sample.
 *
 * @author Simon Hackl
 */
@SuppressWarnings("unused")
public final class Sample extends InfoContainer {

    /**
     * Internal name of the sample.
     */
    public final String name;
    /**
     * .vcf file that is associated with this sample. No permanent storage in MUSIAL session.
     */
    public transient final File vcfFile;
    /**
     * {@link HashMap} that assigns features to alleles (names).
     */
    private final LinkedHashMap<String, String> allele;
    /**
     * {@link HashMap} that assigns features to proteoforms (names).
     */
    private final LinkedHashMap<String, String> proteoform;
    /**
     * Hierarchical map structure to store variant calls wrt. features.
     * <ul>
     *     <li>First level: A {@link Feature#name}.</li>
     *     <li>Second level: Position on feature.</li>
     *     <li>Third level: {@code CALL;DP;ALT:AD,...} with {@code CALL} as one of N (ambiguous), ? (no coverage), or the index of the alternative call. {@code DP, AD} as defined in VCF specification (<a href="https://samtools.github.io/hts-specs/VCFv4.2.pdf">samtools.github.io/hts-specs/VCFv4.2.pdf</a>}).</li>
     * </ul>
     */
    private final LinkedHashMap<String, LinkedHashMap<Integer, String>> calls;
    /**
     * {@link VCFFileReader} instance pointing to the vcf file of this sample. No permanent storage in MUSIAL session.
     */
    public transient VCFFileReader vcfFileReader;

    /**
     * Constructor of {@link Sample}.
     *
     * @param vcfFile {@link File} instances pointing to the .vcf file of this sample.
     * @param name    {@link String} representing the sample name.
     */
    public Sample(File vcfFile, String name, int noFeatures) {
        super();
        this.vcfFile = vcfFile;
        this.name = name;
        this.allele = new LinkedHashMap<>(noFeatures);
        this.proteoform = new LinkedHashMap<>(noFeatures);
        this.calls = new LinkedHashMap<>(noFeatures);
    }

    /**
     * Initialize a {@link VCFFileReader} instance pointing to the .vcf file associated with this instance.
     * <p>
     * Validate that a .tbi index of the .vcf file is present and generates one, if not.
     */
    public void imputeVcfFileReader() throws MusialException {
        try {
            // Check for the existence of a `.tbi.gz` index file of the input `.vcf` file.
            if (!new File(this.vcfFile.getAbsolutePath() + ".tbi").exists()) {
                // If none is present, an index is created and written to the same directory as the input `.vcf` file.
                TabixIndex tabixIndex = IndexFactory.createTabixIndex(
                        this.vcfFile,
                        new VCFCodec(),
                        null
                );
                tabixIndex.write(Path.of(this.vcfFile.getAbsolutePath() + ".tbi"));
            }
            // VCFFileReader can now be initialized.
            this.vcfFileReader = new VCFFileReader(
                    this.vcfFile,
                    new File(this.vcfFile.getAbsolutePath() + ".tbi")
            );
        } catch (Exception e) {
            throw new MusialException("Failed to initialize .vcf file reader for sample " + this.name + "; " + e.getMessage());
        }
    }

    /**
     * Closes and removes the {@link VCFFileReader} of this instance, if imputed.
     */
    public void killVcfFileReader() {
        if (this.vcfFileReader != null) {
            this.vcfFileReader.close();
            this.vcfFileReader = null;
        }
    }

    /**
     * Sets the allele allocation of this sample for feature {@code featureName} to {@code alleleName}.
     *
     * @param featureName {@link Feature#name}.
     * @param alleleName  {@link Form#name}.
     */
    public void setAllele(String featureName, String alleleName) {
        this.allele.put(featureName, alleleName);
    }

    /**
     * Returns the value set for {@code featureName} in {@link #allele} or {@code null} if not set.
     *
     * @param featureName {@link Feature#name} to retrieve {@link Form#name} associated with this sample for.
     */
    public String getAllele(String featureName) {
        return this.allele.getOrDefault(featureName, null);
    }

    /**
     * Sets the proteoform allocation of this sample for feature {@code featureName} to {@code alleleName}.
     *
     * @param featureName    {@link FeatureCoding#name}.
     * @param proteoformName {@link Form#name}.
     */
    public void setProteoform(String featureName, String proteoformName) {
        this.proteoform.put(featureName, proteoformName);
    }

    /**
     * Returns the value set for {@code featureName} in {@link #proteoform} or {@code null} if not set.
     *
     * @param featureName {@link FeatureCoding#name}.
     */
    public String getProteoform(String featureName) {
        return this.proteoform.getOrDefault(featureName, null);
    }

    /**
     * Adds a variant call for the specified feature at the given position.
     *
     * @param featureName The name of the feature.
     * @param pos         The position of the variant call.
     * @param call        A {@link Tuple} representing the primary variant call.
     * @param calls       A {@link List} of {@link Tuple}s representing additional variant calls.
     */
    public void addVariantCall(String featureName, int pos, Tuple<String, Integer> call, List<Tuple<String, Integer>> calls) {
        this.calls.putIfAbsent(featureName, new LinkedHashMap<>(50, 50));
        this.calls.get(featureName).put(
                pos,
                call.a + Constants.FIELD_SEPARATOR_2 + call.b + Constants.FIELD_SEPARATOR_2 + calls.stream().map(c -> c.a + Constants.FIELD_SEPARATOR_1 + c.b).collect(Collectors.joining(","))
        );
    }

    /**
     * Retrieves the variant calls associated with the specified feature.
     *
     * @param featureName The name of the feature.
     * @return A {@link LinkedHashMap} containing variant calls mapped by their positions.
     */
    public LinkedHashMap<Integer, String> getVariantCalls(String featureName) {
        LinkedHashMap<Integer, String> c = this.calls.getOrDefault(featureName, new LinkedHashMap<>());
        LinkedHashMap<Integer, String> r = new LinkedHashMap<>(c.size());
        for (Integer position : c.keySet()) {
            r.put(position, getCallAt(featureName, position));
        }
        return r;
    }

    /**
     * Retrieves the variant call at the specified position for the given feature.
     *
     * @param featureName The name of the feature.
     * @param position    The position of the variant call.
     * @return The variant call at the specified position, or {@link Constants#CALL_INFO_NO_VARIANT} if not present.
     */
    public String getCallAt(String featureName, int position) {
        if (!this.calls.containsKey(featureName) || !this.calls.get(featureName).containsKey(position)) {
            return Constants.CALL_INFO_NO_VARIANT;
        } else {
            String[] fields = this.calls.get(featureName).get(position).split(Constants.FIELD_SEPARATOR_2);
            ArrayList<String> callContents = (ArrayList<String>) Arrays.stream(fields[2].split(",")).map(c -> c.split(Constants.FIELD_SEPARATOR_1)[0]).toList();
            if (fields[0].equals(Constants.CALL_INFO_NO_INFO)) {
                return Constants.CALL_INFO_NO_INFO.repeat(callContents.get(0).length());
            } else if (fields[0].equals(Constants.CALL_INFO_REJECTED)) {
                return Constants.CALL_INFO_REJECTED.repeat(callContents.get(0).length());
            } else {
                return callContents.get(Integer.parseInt(fields[0]));
            }
        }
    }

    /**
     * Retrieves the string representation of the variant call at the specified position for the given feature.
     *
     * @param featureName The name of the feature.
     * @param position    The position of the variant call.
     * @return The string representation of the variant call at the specified position, or {@link Constants#CALL_INFO_NO_VARIANT} if not present.
     */
    public String getCallVariantsStringAt(String featureName, int position) {
        if (!this.calls.containsKey(featureName) || !this.calls.get(featureName).containsKey(position)) {
            return Constants.CALL_INFO_NO_VARIANT;
        } else {
            return this.calls.get(featureName).get(position);
        }
    }

    /**
     * Checks if the variant call at the specified position for the given feature is a reference call.
     *
     * @param featureName The name of the feature.
     * @param position    The position of the variant call.
     * @return {@code true} if the variant call at the specified position is a reference call, otherwise {@code false}.
     */
    public boolean isReferenceCallAt(String featureName, int position) {
        if (!this.calls.containsKey(featureName) || !this.calls.get(featureName).containsKey(position)) {
            return false;
        } else {
            return this.calls.get(featureName).get(position).split(Constants.FIELD_SEPARATOR_2)[0].equals("0");
        }
    }

    /**
     * Checks if the variant call at the specified position for the given feature is ambiguous.
     *
     * @param featureName The name of the feature.
     * @param position    The position of the variant call.
     * @return {@code true} if the variant call at the specified position is ambiguous, otherwise {@code false}.
     */
    public boolean isAmbiguousCallAt(String featureName, int position) {
        if (!this.calls.containsKey(featureName) || !this.calls.get(featureName).containsKey(position)) {
            return false;
        } else {
            return this.calls.get(featureName).get(position).split(Constants.FIELD_SEPARATOR_2)[0].equals(Constants.CALL_INFO_NO_INFO) || this.calls.get(featureName).get(position).split(Constants.FIELD_SEPARATOR_2)[0].equals(Constants.CALL_INFO_REJECTED);
        }
    }

}
