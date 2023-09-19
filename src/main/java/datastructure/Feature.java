package datastructure;

import exceptions.MusialException;
import main.MusialConstants;
import utility.Compression;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentSkipListMap;

/**
 * Container to store representation of a reference sequence location that is subject to analysis.
 * This may be the full genome, a single gene, contigs or plasmids and chromosomes.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.0
 */
@SuppressWarnings("unused")
public class Feature {

    /**
     * The type of this feature; Either coding or non_coding.
     * This value is used to distinguish {@link Feature} from {@link FeatureCoding} instances during parsing an existing MUSIAL storage JSON file.
     */
    protected String type;

    /**
     * The internal name of the feature.
     */
    public final String name;
    /**
     * The chromosome/location of the entry on the reference, i.e., the contig or chromosome the feature is located on.
     */
    public final String chromosome;
    /**
     * The 1-based indexed starting position of the feature.
     */
    public final int start;
    /**
     * The 1-based indexed end position of the feature.
     */
    public final int end;
    /**
     * Indicates if the feature is located on the sense strand.
     */
    public final boolean isSense;
    /**
     * {@link LinkedHashMap} storing all {@link Form} instances, i.e. alleles, associated with this feature.
     */
    protected final LinkedHashMap<String, Form> alleles = new LinkedHashMap<>();
    /**
     * A {@link LinkedHashMap} yielding any meta-information {@link String} key-value pairs about this {@link Feature}.
     */
    protected final LinkedHashMap<String, String> annotations = new LinkedHashMap<>();
    /**
     * Hierarchical map structure to store variants wrt. the nucleotide sequences. The first layer represents the
     * position on the chromosome. The second layer represents the variant content.
     */
    protected final ConcurrentSkipListMap<Integer, ConcurrentSkipListMap<String, VariantAnnotation>>
            nucleotideVariants =
            new ConcurrentSkipListMap<>(Integer::compare);

    /**
     * Constructor of {@link Feature}.
     *
     * @param name       {@link String} representing the internal name of the reference feature to analyze.
     * @param chromosome {@link String} the name of the reference location (contig, chromosome, plasmid) the feature
     *                   is located on.
     * @param start      {@link Integer} The 1-based indexed starting position of the feature on the reference.
     * @param end        {@link Integer} The 1-based indexed end position of the feature on the reference.
     * @throws MusialException If the specified locus is ambiguous.
     */
    public Feature(String name, String chromosome, int start, int end, String type) throws MusialException {
        this.name = name;
        this.chromosome = chromosome;
        this.type = type;
        if (start >= 0 && end > 0 && (end > start)) {
            // CASE: Feature is on sense strand.
            this.isSense = true;
            this.start = start + 1; // +1 is an artifact from bug in GFFParser class.
            this.end = end;
        } else if (start < 0 && end <= 0 && (end < start)) {
            // CASE: Feature is on anti-sense strand.
            this.isSense = false;
            this.start = -start + 1; // +1 is an artifact from bug in GFFParser class.
            this.end = -end;
        } else {
            throw new MusialException(
                    "Failed to add  feature due to faulty position data (start, end)\t" + start + ", " + end);
        }
    }

    /**
     * Checks whether this instance of {@link Feature} is an instance of {@link FeatureCoding}, i.e., declared as a coding feature.
     *
     * @return True if this is instance of {@link FeatureCoding}.
     */
    public boolean isCoding() {
        return this.type.equals("coding");
    }

    /**
     * Adds an allele, i.e., a {@link Form} instance to {@link Feature#alleles} if it not already exists.
     *
     * @param allele {@link Form} representing an allele to add.
     */
    public void addAllele(Form allele) {
        if (!this.alleles.containsKey(allele.name)) {
            this.alleles.put(allele.name, allele);
        }
    }

    /**
     * Removes the {@link Form}, i.e., allele, with alleleName from {@link Feature#alleles} if it exists.
     *
     * @param alleleName {@link String}, the name of the allele to remove.
     */
    public void removeAllele(String alleleName) {
        this.alleles.remove(alleleName);
    }

    /**
     * Returns the {@link Form}, i.e., allele, stored with alleleName from {@link Feature#alleles} or null if it does not exist.
     *
     * @param alleleName {@link String}, the name of the allele to query.
     * @return {@link Form} or null.
     */
    public Form getAllele(String alleleName) {
        return this.alleles.get(alleleName);
    }

    /**
     * Returns the number of alleles stored for this feature.
     *
     * @return Number of alleles.
     */
    public int getAlleleCount() {
        return this.alleles.size();
    }

    /**
     * Checks whether a {@link Form}, i.e., allele is stored with alleleName in {@link Feature#alleles}.
     *
     * @param alleleName {@link String}, the name of the allele to check for.
     * @return True if an {@link Form}, i.e., allele with name alleleName exists.
     */
    public boolean hasAllele(String alleleName) {
        return this.alleles.containsKey(alleleName);
    }

    /**
     * @return {@link String} {@link Iterator} over all available {@link Form} instances in {@link Feature#alleles}.
     */
    public Iterator<String> getAlleleNameIterator() {
        return this.alleles.keySet().iterator();
    }

    /**
     * Adds an annotation specified by the provided key and value to this instances {@link Feature#annotations}.
     *
     * @param key   {@link String} key of the annotation.
     * @param value {@link String} value of the annotation.
     */
    public void addAnnotation(String key, String value) {
        this.annotations.put(key, value);
    }

    /**
     * Removes the annotation stored at key from this instances {@link Feature#annotations}.
     *
     * @param key {@link String} key value of the annotation to remove.
     */
    public void removeAnnotation(String key) {
        this.annotations.remove(key);
    }

    /**
     * Returns the value stored under key in {@link Feature#annotations}.
     *
     * @param key {@link String} key of the annotation.
     * @return The value of the annotation.
     */
    public String getAnnotation(String key) {
        return this.annotations.get(key);
    }

    /**
     * Returns all annotations of this instance.
     *
     * @return Set of all annotations.
     */
    public Set<Map.Entry<String, String>> getAnnotations() {
        return this.annotations.entrySet();
    }

    /**
     * Adds a nucleotide variant to this {@link Feature#nucleotideVariants}.
     *
     * @param position The position of the variant.
     * @param content  The alternative content of the variant.
     */
    public void addNucleotideVariant(int position, String content) {
        if (!this.nucleotideVariants.containsKey(position)) {
            this.nucleotideVariants.put(position, new ConcurrentSkipListMap<>());
        }
        if (!this.nucleotideVariants.get(position).containsKey(content)) {
            this.nucleotideVariants.get(position).put(content, new VariantAnnotation());
        }
    }

    /**
     * Returns the {@link VariantAnnotation} stored at this instances {@link Feature#nucleotideVariants} at position and content.
     *
     * @param position The position of the variant annotation to access.
     * @param content  The alternative content of the variant annotation to access.
     * @return {@link VariantAnnotation} object.
     */
    public VariantAnnotation getNucleotideVariantAnnotation(int position, String content) {
        if (this.nucleotideVariants.containsKey(position)) {
            return this.nucleotideVariants.get(position).get(content);
        } else {
            return null;
        }
    }

    /**
     * Returns all {@link VariantAnnotation}s stored at this instances {@link Feature#nucleotideVariants} at position.
     *
     * @param position The position of the variant annotations to access.
     * @return A map of alternative contents to {@link VariantAnnotation}s.
     */
    public ConcurrentSkipListMap<String, VariantAnnotation> getNucleotideVariants(int position) {
        return this.nucleotideVariants.getOrDefault(position, null);
    }

    /**
     * Returns a map representation of all variants of a single allele of this {@link Feature}.
     *
     * @param alleleName The name of the allele.
     * @return Map of variable positions to alternative base contents.
     * @throws IOException If decoding of the internal variant representation fails.
     */
    public TreeMap<Integer, String> getNucleotideVariants(String alleleName) throws IOException {
        TreeMap<Integer, String> alleleVariants = new TreeMap<>();
        if (!alleleName.equals(MusialConstants.REFERENCE_ID)) {
            String[] variantFields;
            for (String variant : Compression.brotliDecodeString(this.alleles.get(alleleName).getAnnotation(MusialConstants.VARIANTS)).split(MusialConstants.FIELD_SEPARATOR_2)) {
                variantFields = variant.split(MusialConstants.FIELD_SEPARATOR_1);
                alleleVariants.put(Integer.valueOf(variantFields[0]), variantFields[1]);
            }
        }
        return alleleVariants;
    }

    /**
     * Retrieves a set of all variant positions of this {@link Feature}.
     *
     * @return A navigable sorted map of positions at which a variant occurs in any allele of this {@link Feature}.
     */
    public NavigableSet<Integer> getNucleotideVariantPositions() {
        return this.nucleotideVariants.keySet();
    }

}