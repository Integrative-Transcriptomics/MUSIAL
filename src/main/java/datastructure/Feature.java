package datastructure;

import main.Constants;
import org.apache.commons.lang3.tuple.Pair;

import java.util.*;

/**
 * Container to store representation of a reference sequence location that is subject to analysis.
 * This may be the full genome, a single gene, contigs or plasmids and chromosomes.
 *
 * @author Simon Hackl
 */
@SuppressWarnings("unused")
public class Feature extends InfoContainer {

    /**
     * The type of this feature; Either coding or non_coding.
     * This value is used to distinguish {@link Feature} from {@link FeatureCoding} instances during parsing an existing MUSIAL storage JSON file.
     */
    public final String type;

    /**
     * The internal name of the feature.
     */
    public final String name;
    /**
     * The chromosome/location of the entry on the reference, i.e., the contig or chromosome the feature is located on.
     */
    public final String contig;
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
    private final LinkedHashMap<String, Form> alleles = new LinkedHashMap<>();
    /**
     * Hierarchical map structure to store variants wrt. the nucleotide sequences. The first layer represents the
     * position on the chromosome. The second layer represents the variant content.
     */
    private final TreeMap<Integer, LinkedHashMap<String, VariantInformation>>
            nucleotideVariants =
            new TreeMap<>(Integer::compare);

    /**
     * Constructor of {@link Feature}.
     *
     * @param name   {@link String} representing the internal name of the reference feature to analyze.
     * @param contig {@link String} the name of the reference location (contig, chromosome, plasmid) the feature
     *               is located on.
     * @param start  {@link Integer} The 1-based indexed starting position of the feature on the reference.
     * @param end    {@link Integer} The 1-based indexed end position of the feature on the reference.
     * @param strand The strand orientation (+/-) of the coding feature.
     * @param type   The type (coding or non-coding) of this feature.
     */
    public Feature(String name, String contig, int start, int end, char strand, String type) {
        super();
        this.name = name;
        this.contig = contig;
        this.type = type;
        this.start = start;
        this.end = end;
        this.isSense = (strand == '+');
        this.addAllele(new Form(Constants.REFERENCE_FORM_NAME, ""));
    }

    /**
     * Return whether this instance of {@link Feature} is declared as a coding feature.
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
        int c = 0;
        for (String s : this.alleles.keySet()) {
            if (!Objects.equals(s, Constants.REFERENCE_FORM_NAME))
                c += 1;
        }
        return c;
    }

    /**
     * Checks whether a {@link Form}, i.e., allele is stored with alleleName in {@link Feature#alleles}.
     *
     * @param alleleName {@link String}, the name of the allele to check for.
     * @return True if an {@link Form}, i.e., allele with name alleleName exists.
     */
    @SuppressWarnings("BooleanMethodIsAlwaysInverted")
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
     * Adds a nucleotide variant to this {@link Feature#nucleotideVariants}.
     *
     * @param position The position of the variant.
     * @param alt      The alternative content of the variant.
     * @param ref      The reference content of the variant.
     */
    public void addNucleotideVariant(int position, String alt, String ref) {
        this.nucleotideVariants.putIfAbsent(position, new LinkedHashMap<>());
        this.nucleotideVariants.get(position).putIfAbsent(alt, new VariantInformation(ref));
        // Determine variant type.
        String type;
        if (alt.contains(Constants.ANY_NUCLEOTIDE_STRING))
            type = Constants.TYPE_AMBIGUOUS;
        else if (alt.length() == 1)
            type = Constants.TYPE_SUBSTITUTION;
        else {
            if (ref.contains(Constants.DELETION_OR_GAP_STRING))
                type = Constants.TYPE_INSERTION;
            else
                type = Constants.TYPE_DELETION;
        }
        this.nucleotideVariants.get(position).get(alt).addInfo(Constants.TYPE, type);
    }

    /**
     * Returns the {@link VariantInformation} stored at this instances {@link Feature#nucleotideVariants} at position and content.
     *
     * @param pos The position of the variant annotation to access.
     * @param alt The alternative content of the variant annotation to access.
     * @return {@link VariantInformation} object.
     */
    public VariantInformation getNucleotideVariant(int pos, String alt) {
        if (this.nucleotideVariants.containsKey(pos))
            return this.nucleotideVariants.get(pos).getOrDefault(alt, null);
        else
            return null;
    }

    /**
     * Returns {@link VariantInformation}s stored at this instances {@link Feature#nucleotideVariants} at position.
     *
     * @param position    The position of the variant annotations to access.
     * @param onlyPrimary Whether to retrieve only primary or all variants.
     * @return A map of alternative contents to {@link VariantInformation}s.
     */
    public LinkedHashMap<String, VariantInformation> getNucleotideVariantsAt(int position, boolean onlyPrimary) {
        if (onlyPrimary) {
            LinkedHashMap<String, VariantInformation> variants = new LinkedHashMap<>();
            this.nucleotideVariants.getOrDefault(position, null)
                    .entrySet().stream().filter(entry -> Boolean.parseBoolean(entry.getValue().getInfo(Constants.VARIANT_INFO_PRIMARY)))
                    .forEach(entry -> variants.put(entry.getKey(), entry.getValue()));
            return variants;
        } else
            return this.nucleotideVariants.getOrDefault(position, null);
    }

    /**
     * Retrieves a set of all variant positions of this {@link Feature}.
     *
     * @return A navigable sorted map of positions at which a variant occurs in any allele of this {@link Feature}.
     */
    public NavigableSet<Integer> getNucleotideVariantPositions() {
        return this.nucleotideVariants.navigableKeySet();
    }

    /**
     * Retrieves all variants stored for the specified allele/{@link Form} with name {@code alleleName}.
     *
     * @param alleleName The name of the allele of which variants shall be retrieved.
     * @return Navigable map structure with variants of specified allele. Might be empty.
     */
    public TreeMap<Integer, Pair<String, VariantInformation>> getNucleotideVariantsOf(String alleleName) {
        TreeMap<Integer, Pair<String, VariantInformation>> variants = new TreeMap<>();
        int position;
        String alt;
        for (String variantEntry : this.alleles.get(alleleName).variants.split(Constants.ENTRY_SEPARATOR)) {
            String[] variantFields = variantEntry.split(Constants.FIELD_SEPARATOR);
            position = Integer.parseInt(variantFields[0]);
            alt = variantFields[1];
            if (this.nucleotideVariants.containsKey(position) && this.nucleotideVariants.get(position).containsKey(alt))
                variants.put(
                        position,
                        Pair.of(alt, this.nucleotideVariants.get(position).get(alt))
                );
        }
        return variants;
    }

}