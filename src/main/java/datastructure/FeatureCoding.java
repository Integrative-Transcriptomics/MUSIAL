package datastructure;

import exceptions.MusialException;
import main.Constants;
import utility.SequenceOperations;

import java.util.*;

/**
 * Extension of the {@link Feature} class to represent coding features.
 *
 * @author Simon Hackl
 */
@SuppressWarnings("unused")
public class FeatureCoding extends Feature {

    /**
     * {@link LinkedHashMap} storing all {@link Form} instances, i.e. proteoforms, associated with this feature.
     */
    private final LinkedHashMap<String, Form> proteoforms = new LinkedHashMap<>();
    /**
     * Hierarchical map structure to store variants wrt. proteoforms. The first layer represents the
     * position on the feature. The second layer represents the variant content.
     * <p>
     * Stored positions have to be formatted as X+Y; X is the position on the reference, Y is the possibly
     * inserted position.
     */
    private final TreeMap<Integer, LinkedHashMap<String, VariantInformation>> aminoacidVariants =
            new TreeMap<>(Integer::compare);
    /**
     * The translated reference DNA sequence of the entry.
     */
    private String codingSequence = null;

    /**
     * Constructor of {@link FeatureCoding}.
     * <p>
     * Generates a {@link FeatureCoding} with no protein structure associated.
     *
     * @param name       {@link String} representing the internal name of the reference feature to analyze.
     * @param chromosome {@link String} the name of the reference location (contig, chromosome, plasmid) the feature
     *                   is located on.
     * @param start      {@link Integer} The 1-based indexed starting position of the feature on the reference.
     * @param end        {@link Integer} The 1-based indexed end position of the feature on the reference.
     * @throws MusialException If the specified locus is ambiguous.
     */
    public FeatureCoding(String name, String chromosome, int start, int end, String type) throws MusialException {
        super(name, chromosome, start, end, type);
    }

    /**
     * @return {@link FeatureCoding#codingSequence}
     */
    public String getCodingSequence() {
        return this.codingSequence;
    }

    /**
     * Sets the specified string as the coding sequence of this instance.
     * <br>
     * If the sequence is of nucleotide content, the sequence is translated.
     *
     * @param sequence Content to set as coding sequence.
     * @throws MusialException If the translation of a provided nucleotide sequence fails or no nucleotide or aminoacid sequence was provided.
     */
    public void setCodingSequence(String sequence) throws MusialException {
        if (sequence.matches("[ARNDCEQGHILKMFPSTWYVX*]+")) {
            this.codingSequence = sequence;
        } else if (sequence.matches("[ACGTN]")) {
            this.codingSequence = SequenceOperations.translateNucSequence(sequence, true, true, this.isSense);
        } else {
            throw new MusialException("Specified sequence " + sequence.substring(0, Math.max(sequence.length(), 10)) + "... could not be identified as nucleotide or aminoacid sequence.");
        }
    }

    /**
     * Adds a proteoform, i.e., a {@link Form} instance to {@link FeatureCoding#proteoforms} if it not already exists.
     *
     * @param proteoform {@link Form} representing a proteoform to add.
     */
    public void addProteoform(Form proteoform) {
        if (!this.proteoforms.containsKey(proteoform.name)) {
            this.proteoforms.put(proteoform.name, proteoform);
        }
    }

    /**
     * Removes the {@link Form}, i.e., proteoform, with proteoformName from {@link FeatureCoding#proteoforms} if it exists.
     *
     * @param proteoformName {@link String}, the name of the proteoform to remove.
     */
    public void removeProteoform(String proteoformName) {
        this.proteoforms.remove(proteoformName);
    }

    /**
     * Returns the {@link Form}, i.e., proteoform, stored with proteoformName from {@link FeatureCoding#proteoforms} or null if it does not exist.
     *
     * @param proteoformName {@link String}, the name of the proteoform to query.
     * @return {@link Form} or null.
     */
    public Form getProteoform(String proteoformName) {
        return this.proteoforms.get(proteoformName);
    }

    /**
     * Returns the number of proteoforms stored for this feature.
     *
     * @return Number of proteoforms.
     */
    public int getProteoformCount() {
        int c = 0;
        for (String s : this.proteoforms.keySet()) {
            if (!Objects.equals(s, Constants.REFERENCE_FORM_NAME))
                c += 1;
        }
        return c;
    }

    /**
     * Checks whether a {@link Form}, i.e., proteoform is stored with proteoformName in {@link FeatureCoding#proteoforms}.
     *
     * @param proteoformName {@link String}, the name of the proteoform to check for.
     * @return True if an {@link Form}, i.e., proteoform with name proteoformName exists.
     */
    public boolean hasProteoform(String proteoformName) {
        return this.proteoforms.containsKey(proteoformName);
    }

    /**
     * @return {@link String} {@link Iterator} over all available {@link Form} instances in {@link FeatureCoding#proteoforms}.
     */
    public Iterator<String> getProteoformNameIterator() {
        return this.proteoforms.keySet().iterator();
    }

    /**
     * Adds an aminoacid variant to this {@link FeatureCoding#aminoacidVariants}.
     *
     * @param position The position of the variant.
     * @param alt      The alternative content of the variant.
     * @param ref      The reference content of the variant.
     */
    public void addAminoacidVariant(int position, String alt, String ref) {
        this.aminoacidVariants.putIfAbsent(position, new LinkedHashMap<>());
        this.aminoacidVariants.get(position).putIfAbsent(alt, new VariantInformation(ref));
        // Determine variant type.
        StringBuilder typeBuilder = new StringBuilder();
        if (alt.contains(Constants.DELETION_OR_GAP_STRING))
            typeBuilder.append(Constants.TYPE_AMBIGUOUS_PREFIX);
        if (alt.length() == 1)
            typeBuilder.append(Constants.TYPE_SUBSTITUTION);
        else {
            if (ref.contains(Constants.DELETION_OR_GAP_STRING))
                typeBuilder.append(Constants.TYPE_INSERTION);
            else
                typeBuilder.append(Constants.TYPE_DELETION);
        }
        this.aminoacidVariants.get(position).get(alt).addInfo(Constants.TYPE, typeBuilder.toString());
    }

    /**
     * Returns the {@link VariantInformation} stored at this instances {@link FeatureCoding#aminoacidVariants} at position and content.
     *
     * @param position The position of the variant annotation to access.
     * @param content  The alternative content of the variant annotation to access.
     * @return {@link VariantInformation} object.
     */
    public VariantInformation getAminoacidVariant(int position, String content) {
        if (this.aminoacidVariants.containsKey(position)) {
            return this.aminoacidVariants.get(position).get(content);
        } else {
            return null;
        }
    }

    /**
     * Returns all {@link VariantInformation}s stored at this instances {@link FeatureCoding#aminoacidVariants} at position.
     *
     * @param position The position of the variant annotations to access.
     * @return A map of alternative contents to {@link VariantInformation}s.
     */
    public LinkedHashMap<String, VariantInformation> getAminoacidVariantsAt(int position) {
        return this.aminoacidVariants.getOrDefault(position, null);
    }

    /**
     * Retrieves a set of all variant positions of this {@link FeatureCoding}.
     *
     * @return A navigable sorted map of positions at which a variant occurs in any proteoform of this {@link FeatureCoding}.
     */
    public NavigableSet<Integer> getAminoacidVariantPositions() {
        return this.aminoacidVariants.navigableKeySet();
    }
}