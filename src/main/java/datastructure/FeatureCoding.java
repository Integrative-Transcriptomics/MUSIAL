package datastructure;

import exceptions.MusialException;
import main.Musial;
import main.MusialConstants;
import org.apache.commons.io.IOUtils;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;
import utility.Bio;
import utility.Compression;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.NavigableSet;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentSkipListMap;

/**
 * Extension of the {@link Feature} class to represent coding features.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.2
 */
@SuppressWarnings("unused")
public class FeatureCoding extends Feature {

    /**
     * The translated reference DNA sequence of the entry.
     */
    private String codingSequence = null;
    /**
     * The protein's structure in .pdb format, if this feature is a coding sequence.
     */
    private String structure = null;
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
    private final ConcurrentSkipListMap<Integer, ConcurrentSkipListMap<String, VariantAnnotation>> aminoacidVariants =
            new ConcurrentSkipListMap<>(Integer::compare);

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
            this.codingSequence = Bio.translateNucSequence(sequence, true, true, this.isSense);
        } else {
            throw new MusialException("Specified sequence " + sequence.substring(0, Math.max(sequence.length(), 10)) + "... could not be identified as nucleotide or aminoacid sequence.");
        }
    }

    /**
     * @return {@link FeatureCoding#codingSequence}
     */
    public String getCodingSequence() {
        return this.codingSequence;
    }

    /**
     * Sets the protein structure from a pdb format file.
     * <br>
     * The structure is stored as brotli compressed string.
     *
     * @param pdbFile The file object to parse the structure from.
     * @throws IOException If reading of the pdb file or compression of the string content fails.
     */
    public void setStructure(File pdbFile) throws IOException {
        try {
            System.setOut(Musial.EMPTY_STREAM);
            Structure pdbStructure = new PDBFileReader().getStructure(pdbFile);
            this.structure = Compression.brotliEncodeString(pdbStructure.toPDB());
        } finally {
            System.setOut(Musial.ORIGINAL_OUT_STREAM);
        }
    }

    /**
     * @return {@link Structure} re-built from the encoded {@link FeatureCoding#structure}.
     * @throws IOException If building of the structure or decompression fails.
     */
    public Structure getStructure() throws IOException {
        String pdbString = Compression.brotliDecodeString(this.structure);
        return new PDBFileReader().getStructure(IOUtils.toInputStream(pdbString, StandardCharsets.UTF_8));
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
        return this.proteoforms.size();
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
     * @param content  The alternative content of the variant.
     */
    public void addAminoacidVariant(int position, String content) {
        if (!this.aminoacidVariants.containsKey(position)) {
            this.aminoacidVariants.put(position, new ConcurrentSkipListMap<>());
        }
        if (!this.aminoacidVariants.get(position).containsKey(content)) {
            this.aminoacidVariants.get(position).put(content, new VariantAnnotation());
        }
    }

    /**
     * Returns the {@link VariantAnnotation} stored at this instances {@link FeatureCoding#aminoacidVariants} at position and content.
     *
     * @param position The position of the variant annotation to access.
     * @param content  The alternative content of the variant annotation to access.
     * @return {@link VariantAnnotation} object.
     */
    public VariantAnnotation getAminoacidVariantAnnotation(int position, String content) {
        if (this.aminoacidVariants.containsKey(position)) {
            return this.aminoacidVariants.get(position).get(content);
        } else {
            return null;
        }
    }

    /**
     * Returns all {@link VariantAnnotation}s stored at this instances {@link FeatureCoding#aminoacidVariants} at position.
     *
     * @param position The position of the variant annotations to access.
     * @return A map of alternative contents to {@link VariantAnnotation}s.
     */
    public ConcurrentSkipListMap<String, VariantAnnotation> getAminoacidVariants(int position) {
        return this.aminoacidVariants.getOrDefault(position, null);
    }

    /**
     * Returns a map representation of all variants of a single proteoform of this {@link FeatureCoding}.
     *
     * @param proteoformName The name of the proteoform.
     * @return Map of variable positions to alternative base contents.
     * @throws IOException If decoding of the internal variant representation fails.
     */
    public TreeMap<Integer, String> getAminoacidVariants(String proteoformName) throws IOException {
        TreeMap<Integer, String> aminoacidVariants = new TreeMap<>();
        if (this.proteoforms.get(proteoformName).hasAnnotation(MusialConstants.PROTEOFORM_DIFFERENTIAL_SEQUENCE)) {
            return null; // FIXME
        }
        if (!proteoformName.equals(MusialConstants.REFERENCE_ID)) {
            String[] variantFields;
            for (String variant : Compression.brotliDecodeString(this.proteoforms.get(proteoformName).getAnnotation(MusialConstants.VARIANTS)).split(MusialConstants.FIELD_SEPARATOR_2)) {
                variantFields = variant.split(MusialConstants.FIELD_SEPARATOR_1);
                aminoacidVariants.put(Integer.valueOf(variantFields[0]), variantFields[1]);
            }
        }
        return aminoacidVariants;
    }

    /**
     * Retrieves a set of all variant positions of this {@link FeatureCoding}.
     *
     * @return A navigable sorted map of positions at which a variant occurs in any proteoform of this {@link FeatureCoding}.
     */
    public NavigableSet<Integer> getAminoacidVariantPositions() {
        return this.aminoacidVariants.keySet();
    }

}