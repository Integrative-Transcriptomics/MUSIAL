package utility;

import com.google.common.base.Splitter;
import datastructure.ReferenceFeatureEntry;
import datastructure.VariantContent;
import datastructure.VariantContentTable;
import exceptions.MusialBioException;
import org.biojava.nbio.structure.*;
import org.javatuples.Pair;
import org.javatuples.Triplet;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.util.*;
import java.util.concurrent.ConcurrentSkipListMap;

/**
 * Comprises static methods used to solve specifically bioinformatics problems.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class Bio {

    /**
     * Character to indicate a gap/missing information regarding a nucleotide and amino acid sequence.
     */
    public static final char NO_MATCH_CHAR = '#';
    /**
     * Character to indicate a truncation of a sequence.
     */
    public static final char TRUNCATED_CHAR = '*';
    /**
     * Character used as a separator in padded sequences.
     */
    public static final char SEPARATOR_CHAR = '/';
    /**
     * One-letter code character used to indicate a translated stop codon.
     */
    public static final char TERMINATION_AA1 = 'Z';
    /**
     * Three-letter code String used to indicate a translated stop codon.
     */
    public static final String TERMINATION_AA3 = "TER";
    /**
     * One-letter code character used to indicate the translation of a codon containing an N character (i.e. unknown or
     * any nucleotide).
     */
    public static final char UNKNOWN_AA1 = 'X';
    /**
     * Three-letter code string used to indicate the translation of a codon containing an N character (i.e. unknown or
     * any nucleotide).
     */
    public static final String UNKNOWN_AA3 = "ANY";
    /**
     * One-letter code character used to indicate the translation of an incomplete codon.
     */
    public static final char INCOMPLETE_AA1 = 'U';
    /**
     * Three-letter code string used to indicate the translation of an incomplete codon.
     */
    public static final String INCOMPLETE_AA3 = "INC";
    /**
     * Three-letter code string used to indicate the deletion of an amino-acid.
     */
    public static final String DELETION_AA3 = "DEL";
    /**
     * Placeholder string to indicate no content.
     */
    public static final String NONE = "None";
    /**
     * Hash map mapping nucleotide codons to single-letter amino acids. Stop codons map to no character.
     */
    public static final HashMap<String, String> CODON_MAP = new HashMap<>() {{
        put("TTT", "F");
        put("TTC", "F");
        put("TTA", "L");
        put("TTG", "L");
        put("CTT", "L");
        put("CTC", "L");
        put("CTA", "L");
        put("CTG", "L");
        put("ATT", "I");
        put("ATC", "I");
        put("ATA", "I");
        put("ATG", "M");
        put("GTT", "V");
        put("GTC", "V");
        put("GTA", "V");
        put("GTG", "V");
        put("TCT", "S");
        put("TCC", "S");
        put("TCA", "S");
        put("TCG", "S");
        put("CCT", "P");
        put("CCC", "P");
        put("CCA", "P");
        put("CCG", "P");
        put("ACT", "T");
        put("ACC", "T");
        put("ACA", "T");
        put("ACG", "T");
        put("GCT", "A");
        put("GCC", "A");
        put("GCA", "A");
        put("GCG", "A");
        put("TAT", "Y");
        put("TAC", "Y");
        put("TAA", "");
        put("TAG", "");
        put("CAT", "H");
        put("CAC", "H");
        put("CAA", "Q");
        put("CAG", "Q");
        put("AAT", "N");
        put("AAC", "N");
        put("AAA", "K");
        put("AAG", "K");
        put("GAT", "D");
        put("GAC", "D");
        put("GAA", "E");
        put("GAG", "E");
        put("TGT", "C");
        put("TGC", "C");
        put("TGA", "");
        put("TGG", "W");
        put("CGT", "R");
        put("CGC", "R");
        put("CGA", "R");
        put("CGG", "R");
        put("AGT", "S");
        put("AGC", "S");
        put("AGA", "R");
        put("AGG", "R");
        put("GGT", "G");
        put("GGC", "G");
        put("GGA", "G");
        put("GGG", "G");
    }};
    /**
     * Hash map mapping amino-acid one to three letter codes.
     */
    public static final HashMap<String, String> AA1TO3 = new HashMap<>() {{
        put("A", "ALA");
        put("R", "ARG");
        put("N", "ASN");
        put("D", "ASP");
        put("C", "CYS");
        put("E", "GLU");
        put("Q", "GLN");
        put("G", "GLY");
        put("H", "HIS");
        put("I", "ILE");
        put("L", "LEU");
        put("K", "LYS");
        put("M", "MET");
        put("F", "PHE");
        put("P", "PRO");
        put("S", "SER");
        put("T", "THR");
        put("W", "TRP");
        put("Y", "TYR");
        put("V", "VAL");
        put(String.valueOf(UNKNOWN_AA1), UNKNOWN_AA3);
        put(String.valueOf(TERMINATION_AA1), TERMINATION_AA3);
    }};
    /**
     * Hash map mapping amino-acid one letter symbols to indices used to access substitution matrix values.
     */
    public static final HashMap<Character, Integer> AA1_BLOSUM80_INDEX = new HashMap<>() {{
        put('A', 0);
        put('R', 1);
        put('N', 2);
        put('D', 3);
        put('C', 4);
        put('Q', 5);
        put('E', 6);
        put('G', 7);
        put('H', 8);
        put('I', 9);
        put('L', 10);
        put('K', 11);
        put('M', 12);
        put('F', 13);
        put('P', 14);
        put('S', 15);
        put('T', 16);
        put('W', 17);
        put('Y', 18);
        put('V', 19);
        put(UNKNOWN_AA1, 20);
        put(TERMINATION_AA1, 21);
    }};
    /**
     * The BLOSUM80 matrix used for alignment computation.
     * <p>
     * Includes values for 'X'; any amino-acid and 'Z'; termination.
     */
    private static final int[][] BLOSUM80 = {
            {5, -2, -2, -2, -1, -1, -1, 0, -2, -2, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0, -1, -1},
            {-2, 6, -1, -2, -4, 1, -1, -3, 0, -3, -3, 2, -2, -4, -2, -1, -1, -4, -3, -3, -1, -1},
            {-2, -1, 6, 1, -3, 0, -1, -1, 0, -4, -4, 0, -3, -4, -3, 0, 0, -4, -3, -4, -1, -1},
            {-2, -2, 1, 6, -4, -1, 1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4, -1, -1},
            {-1, -4, -3, -4, 9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -1, -1},
            {-1, 1, 0, -1, -4, 6, 2, -2, 1, -3, -3, 1, 0, -4, -2, 0, -1, -3, -2, -3, -1, -1},
            {-1, -1, -1, 1, -5, 2, 6, -3, 0, -4, -4, 1, -2, -4, -2, 0, -1, -4, -3, -3, -1, -1},
            {0, -3, -1, -2, -4, -2, -3, 6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -1, -1},
            {-2, 0, 0, -2, -4, 1, 0, -3, 8, -4, -3, -1, -2, -2, -3, -1, -2, -3, 2, -4, -1, -1},
            {-2, -3, -4, -4, -2, -3, -4, -5, -4, 5, 1, -3, 1, -1, -4, -3, -1, -3, -2, 3, -1, -1},
            {-2, -3, -4, -5, -2, -3, -4, -4, -3, 1, 4, -3, 2, 0, -3, -3, -2, -2, -2, 1, -1, -1},
            {-1, 2, 0, -1, -4, 1, 1, -2, -1, -3, -3, 5, -2, -4, -1, -1, -1, -4, -3, -3, -1, -1},
            {-1, -2, -3, -4, -2, 0, -2, -4, -2, 1, 2, -2, 6, 0, -3, -2, -1, -2, -2, 1, -1, -1},
            {-3, -4, -4, -4, -3, -4, -4, -4, -2, -1, 0, -4, 0, 6, -4, -3, -2, 0, 3, -1, -1, -1},
            {-1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4, 8, -1, -2, -5, -4, -3, -1, -1},
            {1, -1, 0, -1, -2, 0, 0, -1, -1, -3, -3, -1, -2, -3, -1, 5, 1, -4, -2, -2, -1, -1},
            {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2, 1, 5, -4, -2, 0, -1, -1},
            {-3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2, 0, -5, -4, -4, 11, 2, -3, -1, -1},
            {-2, -3, -3, -4, -3, -2, -3, -4, 2, -2, -2, -3, -2, 3, -4, -2, -2, 2, 7, -2, -1, -1},
            {0, -3, -4, -4, -1, -3, -3, -4, -4, 3, 1, -3, 1, -1, -3, -2, 0, -3, -2, 4, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1}
    };
    /**
     * Gap opening penalty used for amino-acid alignment computation.
     */
    private static final int GAP_OPEN_PENALTY = 12;
    /**
     * Gap extension penalty used for amino-acid alignment computation.
     */
    private static final int GAP_EXTENSION_PENALTY = 11;

    /**
     * Constructs and returns a full length (sense orientation) nucleotide sequence as {@link String} for the specified
     * sample by using the specified variants table and reference feature entry.
     * <p>
     * The boolean parameters `deletionsAsGaps` and `addInsertions` indicate whether deletions should be
     * inserted as gaps and insertions of the sample should be included into the constructed sequence.
     *
     * @param table                 {@link VariantContentTable} from which variant contents are retrieved.
     * @param referenceFeatureEntry {@link ReferenceFeatureEntry} internal representation of the reference feature.
     * @param sampleName            {@link String} specifying the sample which variant contents should be incorporated into the
     *                              sequence.
     * @param deletionsAsGaps       {@link Boolean} whether deletions of the specified sample should be considered, i.e. if the
     *                              sample shows a deletion on a specific position this position is skipped for
     *                              constructing the sequence.
     * @param addInsertions         {@link Boolean} whether insertions of the specified sample should be considered.
     * @return {@link String} representing the extracted nucleotide sequence.
     */
    public static String getSequenceFromTable(VariantContentTable table, ReferenceFeatureEntry referenceFeatureEntry,
                                              String sampleName, boolean deletionsAsGaps, boolean addInsertions) {
        if (sampleName.equals("Reference")) {
            return referenceFeatureEntry.getReferenceSequence();
        } else {
            String[] sequence = referenceFeatureEntry.getReferenceSequence().split("");
            String sequenceContent;
            Iterator<VariantContent> contentIterator;
            String referenceFeatureName = referenceFeatureEntry.name;
            int featureStart = referenceFeatureEntry.locationStart;
            Iterator<String> variantPositionIterator = table.getVariantPositions(referenceFeatureEntry.name);
            while (variantPositionIterator.hasNext()) {
                String position = variantPositionIterator.next();
                contentIterator = table.getContent(referenceFeatureName, sampleName, position);
                VariantContent variantContent;
                if (contentIterator != null && contentIterator.hasNext()) {
                    variantContent = contentIterator.next();
                    if (variantContent != null && variantContent.MFA) {
                        int positionValue;
                        int sequenceIndex;
                        if (position.contains("+") && addInsertions) {
                            // CASE: Insertion, add symbols to sequence.
                            positionValue = Integer.parseInt(position.split("\\+")[0]);
                            sequenceIndex = positionValue - featureStart;
                            sequenceContent = sequence[sequenceIndex];
                            sequenceContent += variantContent.content;
                            sequence[sequenceIndex] = sequenceContent;
                        } else if (variantContent.content == VariantContent.DELETION) {
                            // CASE: Deletion, remove symbols from sequence.
                            positionValue = Integer.parseInt(position);
                            sequenceIndex = positionValue - featureStart;
                            if (deletionsAsGaps) {
                                sequence[sequenceIndex] = String.valueOf(VariantContent.DELETION);
                            } else {
                                sequence[sequenceIndex] = "";
                            }
                        } else {
                            if (position.contains("+") && !addInsertions) {
                                continue;
                            }
                            // CASE: SNV, exchange symbol in sequence.
                            positionValue = Integer.parseInt(position);
                            sequenceIndex = positionValue - featureStart;
                            sequence[sequenceIndex] = String.valueOf(variantContent.content);
                        }
                    }
                }
            }
            return String.join("", sequence);
        }
    }

    /**
     * Constructs and returns an iterator over the position annotated full length nucleotide
     * sequence for the specified sample by using the specified variants table and reference feature entry.
     * <p>
     * This method respects the orientation of the passed `referenceFeatureEntry`.
     *
     * @param table                 {@link VariantContentTable} from which variant contents are retrieved.
     * @param referenceFeatureEntry {@link ReferenceFeatureEntry} internal representation of the reference feature.
     * @param sampleName            {@link String} specifying the sample which variant contents should be incorporated into the
     * @return {@link Iterator} over the position annotated reconstructed nucleotide sequence.
     */
    public static Iterator<Pair<String, Character>> getPositionAnnotatedSequenceIteratorFromTable(VariantContentTable table,
                                                                                                  ReferenceFeatureEntry referenceFeatureEntry,
                                                                                                  String sampleName) throws MusialBioException {
        ArrayList<Pair<String, Character>> positionAnnotatedSequenceList = new ArrayList<>();
        String sequence = getSequenceFromTable(table, referenceFeatureEntry, sampleName, true, false);
        if ( !referenceFeatureEntry.isSense ) {
            sequence = getReverseComplement( sequence );
        }
        char[] sequenceArray = sequence.toCharArray();
        HashMap<Integer, String> insertions = getInsertionsOfSample(table, referenceFeatureEntry.name, sampleName,
                table.getAllPositionsSet(referenceFeatureEntry.name), referenceFeatureEntry.isSense);
        char[] insertion;
        int startingPosition = referenceFeatureEntry.isSense ? referenceFeatureEntry.locationStart :
                -referenceFeatureEntry.locationEnd;
        for (int i = 0; i < sequenceArray.length; i++) {
            positionAnnotatedSequenceList.add(new Pair<>(String.valueOf(startingPosition + i), sequenceArray[i]));
            if (insertions.containsKey(startingPosition + i)) {
                insertion = insertions.get(startingPosition + i).toCharArray();
                for (int j = 0; j < insertion.length; j++) {
                    positionAnnotatedSequenceList.add(new Pair<>(startingPosition + i + "+" + (j + 1),
                            insertion[j]));
                }
            }
        }
        return positionAnnotatedSequenceList.iterator();
    }

    /**
     * Constructs and returns a {@link HashMap} mapping {@link String} representations of genome positions to
     * {@link String}s representing inserted nucleotide sequences with regard to one analyzed sample.
     * <p>
     * The positions are thereby chosen as the position after which the respective insertion is observed.
     * <p>
     * FIXME: In the case of sense orientation and insertions at the very start and anti-sense orientation and
     * insertions at the very end, no position can be selected to assign the insertion to.
     *
     * @param table                {@link VariantContentTable} from which variant contents are retrieved.
     * @param referenceFeatureName {@link String} specifying the reference feature to which regard the insertions are
     *                             extracted.
     * @param sampleName           {@link String} specifying the sample for which insertions shall be extracted.
     * @param positions            {@link Iterator<String>} specifying which positions are scanned for insertions; this
     *                             should be `table.getAllPositionsIterator(referenceFeatureName)` in order
     *                             to obtain a biologically meaningful result.
     * @param asSense              {@link Boolean} specifying whether the insertions should be extracted with respect to the sense
     *                             or anti-sense strand; i.e. if set to false the stored insertions are the reverse
     *                             complement sequences of the determined insertions and the position is set to -
     *                             (i+1) if i is the position of the insertion in sense orientation.
     * @return {@link HashMap} mapping genome positions to inserted nucleotide sequences.
     * @throws MusialBioException If, during reverse complement computation in the case of anti-sense representation,
     *                            any base can not be translated.
     */
    public static HashMap<Integer, String> getInsertionsOfSample(VariantContentTable table, String referenceFeatureName,
                                                                 String sampleName, NavigableSet<String> positions,
                                                                 boolean asSense) throws MusialBioException {
        HashMap<Integer, String> insertions = new HashMap<>();
        StringBuilder insertionBuilder = new StringBuilder();
        String insertionPosition = null;
        for (String position : positions) {
            if (position.contains("+")) {
                insertionPosition = position;
                Iterator<VariantContent> sampleContentIterator = table.getContent(referenceFeatureName, sampleName, position);
                VariantContent sampleContent = null;
                if (sampleContentIterator != null && sampleContentIterator.hasNext()) {
                    sampleContent = sampleContentIterator.next();
                }
                if (sampleContent != null && sampleContent.MFA) {
                    insertionBuilder.append(sampleContent.content);
                }
            } else {
                if (insertionBuilder.length() != 0 && insertionPosition != null) {
                    String insertionSequence = insertionBuilder.toString();
                    if (asSense) {
                        insertions.put(
                                Integer.valueOf(insertionPosition.substring(0, insertionPosition.indexOf("+"))),
                                insertionSequence
                        );
                    } else {
                        insertions.put(
                                -(Integer.parseInt(insertionPosition.substring(0, insertionPosition.indexOf("+"))) + 1),
                                getReverseComplement(insertionSequence));
                    }
                    insertionBuilder = new StringBuilder();
                }
            }
        }
        return insertions;
    }

    /**
     * Matches a amino-acid and nucleotide sequence.
     * <p>
     * This is done by translating the nucleotide sequence, for every possible frameshift if applicable, into an amino
     * acid sequence and aligning the translated nucleotide sequence with the amino-acid sequence globally.
     * <p>
     * The alignment is used to deduce a mapping of nucleotide codons to the amino acid sequence.
     *
     * @param aminoAcidSequence  {@link String} representing a amino-acid sequence.
     * @param nucleotideSequence {@link String} representing a nucleotide sequence.
     * @return {@link ArrayList<String>} comprising the padded amino acid (at index 0) and nucleotide sequence (at
     * index 1) that were derived as the best alignment.
     * @throws MusialBioException If an error occurs during the translation of the nucleotide sequence, i.e. any codon
     *                            with a length different than three is detected or any codon can not be translated.
     */
    public static Pair<String, String> getNucleotideToAminoAcidSequenceMatching(String aminoAcidSequence,
                                                                                String nucleotideSequence)
            throws MusialBioException {
        // Initialize map to store matched sequences and variable to store matching results.
        HashMap<Integer, Pair<String, String>> matchedSequences = new HashMap<>();
        Triplet<Integer, String, String> matchingResult;
        // Remove all gap symbols in the nucleotide sequence (this would rise errors in the sequence
        // alignment process).
        nucleotideSequence = removeGapsFromSequence(nucleotideSequence);
        // Check if nucleotide sequence is fully coding sequence, i.e. divisible by three.
        int residual = nucleotideSequence.length() % 3;
        int length = nucleotideSequence.length();
        // Based on the residual assign frame shift and split nucleotide sequence.
        String translatedNucSequence;
        ArrayList<String> splitNucSequence;
        if (residual == 0) {
            // Truncate no characters.
            splitNucSequence = new ArrayList<>();
            for (String s : Splitter.fixedLength(3).split(nucleotideSequence)) {
                splitNucSequence.add(s);
            }
            translatedNucSequence = translateNucSequence(splitNucSequence, true, false);
            matchingResult = globallyMatchNucleotideToAminoAcidSequence(aminoAcidSequence, translatedNucSequence,
                    splitNucSequence, "", "");
            matchedSequences
                    .put(matchingResult.getValue0(), new Pair<>(matchingResult.getValue1(), matchingResult.getValue2()));
        } else if (residual == 1) {
            // Truncate character at start.
            splitNucSequence = new ArrayList<>();
            for (String s : Splitter.fixedLength(3).split(nucleotideSequence.substring(1))) {
                splitNucSequence.add(s);
            }
            translatedNucSequence = translateNucSequence(splitNucSequence, true, false);
            matchingResult = globallyMatchNucleotideToAminoAcidSequence(aminoAcidSequence, translatedNucSequence,
                    splitNucSequence, nucleotideSequence.substring(0, 1), "");
            matchedSequences
                    .put(matchingResult.getValue0(), new Pair<>(matchingResult.getValue1(), matchingResult.getValue2()));
            // Truncate character at end.
            splitNucSequence = new ArrayList<>();
            for (String s : Splitter.fixedLength(3).split(nucleotideSequence.substring(0, length - 1))) {
                splitNucSequence.add(s);
            }
            translatedNucSequence = translateNucSequence(splitNucSequence, true, false);
            matchingResult = globallyMatchNucleotideToAminoAcidSequence(aminoAcidSequence, translatedNucSequence,
                    splitNucSequence, "", nucleotideSequence.substring(length - 1));
            matchedSequences
                    .put(matchingResult.getValue0(), new Pair<>(matchingResult.getValue1(), matchingResult.getValue2()));
        } else {
            // Truncate two characters at start.
            splitNucSequence = new ArrayList<>();
            for (String s : Splitter.fixedLength(3).split(nucleotideSequence.substring(2))) {
                splitNucSequence.add(s);
            }
            translatedNucSequence = translateNucSequence(splitNucSequence, true, false);
            matchingResult = globallyMatchNucleotideToAminoAcidSequence(aminoAcidSequence, translatedNucSequence,
                    splitNucSequence, nucleotideSequence.substring(0, 2), "");
            matchedSequences
                    .put(matchingResult.getValue0(), new Pair<>(matchingResult.getValue1(), matchingResult.getValue2()));
            // Truncate two characters at end.
            splitNucSequence = new ArrayList<>();
            for (String s : Splitter.fixedLength(3).split(nucleotideSequence.substring(0, length - 2))) {
                splitNucSequence.add(s);
            }
            translatedNucSequence = translateNucSequence(splitNucSequence, true, false);
            matchingResult = globallyMatchNucleotideToAminoAcidSequence(aminoAcidSequence, translatedNucSequence,
                    splitNucSequence, "", nucleotideSequence.substring(length - 2));
            matchedSequences
                    .put(matchingResult.getValue0(), new Pair<>(matchingResult.getValue1(), matchingResult.getValue2()));
            // Truncate one character at start and end.
            splitNucSequence = new ArrayList<>();
            for (String s : Splitter.fixedLength(3).split(nucleotideSequence.substring(1, length - 1))) {
                splitNucSequence.add(s);
            }
            translatedNucSequence = translateNucSequence(splitNucSequence, true, false);
            matchingResult = globallyMatchNucleotideToAminoAcidSequence(aminoAcidSequence, translatedNucSequence,
                    splitNucSequence, nucleotideSequence.substring(0, 1),
                    nucleotideSequence.substring(length - 1));
            matchedSequences
                    .put(matchingResult.getValue0(), new Pair<>(matchingResult.getValue1(), matchingResult.getValue2()));
        }
        // The sequences with the highest alignment score are returned as best matching.
        int maxScore = Integer.MIN_VALUE;
        Pair<String, String> bestMatching = new Pair<>("", "");
        for (Map.Entry<Integer, Pair<String, String>> entry : matchedSequences.entrySet()) {
            if (entry.getKey() > maxScore) {
                bestMatching = new Pair<>(entry.getValue().getValue0(), entry.getValue().getValue1());
            }
            maxScore = entry.getKey();
        }
        return bestMatching;
    }

    /**
     * Matches a nucleotide sequence with an amino-acid sequence based on a global alignment of the amino acid sequence
     * with the translated nucleotide sequence.
     * <p>
     * The global gap-affine alignment computation uses the BLOSUM80 scoring matrix, can handle any amino-acid
     * symbolized by an 'X' and uses an gap-opening and extension penalty of 6 and 4, respectively.
     * <p>
     * The aligned translated nucleotide sequence is mapped back to a nucleotide sequence and possible truncated
     * nucleotides are re-appended. The aligned amino-acid sequence is represented as three-letter code amino-acids
     * separated by a separator symbol. The aligned nucleotide sequence is represented by blocks of length three, i.e.
     * codons, separated by a separator symbol. Instead of the common alignment gap symbol '-', the symbol '_' is used
     * as '-' is internally reserved for deletions wrt. nucleotide sequences.
     *
     * @param aaSeq                   {@link String} representing a non-nucleotide-derived amino acid sequence.
     * @param transAASeq              {@link String} representing a amino acid sequence that originates from translating a
     *                                nucleotide sequence.
     * @param splitNucleotideSequence {@link String} the original nucleotide sequence split into chunks of length three, i.e.
     *                                codons.
     * @param truncatedStart          {@link String} nucleotides that were truncated from the start of the nucleotide sequence
     *                                due to a frameshift.
     * @param truncatedEnd            {@link String} nucleotides that were truncated from the end of the nucleotide sequence due
     *                                to a frameshift.
     * @return {@link Triplet} yielding the alignment score in the first field, the matched amino-acid sequence in the
     * second field and the matched nucleotide sequence in the third field.
     */
    private static Triplet<Integer, String, String> globallyMatchNucleotideToAminoAcidSequence(String aaSeq,
                                                                                               String transAASeq,
                                                                                               ArrayList<String> splitNucleotideSequence,
                                                                                               String truncatedStart,
                                                                                               String truncatedEnd) {
    /*
    Notes on the alignment matrices:
    - The y-axis will yield the proteins amino acid sequence.
    - The x-axis will yield the translated reference gene amino acid sequence.
    - Indels are wrt. the proteins amino acid sequence; thus insertions are traced back by walking vertically and
    deletions are traced back by walking horizontally in the matrix.
     */
        int[][] alignmentMatrix = new int[aaSeq.length() + 1][transAASeq.length() + 1];
        int[][] matchScoreMatrix = new int[aaSeq.length() + 1][transAASeq.length() + 1];
        int[][] insertionScoreMatrix = new int[aaSeq.length() + 1][transAASeq.length() + 1];
        int[][] deletionScoreMatrix = new int[aaSeq.length() + 1][transAASeq.length() + 1];
        char[][] tracebackMatrix = new char[aaSeq.length() + 1][transAASeq.length() + 1];
        char[] aaSeqArray = aaSeq.toCharArray();
        char[] transAASeqArray = transAASeq.toCharArray();
        int alignmentScore;
        StringBuilder matchedNucleotideSequenceBuilder = new StringBuilder();
        StringBuilder matchedAminoAcidSequenceBuilder = new StringBuilder();
    /*
    (1) Compute global sequence alignment.
     */
        alignmentMatrix[0][0] = 0;
        matchScoreMatrix[0][0] = 0;
        insertionScoreMatrix[0][0] = 0;
        deletionScoreMatrix[0][0] = 0;
        for (int i = 0; i < aaSeq.length() + 1; i++) {
            alignmentMatrix[i][0] = -GAP_OPEN_PENALTY - (i - 1) * GAP_EXTENSION_PENALTY;
            matchScoreMatrix[i][0] = -GAP_OPEN_PENALTY - (i - 1) * GAP_EXTENSION_PENALTY;
            insertionScoreMatrix[i][0] = -GAP_OPEN_PENALTY - (i - 1) * GAP_EXTENSION_PENALTY;
            deletionScoreMatrix[i][0] = -GAP_OPEN_PENALTY - (i - 1) * GAP_EXTENSION_PENALTY;
            tracebackMatrix[i][0] = 'I';
        }
        for (int j = 0; j < transAASeq.length() + 1; j++) {
            alignmentMatrix[0][j] = -GAP_OPEN_PENALTY - (j - 1) * GAP_EXTENSION_PENALTY;
            matchScoreMatrix[0][j] = -GAP_OPEN_PENALTY - (j - 1) * GAP_EXTENSION_PENALTY;
            insertionScoreMatrix[0][j] = -GAP_OPEN_PENALTY - (j - 1) * GAP_EXTENSION_PENALTY;
            deletionScoreMatrix[0][j] = -GAP_OPEN_PENALTY - (j - 1) * GAP_EXTENSION_PENALTY;
            tracebackMatrix[0][j] = 'D';
        }
        double max;
        for (int i = 1; i < aaSeq.length() + 1; i++) {
            for (int j = 1; j < transAASeq.length() + 1; j++) {
                matchScoreMatrix[i][j] =
                        alignmentMatrix[i - 1][j - 1]
                                + BLOSUM80[AA1_BLOSUM80_INDEX.get(aaSeqArray[i - 1])][AA1_BLOSUM80_INDEX.get(transAASeqArray[j - 1])];
                insertionScoreMatrix[i][j] =
                        Math.max(
                                alignmentMatrix[i - 1][j] - GAP_OPEN_PENALTY,
                                insertionScoreMatrix[i - 1][j] - GAP_EXTENSION_PENALTY
                        );
                deletionScoreMatrix[i][j] =
                        Math.max(
                                alignmentMatrix[i][j - 1] - GAP_OPEN_PENALTY,
                                deletionScoreMatrix[i][j - 1] - GAP_EXTENSION_PENALTY
                        );
                max = Double.NEGATIVE_INFINITY;
                if (insertionScoreMatrix[i][j] >= max) {
                    max = insertionScoreMatrix[i][j];
                    alignmentMatrix[i][j] = (int) max;
                    tracebackMatrix[i][j] = 'I';
                }
                if (deletionScoreMatrix[i][j] >= max) {
                    max = deletionScoreMatrix[i][j];
                    alignmentMatrix[i][j] = (int) max;
                    tracebackMatrix[i][j] = 'D';
                }
                if (matchScoreMatrix[i][j] >= max) {
                    max = matchScoreMatrix[i][j];
                    alignmentMatrix[i][j] = (int) max;
                    tracebackMatrix[i][j] = 'M';
                }
            }
        }
        alignmentScore = alignmentMatrix[aaSeq.length()][transAASeq.length()];
    /*
    (2) Compute traceback path from global alignment.
    */
        LinkedList<Character> tracebackPath = new LinkedList<>();
        boolean constructTracebackPath = true;
        int i = aaSeq.length();
        int j = transAASeq.length();
        char tracebackDirection = tracebackMatrix[i][j];
        while (constructTracebackPath) {
            if (tracebackDirection == 'M') {
                tracebackPath.add('M');
                i = i - 1;
                j = j - 1;
            } else if (tracebackDirection == 'D') {
                tracebackPath.add('D');
                j = j - 1;
            } else if (tracebackDirection == 'I') {
                tracebackPath.add('I');
                i = i - 1;
            }
            if (i == 0 && j == 0) {
                constructTracebackPath = false;
            } else {
                tracebackDirection = tracebackMatrix[i][j];
            }
        }
        Collections.reverse(tracebackPath);
    /*
    (3) Characters truncated at the start are re-appended to nucleotide sequence.
     */
        if (truncatedStart.length() != 0) {
            matchedNucleotideSequenceBuilder.append(truncatedStart).append(SEPARATOR_CHAR);
            matchedAminoAcidSequenceBuilder.append(String.valueOf(TRUNCATED_CHAR).repeat(truncatedStart.length()))
                    .append(SEPARATOR_CHAR);
        }
    /*
    (4) Iterate over traceback path and construct aligned nucleotide and amino-acid sequence.
    */
        int splitNucleotideSequenceIndex = 0;
        int proteinAminoAcidSequenceIndex = 0;
        String currentCodon;
        for (Character character : tracebackPath) {
            currentCodon = splitNucleotideSequence.get(splitNucleotideSequenceIndex);
            // Start walking along the traceback path.
            tracebackDirection = character;
            if (tracebackDirection == 'M') {
                // CASE: Match or mismatch of nucleotide and amino-acid sequence.
                // Add content to nucleotide sequence.
                matchedNucleotideSequenceBuilder.append(currentCodon).append(SEPARATOR_CHAR);
                splitNucleotideSequenceIndex += 1;
                // Add content to protein sequence.
                matchedAminoAcidSequenceBuilder.append(AA1TO3.get(String.valueOf(aaSeqArray[proteinAminoAcidSequenceIndex])))
                        .append(SEPARATOR_CHAR);
                proteinAminoAcidSequenceIndex += 1;
            } else if (tracebackDirection == 'D') {
                // CASE: Deletion wrt. to amino-acid sequence (of protein).
                // Add content to nucleotide sequence.
                matchedNucleotideSequenceBuilder.append(currentCodon).append(SEPARATOR_CHAR);
                splitNucleotideSequenceIndex += 1;
                // Add gap to protein sequence.
                matchedAminoAcidSequenceBuilder.append(String.valueOf(NO_MATCH_CHAR).repeat(3)).append(SEPARATOR_CHAR);
            } else if (tracebackDirection == 'I') {
                // CASE: Insertion wrt. to amino-acid sequence (of protein).
                // Add gap to nucleotide sequence.
                matchedNucleotideSequenceBuilder.append(String.valueOf(NO_MATCH_CHAR).repeat(3)).append(SEPARATOR_CHAR);
                // Add content to protein sequence.
                matchedAminoAcidSequenceBuilder.append(AA1TO3.get(String.valueOf(aaSeqArray[proteinAminoAcidSequenceIndex])))
                        .append(SEPARATOR_CHAR);
                proteinAminoAcidSequenceIndex += 1;
            }
        }
    /*
    (5) Characters truncated at the end are re-appended.
     */
        if (truncatedEnd.length() != 0) {
            matchedNucleotideSequenceBuilder.append(truncatedEnd).append(SEPARATOR_CHAR);
            matchedAminoAcidSequenceBuilder.append(String.valueOf(TRUNCATED_CHAR).repeat(truncatedEnd.length()))
                    .append(SEPARATOR_CHAR);
        }
    /*
    (6) Insert results into Triplet.
     */
        matchedNucleotideSequenceBuilder.setLength(matchedNucleotideSequenceBuilder.length() - 1);
        String matchedNucleotideSequence = matchedNucleotideSequenceBuilder.toString();
        matchedAminoAcidSequenceBuilder.setLength(matchedAminoAcidSequenceBuilder.length() - 1);
        String matchedAminoAcidSequence = matchedAminoAcidSequenceBuilder.toString();
        return new Triplet<>(alignmentScore, matchedAminoAcidSequence, matchedNucleotideSequence);
    }

    /**
     * Allocates sample variants data and nucleotide reference data to a specified reference protein structures chain.
     *
     * @param variantsTable         {@link VariantContentTable} instance yielding all sample variants information.
     * @param referenceFeatureEntry {@link ReferenceFeatureEntry} internal representation of the reference feature.
     * @param proteinStructure      {@link Structure} yielding the protein structure that shall be used as the with the
     *                              reference feature associated protein structure.
     * @param chainId               {@link String} chain identifier to access protein data from {@link Structure} objects and naming.
     * @param chainAsymId           {@link String} internal chain identifier to access protein data from {@link Structure} objects.
     * @param isSense               {@link Boolean} indicating whether the feature to analyze is located on the sense strand.
     * @return {@link JSONObject} yielding sample and reference nucleotide information allocated to protein positions.
     * @throws MusialBioException If any invoked method call (i.e. codon translation or base reversal) throws this
     *                            exception.
     */
    @SuppressWarnings("unchecked")
    public static JSONObject superimposeDataToStructure(VariantContentTable variantsTable,
                                                        ReferenceFeatureEntry referenceFeatureEntry,
                                                        Structure proteinStructure, String chainId, String chainAsymId,
                                                        boolean isSense)
            throws MusialBioException {
        JSONObject superpositionData = new JSONObject(); // Final result to be returned.
        ConcurrentSkipListMap<String, JSONObject> superpositionSortedMap =
                new ConcurrentSkipListMap<>(new VariantPositionComparator()); // Superpositions; may contain
        /*
        (1) Access reference feature and amino-acid sequence and match/align both sequences.
         */
        String referenceSequence = getSequenceFromTable(variantsTable, referenceFeatureEntry, "Reference",
                false, false);
        String referenceFeatureName = referenceFeatureEntry.name;
        if (!isSense) {
            referenceSequence = Bio.getReverseComplement(referenceSequence);
        }
        String aminoAcidSequence = getSequencesFromPdb(proteinStructure).get(chainId);
        Pair<String, String> matchedSequences = Bio.getNucleotideToAminoAcidSequenceMatching(aminoAcidSequence,
                referenceSequence);
        Iterator<String> matchedAminoAcidSequenceIterator =
                Arrays.stream(matchedSequences.getValue0().split(String.valueOf(SEPARATOR_CHAR))).iterator();
        Iterator<String> matchedReferenceSequenceIterator =
                Arrays.stream(matchedSequences.getValue1().split(String.valueOf(SEPARATOR_CHAR))).iterator();
        /*
        (2) Insert reference feature and amino-acid information into superposition data object.
        */
        superpositionData.put("ReferenceProteinLength", aminoAcidSequence.length());
        // Get first position of protein structure from .pdb file.
        int proteinPosition = Bio.getFirstResidueNumberFromProteinStructure(proteinStructure, chainAsymId);
        // Get first reference feature position.
        int referenceFeatureStartPosition;
        if (isSense) {
            referenceFeatureStartPosition = referenceFeatureEntry.locationStart;
        } else {
            // If feature is on anti-sense strand we multiply by -1 to indicate strandedness.
            referenceFeatureStartPosition = -referenceFeatureEntry.locationEnd;
        }
        int referenceFeaturePosition = referenceFeatureStartPosition;
        // Index is used as common coordinate system, i.e. the superimposed position.
        int index = 1;
        while (matchedAminoAcidSequenceIterator.hasNext() && matchedReferenceSequenceIterator.hasNext()) {
            // Get next content of matched reference sequences.
            String referenceAminoAcidContent = matchedAminoAcidSequenceIterator.next();
            String referenceNucleotideContent = matchedReferenceSequenceIterator.next();
            // For each index a JSON object is generated that assigns relevant information to that superimposed position.
            JSONObject indexData = new JSONObject();
            // Add reference feature information.
            ArrayList<Integer> referenceFeaturePositions = new ArrayList<>();
            if (referenceNucleotideContent.contains(String.valueOf(NO_MATCH_CHAR))) {
                indexData.put("ReferenceFeatureContent", NONE);
                indexData.put("ReferenceFeatureTranslatedContent", NONE);
                referenceFeaturePositions.add(-1);
                indexData.put("ReferenceFeaturePositions", referenceFeaturePositions);
            } else {
                if (referenceNucleotideContent.length() == 3) {
                    // If the reference content does correspond to a full codon.
                    indexData.put("ReferenceFeatureContent", referenceNucleotideContent);
                    indexData
                            .put("ReferenceFeatureTranslatedContent", translateCodon(referenceNucleotideContent, true
                                    , true, false));
                    for (int i = 0; i < referenceNucleotideContent.toCharArray().length; i++) {
                        referenceFeaturePositions.add(referenceFeaturePosition);
                        referenceFeaturePosition += 1;
                    }
                    indexData.put("ReferenceFeaturePositions", referenceFeaturePositions);
                } else {
                    indexData.put("ReferenceFeatureContent", referenceNucleotideContent);
                    indexData.put("ReferenceFeatureTranslatedContent", NONE);
                    for (int i = 0; i < referenceNucleotideContent.toCharArray().length; i++) {
                        referenceFeaturePositions.add(referenceFeaturePosition);
                        referenceFeaturePosition += 1;
                    }
                    indexData.put("ReferenceFeaturePositions", referenceFeaturePositions);
                }
            }
            // Add protein amino-acid information.
            if (referenceAminoAcidContent.contains(String.valueOf(NO_MATCH_CHAR))) {
                indexData.put("ProteinContent", NONE);
                indexData.put("ProteinPosition", -1);
            } else {
                indexData.put("ProteinContent", referenceAminoAcidContent);
                indexData.put("ProteinPosition", proteinPosition);
                proteinPosition += 1;
            }
            // Prepare containers for variant data.
            indexData.put("PerGenotypeNucleotideContent", new JSONObject());
            indexData.put("PerGenotypePositions", new JSONObject());
            indexData.put("PerProteoformAminoAcidContent", new JSONObject());
            superpositionSortedMap.put(String.valueOf(index), indexData);
            index += 1;
        }

        /*
        (3) Extract sample information; i.e. insertions and sequences for downstream processing and to infer
        genotypes and proteoforms.
         */
        String sampleSequence;
        ArrayList<String> splitSampleSequence;
        String translatedSampleSequence;
        String genotypeSequence;
        String genotypeName;
        String proteoformSequence;
        String proteoformName;
        HashMap<String, ArrayList<String>> genotypesSamplesBySequence = new HashMap<>();
        ArrayList<String> genotypeSampleNames;
        JSONObject genotypesJSONObject = new JSONObject();
        JSONObject genotypeJSONObject;
        HashMap<String, String> genotypeSequences = new HashMap<>();
        HashMap<String, ArrayList<String>> proteoformGenotypesBySequence = new HashMap<>();
        ArrayList<String> proteoformGenotypeNames;
        JSONObject proteoformsJSONObject = new JSONObject();
        JSONObject proteoformJSONObject;
        HashMap<String, String> proteoformSequences = new HashMap<>();
        HashMap<String, Iterator<Pair<String, Character>>> genotypePositionAnnotatedSequences = new HashMap<>();
        HashMap<String, String> proteoformRepresentatives = new HashMap<>();
        int genotypeCounter = 1;
        int proteoformCounter = 1;
        for (String sampleName : variantsTable.getSampleNames(referenceFeatureName)) {
            sampleSequence = getSequenceFromTable(variantsTable, referenceFeatureEntry, sampleName, false, true);
            if (!isSense) {
                sampleSequence = getReverseComplement(sampleSequence);
            }
            if (genotypesSamplesBySequence.containsKey(sampleSequence)) {
                genotypesSamplesBySequence.get(sampleSequence).add(sampleName);
            } else {
                genotypeName = "Genotype" + genotypeCounter;
                genotypeSampleNames = new ArrayList<>();
                genotypeSampleNames.add(sampleName);
                genotypesSamplesBySequence.put(sampleSequence, genotypeSampleNames);
                genotypePositionAnnotatedSequences.put(genotypeName, getPositionAnnotatedSequenceIteratorFromTable(variantsTable, referenceFeatureEntry, sampleName));
                genotypeSequences.put(genotypeName, sampleSequence);
                genotypeCounter += 1;
                splitSampleSequence = new ArrayList<>();
                for (String s : Splitter.fixedLength(3).split(sampleSequence)) {
                    splitSampleSequence.add(s);
                }
                translatedSampleSequence = translateNucSequence(splitSampleSequence, true, true);
                if (proteoformGenotypesBySequence.containsKey(translatedSampleSequence)) {
                    proteoformGenotypesBySequence.get(translatedSampleSequence).add(genotypeName);
                } else {
                    proteoformName = "Proteoform" + proteoformCounter;
                    proteoformGenotypeNames = new ArrayList<>();
                    proteoformGenotypeNames.add(genotypeName);
                    proteoformGenotypesBySequence.put(translatedSampleSequence, proteoformGenotypeNames);
                    proteoformSequences.put(proteoformName, translatedSampleSequence);
                    proteoformRepresentatives.put(genotypeName, proteoformName);
                    proteoformCounter += 1;
                }
            }
        }
        for (Map.Entry<String, String> entry : genotypeSequences.entrySet()) {
            genotypeName = entry.getKey();
            genotypeSequence = entry.getValue();
            genotypeJSONObject = new JSONObject();
            genotypeSampleNames = genotypesSamplesBySequence.get(genotypeSequence);
            genotypeJSONObject.put("Samples", genotypeSampleNames);
            genotypeJSONObject.put("NucleotideSequence", genotypeSequence);
            genotypesJSONObject.put(genotypeName, genotypeJSONObject);
        }
        superpositionData.put("Genotypes", genotypesJSONObject);
        for (Map.Entry<String, String> entry : proteoformSequences.entrySet()) {
            proteoformName = entry.getKey();
            proteoformSequence = entry.getValue();
            proteoformJSONObject = new JSONObject();
            proteoformGenotypeNames = proteoformGenotypesBySequence.get(proteoformSequence);
            proteoformJSONObject.put("Genotypes", proteoformGenotypeNames);
            proteoformJSONObject.put("AminoAcidSequence", proteoformSequence);
            proteoformsJSONObject.put(proteoformName, proteoformJSONObject);
        }
        superpositionData.put("Proteoforms", proteoformsJSONObject);

        /*
        (4) For each genotype, superimpose variant content with reference data.
         */
        Iterator<Pair<String, Character>> positionAnnotatedSequence;
        Pair<String, Character> positionAnnotatedContent;
        String position;
        char content;
        int imbalance;
        int insertionCount;
        StringBuilder codonsBuilder = new StringBuilder();
        ArrayList<String> codons;
        String codon;
        String translatedCodon;
        ArrayList<String> codonsPositions = new ArrayList<>();
        JSONObject indexData;
        for (Map.Entry<String, Iterator<Pair<String, Character>>> entry : genotypePositionAnnotatedSequences.entrySet()) {
            genotypeName = entry.getKey();
            positionAnnotatedSequence = entry.getValue();
            imbalance = 0;
            insertionCount = 1;
            index = 0;
            while (positionAnnotatedSequence.hasNext()) {
                positionAnnotatedContent = positionAnnotatedSequence.next();
                position = positionAnnotatedContent.getValue0();
                content = positionAnnotatedContent.getValue1();
                if (position.contains("+")) {
                    imbalance += 1;
                } else {
                    insertionCount = 1;
                }
                if (content == VariantContent.DELETION) {
                    imbalance -= 1;
                } else {
                    codonsPositions.add(position);
                    codonsBuilder.append(content);
                }
                if (((codonsBuilder.length() % 3 == 0) && (Math.abs(imbalance) % 3 == 0)) || !positionAnnotatedSequence.hasNext()) {
                    int deletions = 0;
                    boolean ambiguousDeletions = false;
                    if (codonsBuilder.length() == 0 && imbalance == -3) {
                        // Case in which only one in-frame deletion was observed.
                        deletions = 1;
                    } else {
                        codons = new ArrayList<>();
                        Splitter.fixedLength(3).split(codonsBuilder.toString()).iterator().forEachRemaining(codons::add);
                        int codonIndex = 0;
                        for (String c : codons) {
                            codon = c;
                            translatedCodon = translateCodon(codon, true, true, true);
                            if (imbalance == 0) {
                                // If imbalance is zero, indel events have neutralized each other or no such event was observed.
                                index += 1;
                                indexData = superpositionSortedMap.get(String.valueOf(index));
                                ((JSONObject) indexData.get("PerGenotypeNucleotideContent")).put(genotypeName,
                                        codon);
                                JSONArray pos = new JSONArray();
                                for (int i = 0; i < 3; i++) {
                                    if ( !codonsPositions.isEmpty() ){
                                        pos.add(codonsPositions.remove(0));
                                    }
                                }
                                ((JSONObject) indexData.get("PerGenotypePositions")).put(genotypeName, pos);
                                if (proteoformRepresentatives.containsKey(genotypeName)) {
                                    ((JSONObject) indexData.get("PerProteoformAminoAcidContent")).put(proteoformRepresentatives.get(genotypeName),
                                            translatedCodon);
                                }
                            } else if (imbalance > 0) {
                                // If imbalance is greater than zero, at least one insertion was observed.
                                if (codons.size() == 1) {
                                    // If exactly one codon was parsed, an in-frame insertion was observed.
                                    if (!superpositionSortedMap.containsKey(index + "+" + insertionCount)) {
                                        indexData = new JSONObject();
                                        indexData.put("ReferenceFeatureContent", NONE);
                                        indexData.put("ReferenceFeatureTranslatedContent", NONE);
                                        JSONArray pos = new JSONArray();
                                        pos.add(-1);
                                        pos.add(-1);
                                        pos.add(-1);
                                        indexData.put("ReferenceFeaturePositions", pos);
                                        indexData.put("ProteinContent", NONE);
                                        indexData.put("ProteinPosition", -1);
                                        indexData.put("PerGenotypeNucleotideContent", new JSONObject());
                                        indexData.put("PerGenotypePositions", new JSONObject());
                                        indexData.put("PerProteoformAminoAcidContent", new JSONObject());
                                        ((JSONObject) indexData.get("PerGenotypeNucleotideContent")).put(genotypeName,
                                                codon);
                                        pos = new JSONArray();
                                        for (int i = 0; i < 3; i++) {
                                            if ( !codonsPositions.isEmpty() ){
                                                pos.add(codonsPositions.remove(0));
                                            }
                                        }
                                        ((JSONObject) indexData.get("PerGenotypePositions")).put(genotypeName, pos);
                                        if (proteoformRepresentatives.containsKey(genotypeName)) {
                                            ((JSONObject) indexData.get("PerProteoformAminoAcidContent")).put(proteoformRepresentatives.get(genotypeName),
                                                    translatedCodon);
                                        }
                                        superpositionSortedMap.put(index + "+" + insertionCount, indexData);
                                    } else {
                                        indexData = superpositionSortedMap.get(index + "+" + insertionCount);
                                        ((JSONObject) indexData.get("PerGenotypeNucleotideContent")).put(genotypeName,
                                                codon);
                                        JSONArray pos = new JSONArray();
                                        for (int i = 0; i < 3; i++) {
                                            if ( !codonsPositions.isEmpty() ){
                                                pos.add(codonsPositions.remove(0));
                                            }
                                        }
                                        ((JSONObject) indexData.get("PerGenotypePositions")).put(genotypeName, pos);
                                        if (proteoformRepresentatives.containsKey(genotypeName)) {
                                            ((JSONObject) indexData.get("PerProteoformAminoAcidContent")).put(proteoformRepresentatives.get(genotypeName),
                                                    translatedCodon);
                                        }
                                    }
                                    insertionCount += 1;
                                } else {
                                    // Otherwise, multiple events led to the observed imbalance.
                                    int insertions = (int) Math.ceil( imbalance / 3.0 );
                                    translatedCodon = translatedCodon.toLowerCase(Locale.ROOT); // Indicate ambiguity.
                                    if (codonIndex < (codons.size() - insertions)) {
                                        index += 1;
                                        indexData = superpositionSortedMap.get(String.valueOf(index));
                                        ((JSONObject) indexData.get("PerGenotypeNucleotideContent")).put(genotypeName,
                                                codon);
                                        JSONArray pos = new JSONArray();
                                        for (int i = 0; i < 3; i++) {
                                            if ( !codonsPositions.isEmpty() ){
                                                pos.add(codonsPositions.remove(0));
                                            }
                                        }
                                        ((JSONObject) indexData.get("PerGenotypePositions")).put(genotypeName, pos);
                                        if (proteoformRepresentatives.containsKey(genotypeName)) {
                                            ((JSONObject) indexData.get("PerProteoformAminoAcidContent")).put(proteoformRepresentatives.get(genotypeName),
                                                    translatedCodon);
                                        }
                                    } else {
                                        if (!superpositionSortedMap.containsKey(index + "+" + insertionCount)) {
                                            indexData = new JSONObject();
                                            indexData.put("ReferenceFeatureContent", NONE);
                                            indexData.put("ReferenceFeatureTranslatedContent", NONE);
                                            JSONArray pos = new JSONArray();
                                            pos.add(-1);
                                            pos.add(-1);
                                            pos.add(-1);
                                            indexData.put("ReferenceFeaturePositions", pos);
                                            indexData.put("ProteinContent", NONE);
                                            indexData.put("ProteinPosition", -1);
                                            indexData.put("PerGenotypeNucleotideContent", new JSONObject());
                                            indexData.put("PerGenotypePositions", new JSONObject());
                                            indexData.put("PerProteoformAminoAcidContent", new JSONObject());
                                            ((JSONObject) indexData.get("PerGenotypeNucleotideContent")).put(genotypeName,
                                                    codon);
                                            pos = new JSONArray();
                                            for (int i = 0; i < 3; i++) {
                                                if ( !codonsPositions.isEmpty() ){
                                                    pos.add(codonsPositions.remove(0));
                                                }
                                            }
                                            ((JSONObject) indexData.get("PerGenotypePositions")).put(genotypeName, pos);
                                            if (proteoformRepresentatives.containsKey(genotypeName)) {
                                                ((JSONObject) indexData.get("PerProteoformAminoAcidContent")).put(proteoformRepresentatives.get(genotypeName),
                                                        translatedCodon);
                                            }
                                            superpositionSortedMap.put(index + "+" + insertionCount, indexData);
                                        } else {
                                            indexData = superpositionSortedMap.get(index + "+" + insertionCount);
                                            ((JSONObject) indexData.get("PerGenotypeNucleotideContent")).put(genotypeName,
                                                    codon);
                                            JSONArray pos = new JSONArray();
                                            for (int i = 0; i < 3; i++) {
                                                if ( !codonsPositions.isEmpty() ){
                                                    pos.add(codonsPositions.remove(0));
                                                }
                                            }
                                            ((JSONObject) indexData.get("PerGenotypePositions")).put(genotypeName, pos);
                                            if (proteoformRepresentatives.containsKey(genotypeName)) {
                                                ((JSONObject) indexData.get("PerProteoformAminoAcidContent")).put(proteoformRepresentatives.get(genotypeName),
                                                        translatedCodon);
                                            }
                                        }
                                        insertionCount += 1;
                                    }
                                }
                            } else {
                                // If imbalance is below zero, at least one deletion was observed.
                                deletions = Math.abs(imbalance) / 3;
                                ambiguousDeletions = true;
                                translatedCodon = translatedCodon.toLowerCase(Locale.ROOT); // Indicate ambiguity.
                                index += 1;
                                indexData = superpositionSortedMap.get(String.valueOf(index));
                                ((JSONObject) indexData.get("PerGenotypeNucleotideContent")).put(genotypeName,
                                        codon);
                                JSONArray pos = new JSONArray();
                                for (int i = 0; i < 3; i++) {
                                    if ( !codonsPositions.isEmpty() ){
                                        pos.add(codonsPositions.remove(0));
                                    }
                                }
                                ((JSONObject) indexData.get("PerGenotypePositions")).put(genotypeName, pos);
                                if (proteoformRepresentatives.containsKey(genotypeName)) {
                                    ((JSONObject) indexData.get("PerProteoformAminoAcidContent")).put(proteoformRepresentatives.get(genotypeName),
                                            translatedCodon);
                                }
                            }
                            codonIndex += 1;
                        }
                    }
                    while (deletions != 0) {
                        // For each observed deletion insert one dummy deletion.
                        index += 1;
                        translatedCodon = DELETION_AA3;
                        indexData = superpositionSortedMap.get(String.valueOf(index));
                        ((JSONObject) indexData.get("PerGenotypeNucleotideContent")).put(genotypeName,
                                NONE);
                        JSONArray pos = new JSONArray();
                        pos.add(-1);
                        pos.add(-1);
                        pos.add(-1);
                        ((JSONObject) indexData.get("PerGenotypePositions")).put(genotypeName, pos);
                        if (proteoformRepresentatives.containsKey(genotypeName)) {
                            if (ambiguousDeletions) {
                                translatedCodon = translatedCodon.toLowerCase(Locale.ROOT);
                            }
                            ((JSONObject) indexData.get("PerProteoformAminoAcidContent")).put(proteoformRepresentatives.get(genotypeName),
                                    translatedCodon);
                        }
                        deletions -= 1;
                    }
                    codonsBuilder = new StringBuilder();
                    codonsPositions = new ArrayList<>();
                    imbalance = 0;
                }
            }
        }

        // (5) Add superimposed data to JSON object.
        int shift = 0;
        for (
                Map.Entry<String, JSONObject> entry : superpositionSortedMap.entrySet()) {
            position = entry.getKey();
            if (position.contains("+")) {
                shift += 1;
                position = String.valueOf(Integer.parseInt(position.substring(0, position.indexOf("+"))) + shift);
            } else {
                position = String.valueOf(Integer.parseInt(position) + shift);
            }
            indexData = entry.getValue();
            for (String pfName : proteoformSequences.keySet()) {
                if ( !((JSONObject) indexData.get("PerProteoformAminoAcidContent")).containsKey(pfName) ) {
                    ((JSONObject) indexData.get("PerProteoformAminoAcidContent")).put(pfName, NONE);
                }
            }
            for (String gtName : genotypeSequences.keySet()) {
                if ( !((JSONObject) indexData.get("PerGenotypeNucleotideContent")).containsKey(gtName) ) {
                    ((JSONObject) indexData.get("PerGenotypeNucleotideContent")).put(gtName, NONE);
                }
                if ( !((JSONObject) indexData.get("PerGenotypePositions")).containsKey(gtName) ) {
                    JSONArray pos = new JSONArray();
                    pos.add(-1);
                    pos.add(-1);
                    pos.add(-1);
                    ((JSONObject) indexData.get("PerGenotypePositions")).put(gtName, pos);
                }
            }
            superpositionData.put(position, indexData);
        }
        return superpositionData;
    }

    /**
     * Translates a nucleotide sequence, split into codons of length three, into a single-letter amino acid sequence.
     *
     * @param splitNucSequence   {@link ArrayList<String>} representing a nucleotide sequence that was split into codons
     *                           of length 3.
     * @param includeTermination {@link Boolean} whether a representation for translation termination (in the case of
     *                           three-letter code 'TER' and in the case of one-letter code 'Z') shall
     *                           be returned or not, i.e. an empty String is returned instead.
     * @param includeIncomplete  {@link Boolean} whether an incomplete amino-acid should be added to the end if the
     *                           sequence contains an incomplete codon at the end.
     * @return {@link String} representing the translated nucleotide sequence.
     * @throws MusialBioException If any codon with a length different than three is detected.
     */
    public static String translateNucSequence(ArrayList<String> splitNucSequence, boolean includeTermination,
                                              boolean includeIncomplete) throws MusialBioException {
        StringBuilder translatedNucSequenceBuilder = new StringBuilder();
        for (String s : splitNucSequence) {
            if (s.length() != 3) {
                if (!includeIncomplete) {
                    throw new MusialBioException("Failed to translate nucleotide sequence containing codon of length unequal 3.");
                }
            } else {
                translatedNucSequenceBuilder.append(translateCodon(s, false, includeTermination, includeIncomplete));
            }
        }
        return translatedNucSequenceBuilder.toString();
    }

    /**
     * Returns the reverse complement of the passed nucleotide sequence.
     *
     * @param sequence {@link String} representing a nucleotide sequence.
     * @return {@link String} representing the reverse complement of the passed sequence.
     * @throws MusialBioException If any base occurs in the sequence that can not be reversed, i.e. for which no
     *                            complement base is known.
     */
    public static String getReverseComplement(String sequence) throws MusialBioException {
        StringBuilder reverseComplementBuilder = new StringBuilder();
        char[] sequenceCharArray = sequence.toCharArray();
        for (int i = sequenceCharArray.length - 1; i >= 0; i--) {
            char sequenceAtI = sequenceCharArray[i];
            reverseComplementBuilder.append(invertBase(sequenceAtI));
        }
        return reverseComplementBuilder.toString();
    }

    /**
     * Returns the complement of a nucleotide base, i.e. switches A to T, C to G and vice versa. For non nucleotide
     * symbols, e.g. gaps (-) or any nucleotides (N), the identity is returned.
     *
     * @param base {@link Character} base to invert.
     * @return {@link Character} the inverted base.
     * @throws MusialBioException If no complement base for the specified base is known.
     */
    public static Character invertBase(char base) throws MusialBioException {
        return switch (base) {
            case VariantContent.ALT_A -> VariantContent.ALT_T;
            case VariantContent.ALT_C -> VariantContent.ALT_G;
            case VariantContent.ALT_G -> VariantContent.ALT_C;
            case VariantContent.ALT_T -> VariantContent.ALT_A;
            case VariantContent.REFERENCE -> VariantContent.REFERENCE;
            case VariantContent.NO_CALL -> VariantContent.NO_CALL;
            case VariantContent.DELETION -> VariantContent.DELETION;
            case NO_MATCH_CHAR -> NO_MATCH_CHAR;
            case VariantContent.INSERTION_DUMMY -> VariantContent.INSERTION_DUMMY;
            default -> throw new MusialBioException("Unable to invert base symbol " + base + ".");
        };
    }

    /**
     * Returns a amino-acid representation of the passed `codon`.
     *
     * @param codon              {@link String} representing the nucleotide codon to translate.
     * @param asAA3              {@link Boolean} whether the codon should be translated to the amino-acid three-letter code or not.
     * @param includeTermination {@link Boolean} whether a representation for translation termination (in the case of
     *                           three-letter code 'TER' and in the case of one-letter code 'Z') shall
     *                           be returned or not, i.e. an empty String is returned instead.
     * @param includeIncomplete  {@link Boolean} whether incomplete codons, i.e. with a length other than 3, should be
     *                           translated as incomplete amino-acids.
     * @return {@link String} representing a translated amino-acid.
     * @throws MusialBioException If the passed codon has a length other than three and `includeIncomplete` is false.
     */
    public static String translateCodon(String codon, boolean asAA3, boolean includeTermination, boolean includeIncomplete)
            throws MusialBioException {
        if (codon.length() != 3) {
            if (!includeIncomplete) {
                throw new MusialBioException("Unable to translate codon " + codon + " with length different from three.");
            } else {
                if (asAA3) {
                    return INCOMPLETE_AA3;
                } else {
                    return String.valueOf(INCOMPLETE_AA1);
                }
            }
        } else if (codon.contains(String.valueOf(VariantContent.NO_CALL))) {
            if (asAA3) {
                return UNKNOWN_AA3;
            } else {
                return String.valueOf(UNKNOWN_AA1);
            }
        } else if (CODON_MAP.get(codon).equals("") && includeTermination) {
            if (asAA3) {
                return TERMINATION_AA3;
            } else {
                return String.valueOf(TERMINATION_AA1);
            }
        } else {
            if (CODON_MAP.get(codon).equals("")) {
                return "";
            } else {
                if (asAA3) {
                    return AA1TO3.get(CODON_MAP.get(codon));
                } else {
                    return CODON_MAP.get(codon);
                }
            }
        }
    }

    /**
     * Removes all gap symbols ("-") from the passed sequence.
     *
     * @param sequence The sequence from which gaps should be removed.
     * @return {@link String}, the passed sequence without any gap symbols.
     */
    public static String removeGapsFromSequence(String sequence) {
        return sequence.replace(String.valueOf(VariantContent.DELETION), "");
    }

    /**
     * Returns the chains sequences from a {@link Structure} instance parsed from a .pdb file.
     *
     * @param structure {@link Structure} instance to extract the chains sequences from.
     * @return {@link HashMap} storing the .pdb files amino-acid sequences per chain.
     */
    public static HashMap<String, String> getSequencesFromPdb(Structure structure) {
        // Initialize results list.
        HashMap<String, String> pdbChainsSequences = new HashMap<>();
        // Parse nucleotide information from chains.
        int modelNr = 0;
        List<Chain> chains = structure.getChains(modelNr);
        for (Chain chain : chains) {
            pdbChainsSequences.put(chain.getName(), chain.getAtomSequence());
        }
        return pdbChainsSequences;
    }

    /**
     * Returns the first residue number from a {@link Structure} instances specified chain (via chainId).
     *
     * @param structure   {@link Structure} instance to access.
     * @param chainAsymId {@link String} representing the AsymId of the chain to access.
     * @return The first residue number of the specified structure and chain.
     */
    public static int getFirstResidueNumberFromProteinStructure(Structure structure, String chainAsymId) {
        // Parse nucleotide information from chains.
        Chain chain = structure.getChain(chainAsymId);
        ResidueNumber residueNumber = chain.getAtomGroups().get(0).getResidueNumber();
        return residueNumber.getSeqNum();
    }

    /**
     * Returns a {@link HashMap} mapping residue number and chain names as (concatenated) keys to the respective
     * {@link Group} instances.
     *
     * @param structure {@link Structure} instance to extract the residue {@link Group} objects from.
     * @return The mapping of residue identifiers to residues.
     */
    public static HashMap<String, Group> getResiduesFromStructure(Structure structure) {
        // Initialize results map.
        HashMap<String, Group> residueMap = new HashMap<>();
        for (Chain chain : structure.getChains()) {
            for (Group atomGroup : chain.getAtomGroups()) {
                ResidueNumber seqResNumber = atomGroup.getResidueNumber();
                residueMap.put(seqResNumber.getSeqNum() + seqResNumber.getChainName(), atomGroup);
            }
        }
        return residueMap;
    }

    /**
     * Computes the side-chains center of mass for a specified set of amino-acid residues represented as {@link Group}
     * objects.
     *
     * @param residues {@link HashMap} mapping from residue identifiers to {@link Group} instances.
     * @return The mapping of residue identifiers to residue side-chains center of mass.
     */
    public static HashMap<String, ArrayList<Double>> getSidechainCentersOfMass(HashMap<String, Group> residues) {
        // Initialize results map.
        HashMap<String, ArrayList<Double>> residueSidechainCentersOfMass = new HashMap<>();
        for (String identifier : residues.keySet()) {
            Group residue = residues.get(identifier);
            ArrayList<Double> sidechainCenterOfMassCoords = new ArrayList<>();
            // First exclude all non-sidechain atoms.
            ArrayList<Atom> sidechainAtoms = new ArrayList<>();
            for (Atom atom : residue.getAtoms()) {
                String atomName = atom.getName();
                if (!(atomName.equals("N") || atomName.equals("O") || atomName.equals("C") || atomName.equals("CA") ||
                        atomName.equals("H") || atomName.equals("HA") || atomName.equals("1HA") || atomName.equals("2HA"))) {
                    sidechainAtoms.add(atom);
                }
            }
            // If no sidechain was detected use C-alpha atom.
            if (sidechainAtoms.size() == 0) {
                for (Atom atom : residue.getAtoms()) {
                    String atomName = atom.getName();
                    if (atomName.equals("CA")) {
                        sidechainAtoms.add(atom);
                    }
                }
            }
            // Compute the side-chains center of mass:
            // (1) Compute the total mass of the sidechain atoms.
            double M = sidechainAtoms.stream().mapToDouble(a -> a.getElement().getAtomicMass()).sum();
            // (2) Compute x-coord.
            sidechainCenterOfMassCoords.add(
                    sidechainAtoms.stream().mapToDouble(a -> (a.getElement().getAtomicMass() * a.getX())).sum() / M);
            // (3) Compute y-coord.
            sidechainCenterOfMassCoords.add(
                    sidechainAtoms.stream().mapToDouble(a -> (a.getElement().getAtomicMass() * a.getY())).sum() / M);
            // (4) Compute z-coord.
            sidechainCenterOfMassCoords.add(
                    sidechainAtoms.stream().mapToDouble(a -> (a.getElement().getAtomicMass() * a.getZ())).sum() / M);
            residueSidechainCentersOfMass.put(identifier, sidechainCenterOfMassCoords);
        }
        return residueSidechainCentersOfMass;
    }

}
