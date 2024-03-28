package utility;

import com.google.common.base.Splitter;
import exceptions.MusialException;
import main.Constants;
import org.apache.commons.lang3.tuple.Triple;
import org.javatuples.Triplet;

import java.util.*;

/**
 * Implements static methods for sequence operations, such as inversion, translation, and alignment.
 *
 * @author Simon Hackl
 * @version 2.3
 * @since 2.0
 */
public final class SequenceOperations {

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
     * Hash map mapping amino-acid one to three-letter codes.
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
        put(Constants.ANY_AMINOACID_STRING, Constants.ANY_AMINOACID3_STRING);
        put(Constants.TERMINATION_AMINOACID_STRING, Constants.TERMINATION_AMINOACID3_STRING);
        put(Constants.DELETION_OR_GAP_STRING, Constants.DELETION_OR_GAP_3_STRING);
    }};
    /**
     * Hash map mapping amino-acid one letter symbols to indices used to access substitution matrix values.
     */
    public static final HashMap<Character, Integer> AA1_PAM250_INDEX = new HashMap<>() {{
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
        put('B', 20);
        put('J', 21);
        put('Z', 22);
        put(Constants.ANY_AMINOACID_CHAR, 23);
        put(Constants.TERMINATION_AMINOACID_CHAR, 24);
    }};
    /**
     * The PAM250 matrix used for alignment computation.
     * <p>
     * Source <a href="https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/PAM250">ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/PAM250</a>
     */
    private static final int[][] PAM250 = {
            {2, -2, 0, 0, -2, 0, 0, 1, -1, -1, -2, -1, -1, -3, 1, 1, 1, -6, -3, 0, 0, -1, 0, -1, -8},
            {-2, 6, 0, -1, -4, 1, -1, -3, 2, -2, -3, 3, 0, -4, 0, 0, -1, 2, -4, -2, -1, -3, 0, -1, -8},
            {0, 0, 2, 2, -4, 1, 1, 0, 2, -2, -3, 1, -2, -3, 0, 1, 0, -4, -2, -2, 2, -3, 1, -1, -8},
            {0, -1, 2, 4, -5, 2, 3, 1, 1, -2, -4, 0, -3, -6, -1, 0, 0, -7, -4, -2, 3, -3, 3, -1, -8},
            {-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3, 0, -2, -8, 0, -2, -4, -5, -5, -1, -8},
            {0, 1, 1, 2, -5, 4, 2, -1, 3, -2, -2, 1, -1, -5, 0, -1, -1, -5, -4, -2, 1, -2, 3, -1, -8},
            {0, -1, 1, 3, -5, 2, 4, 0, 1, -2, -3, 0, -2, -5, -1, 0, 0, -7, -4, -2, 3, -3, 3, -1, -8},
            {1, -3, 0, 1, -3, -1, 0, 5, -2, -3, -4, -2, -3, -5, 0, 1, 0, -7, -5, -1, 0, -4, 0, -1, -8},
            {-1, 2, 2, 1, -3, 3, 1, -2, 6, -2, -2, 0, -2, -2, 0, -1, -1, -3, 0, -2, 1, -2, 2, -1, -8},
            {-1, -2, -2, -2, -2, -2, -2, -3, -2, 5, 2, -2, 2, 1, -2, -1, 0, -5, -1, 4, -2, 3, -2, -1, -8},
            {-2, -3, -3, -4, -6, -2, -3, -4, -2, 2, 6, -3, 4, 2, -3, -3, -2, -2, -1, 2, -3, 5, -3, -1, -8},
            {-1, 3, 1, 0, -5, 1, 0, -2, 0, -2, -3, 5, 0, -5, -1, 0, 0, -3, -4, -2, 1, -3, 0, -1, -8},
            {-1, 0, -2, -3, -5, -1, -2, -3, -2, 2, 4, 0, 6, 0, -2, -2, -1, -4, -2, 2, -2, 3, -2, -1, -8},
            {-3, -4, -3, -6, -4, -5, -5, -5, -2, 1, 2, -5, 0, 9, -5, -3, -3, 0, 7, -1, -4, 2, -5, -1, -8},
            {1, 0, 0, -1, -3, 0, -1, 0, 0, -2, -3, -1, -2, -5, 6, 1, 0, -6, -5, -1, -1, -2, 0, -1, -8},
            {1, 0, 1, 0, 0, -1, 0, 1, -1, -1, -3, 0, -2, -3, 1, 2, 1, -2, -3, -1, 0, -2, 0, -1, -8},
            {1, -1, 0, 0, -2, -1, 0, 0, -1, 0, -2, 0, -1, -3, 0, 1, 3, -5, -3, 0, 0, -1, -1, -1, -8},
            {-6, 2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4, 0, -6, -2, -5, 17, 0, -6, -5, -3, -6, -1, -8},
            {-3, -4, -2, -4, 0, -4, -4, -5, 0, -1, -1, -4, -2, 7, -5, -3, -3, 0, 10, -2, -3, -1, -4, -1, -8},
            {0, -2, -2, -2, -2, -2, -2, -1, -2, 4, 2, -2, 2, -1, -1, -1, 0, -6, -2, 4, -2, 2, -2, -1, -8},
            {0, -1, 2, 3, -4, 1, 3, 0, 1, -2, -3, 1, -2, -4, -1, 0, 0, -5, -3, -2, 3, -3, 2, -1, -8},
            {-1, -3, -3, -3, -5, -2, -3, -4, -2, 3, 5, -3, 3, 2, -2, -2, -1, -3, -1, 2, -3, 5, -2, -1, -8},
            {0, 0, 1, 3, -5, 3, 3, 0, 2, -2, -3, 0, -2, -5, 0, 0, -1, -6, -4, -2, 2, -2, 3, -1, -8},
            {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -8},
            {-8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, 1}
    };

    /**
     * Returns an amino-acid representation of the passed `codon`.
     *
     * @param codon              {@link String} representing the nucleotide codon to translate.
     * @param asAA3              {@link Boolean} whether the codon should be translated to the amino-acid three-letter code or not.
     * @param includeTermination {@link Boolean} whether a representation for translation termination (in the case of
     *                           three-letter code 'TER' and in the case of one-letter code 'Z') shall
     *                           be returned or not, i.e. an empty String is returned instead.
     * @param includeIncomplete  {@link Boolean} whether incomplete codons, i.e. with a length other than 3, should be
     *                           translated as incomplete amino-acids.
     * @return {@link String} representing a translated amino-acid.
     * @throws MusialException If the passed codon has a length other than three and `includeIncomplete` is false.
     */
    public static String translateCodon(String codon, boolean asAA3, boolean includeTermination,
                                        boolean includeIncomplete)
            throws MusialException {
        if (codon.length() != 3) {
            if (!includeIncomplete) {
                throw new MusialException("Unable to translate codon " + codon + " with length different from three.");
            } else {
                if (asAA3) {
                    return Constants.ANY_AMINOACID3_STRING;
                } else {
                    return Constants.ANY_AMINOACID_STRING;
                }
            }
        } else if (codon.contains(Constants.ANY_NUCLEOTIDE_STRING)) {
            if (asAA3) {
                return Constants.ANY_AMINOACID3_STRING;
            } else {
                return Constants.ANY_AMINOACID_STRING;
            }
        } else if (CODON_MAP.get(codon).equals("") && includeTermination) {
            if (asAA3) {
                return Constants.TERMINATION_AMINOACID3_STRING;
            } else {
                return Constants.TERMINATION_AMINOACID_STRING;
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
     * @throws MusialException If any codon with a length different from three is detected.
     */
    public static String translateNucSequence(Iterable<String> splitNucSequence, boolean includeTermination,
                                              boolean includeIncomplete) throws MusialException {
        StringBuilder translatedNucSequenceBuilder = new StringBuilder();
        for (String s : splitNucSequence) {
            if (s.length() != 3) {
                if (!includeIncomplete) {
                    throw new MusialException("Failed to translate nucleotide sequence containing codon of length unequal 3.");
                }
            } else {
                translatedNucSequenceBuilder.append(translateCodon(s, false, includeTermination, includeIncomplete));
            }
        }
        return translatedNucSequenceBuilder.toString();
    }

    /**
     * Translates a nucleotide sequence into a single-letter amino acid sequence.
     *
     * @param nucSequence        {@link String} representing a nucleotide sequence.
     * @param includeTermination {@link Boolean} whether a representation for translation termination (in the case of
     *                           three-letter code 'TER' and in the case of one-letter code 'Z') shall
     *                           be returned or not, i.e. an empty String is returned instead.
     * @param includeIncomplete  {@link Boolean} whether an incomplete amino-acid should be added to the end if the
     *                           sequence contains an incomplete codon at the end.
     * @param asSense            {@link Boolean} whether the sequence shall be translated as sense or anti-sense.
     * @return {@link String} representing the translated nucleotide sequence.
     * @throws MusialException If any codon with a length different from three is detected.
     */
    public static String translateNucSequence(String nucSequence, boolean includeTermination,
                                              boolean includeIncomplete, boolean asSense) throws MusialException {
        if (!asSense) {
            nucSequence = reverseComplement(nucSequence);
        }
        Iterable<String> splitNucSequence = Splitter.fixedLength(3).split(nucSequence);
        return translateNucSequence(splitNucSequence, includeTermination, includeIncomplete);
    }

    /**
     * Returns the reverse complement of the passed nucleotide sequence.
     *
     * @param sequence {@link String} representing a nucleotide sequence.
     * @return {@link String} representing the reverse complement of the passed sequence.
     */
    public static String reverseComplement(String sequence) {
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
     * symbols the identity is returned.
     *
     * @param base {@link Character} base to invert.
     * @return {@link Character} the inverted base.
     */
    public static Character invertBase(char base) {
        return switch (base) {
            case 'A' -> 'T';
            case 'C' -> 'G';
            case 'G' -> 'C';
            case 'T' -> 'A';
            default -> base;
        };
    }

    /**
     * Invokes the {@link SequenceOperations#globalSequenceAlignment} method with pre-specified parameters for nucleotide sequence
     * alignment.
     * <p>
     * - Uses a simple substitution matrix that scores matches with 1 and mismatches with -1.
     * - Uses a gap open and extension penalty of -2 and -1, respectively.
     *
     * @param nucSeq1    {@link String} representation of the first nucleotide sequence for alignment.
     * @param nucSeq2    {@link String} representation of the second nucleotide sequence for alignment.
     * @param left_mode  {@link GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES} value to indicate how to handle left-marginal gaps.
     * @param right_mode {@link GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES} value to indicate how to handle right-marginal gaps.
     * @param bandWidth  {@link Integer} specifying the band-width for banded alignment or null for non-banded alignment.
     * @return {@link Triplet} storing the alignment score, the aligned first sequence and the aligned second sequence.
     */
    public static Triplet<Integer, String, String> globalNucleotideSequenceAlignment(String nucSeq1, String nucSeq2,
                                                                                     GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES left_mode,
                                                                                     GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES right_mode,
                                                                                     Integer bandWidth) {
        HashMap<Character, Integer> simpleNucleotideScoringMatrixIndexMap = new HashMap<>() {{
            put('A', 0);
            put('C', 1);
            put('G', 2);
            put('T', 3);
            put('N', 4);
        }};
        int[][] simpleNucleotideScoringMatrix = {
                {1, -1, -1, -1, -1},
                {-1, 1, -1, -1, -1},
                {-1, -1, 1, -1, -1},
                {-1, -1, -1, 1, -1},
                {-1, -1, -1, -1, -1},
        };
        return globalSequenceAlignment(nucSeq1, nucSeq2, simpleNucleotideScoringMatrixIndexMap,
                simpleNucleotideScoringMatrix,
                2, 1, left_mode, right_mode, bandWidth);
    }

    /**
     * Invokes the {@link SequenceOperations#globalSequenceAlignment} method with pre-specified parameters for amino-acid sequence
     * alignment.
     * <p>
     * Uses the PAM250 substitution matrix.
     *
     * @param aaSeq1           {@link String} representation of the first amino-acid sequence for alignment.
     * @param aaSeq2           {@link String} representation of the second amino-acid sequence for alignment.
     * @param gapOpenPenalty   {@link Integer} (positive!) to use as gap open penalty.
     * @param gapExtendPenalty {@link Integer} (positive!) to use as gap extension penalty.
     * @param left_mode        {@link GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES} value to indicate how to handle left-marginal gaps.
     * @param right_mode       {@link GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES} value to indicate how to handle right-marginal gaps.
     * @param bandWidth        {@link Integer} specifying the band-width for banded alignment or null for non-banded alignment.
     * @return {@link Triplet} storing the alignment score, the aligned first sequence and the aligned second sequence.
     */
    public static Triplet<Integer, String, String> globalAminoAcidSequenceAlignment(String aaSeq1, String aaSeq2,
                                                                                    int gapOpenPenalty,
                                                                                    int gapExtendPenalty,
                                                                                    GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES left_mode,
                                                                                    GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES right_mode,
                                                                                    Integer bandWidth) {
        return globalSequenceAlignment(aaSeq1, aaSeq2, AA1_PAM250_INDEX, PAM250, gapOpenPenalty, gapExtendPenalty,
                left_mode, right_mode, bandWidth);
    }

    /**
     * Computes a global sequence alignment using the gap-affine Needleman-Wunsch algorithm.
     *
     * @param seq1                  {@link String} representation of the first sequence for alignment.
     * @param seq2                  {@link String} representation of the second sequence for alignment.
     * @param scoringMatrixIndexMap {@link HashMap} mapping characters that may occur in the aligned sequences to integer indices of the scoring matrix.
     * @param scoringMatrix         {@link Integer[][]} containing the scoring values for each pair of characters that may occur in the aligned sequences.
     * @param gapOpenPenalty        {@link Integer} value used to penalize the opening of a gap.
     * @param gapExtendPenalty      {@link Integer} value used to penalize the extension of a gap.
     * @param left_mode             {@link GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES} value to indicate how to handle left-marginal gaps.
     * @param right_mode            {@link GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES} value to indicate how to handle right-marginal gaps.
     * @param bandWidth             {@link Integer} specifying the band-width for banded alignment or null for non-banded alignment.
     * @return {@link Triplet} storing the alignment score, the aligned first sequence and the aligned second sequence.
     */
    private static Triplet<Integer, String, String> globalSequenceAlignment(String seq1, String seq2,
                                                                            HashMap<Character, Integer> scoringMatrixIndexMap,
                                                                            int[][] scoringMatrix, int gapOpenPenalty,
                                                                            int gapExtendPenalty,
                                                                            GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES left_mode,
                                                                            GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES right_mode,
                                                                            Integer bandWidth) {
    /*
    Notes on the alignment matrices:
    - The y-axis will yield seq1.
    - The x-axis will yield seq2
    - Indels are wrt. seq1; thus insertions are traced back by walking vertically and deletions are traced back by
    walking horizontally in the matrix.
     */
        @SuppressWarnings("DuplicatedCode")
        double[][] alignmentMatrix = new double[seq1.length() + 1][seq2.length() + 1];
        Arrays.stream(alignmentMatrix).forEach(row -> Arrays.fill(row, Double.NEGATIVE_INFINITY));
        double[][] matchScoreMatrix = new double[seq1.length() + 1][seq2.length() + 1];
        Arrays.stream(matchScoreMatrix).forEach(row -> Arrays.fill(row, Double.NEGATIVE_INFINITY));
        double[][] insertionScoreMatrix = new double[seq1.length() + 1][seq2.length() + 1];
        Arrays.stream(insertionScoreMatrix).forEach(row -> Arrays.fill(row, Double.NEGATIVE_INFINITY));
        double[][] deletionScoreMatrix = new double[seq1.length() + 1][seq2.length() + 1];
        Arrays.stream(deletionScoreMatrix).forEach(row -> Arrays.fill(row, Double.NEGATIVE_INFINITY));
        char[][] tracebackMatrix = new char[seq1.length() + 1][seq2.length() + 1];
        char[] seq1Array = seq1.toCharArray();
        char[] seq2Array = seq2.toCharArray();
        int alignmentScore;
        int capacity = Math.max(seq1.length(), seq2.length());
        StringBuilder seq1Builder = new StringBuilder(capacity);
        StringBuilder seq2Builder = new StringBuilder(capacity);
        /*
    (1) Compute global sequence alignment.
     */
        alignmentMatrix[0][0] = 0;
        matchScoreMatrix[0][0] = 0;
        insertionScoreMatrix[0][0] = 0;
        deletionScoreMatrix[0][0] = 0;
        int gapCost;
        // i -> PREFIX
        for (int i = 1; i < (Objects.isNull(bandWidth) ? seq1.length() + 1 : bandWidth + 1); i++) {
            gapCost = switch (left_mode) {
                case FREE -> 0;
                case PENALIZE -> -gapOpenPenalty - (i - 1) * gapExtendPenalty;
                case FORBID -> -gapOpenPenalty * seq1.length();
            };
            alignmentMatrix[i][0] = gapCost;
            matchScoreMatrix[i][0] = gapCost;
            insertionScoreMatrix[i][0] = gapCost;
            deletionScoreMatrix[i][0] = gapCost;
            tracebackMatrix[i][0] = 'I';
        }
        // j -> SUFFIX
        for (int j = 1; j < (Objects.isNull(bandWidth) ? seq2.length() + 1 : bandWidth + 1); j++) {
            gapCost = switch (right_mode) {
                case FREE -> 0;
                case PENALIZE -> -gapOpenPenalty - (j - 1) * gapExtendPenalty;
                case FORBID -> -gapOpenPenalty * seq2.length();
            };
            alignmentMatrix[0][j] = gapCost;
            matchScoreMatrix[0][j] = gapCost;
            insertionScoreMatrix[0][j] = gapCost;
            deletionScoreMatrix[0][j] = gapCost;
            tracebackMatrix[0][j] = 'D';
        }
        double max;
        for (int i = 1; i < seq1.length() + 1; i++) {
            int jLeftBound = Math.max(
                    1,
                    (Objects.isNull(bandWidth) ? 1 : i - bandWidth)
            );
            int jRightBound = Math.min(
                    seq2.length(),
                    (Objects.isNull(bandWidth) ? seq2.length() : i + bandWidth)
            );
            for (int j = jLeftBound; j < jRightBound + 1; j++) {
                matchScoreMatrix[i][j] =
                        alignmentMatrix[i - 1][j - 1]
                                + scoringMatrix[scoringMatrixIndexMap.get(seq1Array[i - 1])][scoringMatrixIndexMap
                                .get(seq2Array[j - 1])];
                insertionScoreMatrix[i][j] =
                        Math.max(
                                alignmentMatrix[i - 1][j] - gapOpenPenalty,
                                insertionScoreMatrix[i - 1][j] - gapExtendPenalty
                        );
                deletionScoreMatrix[i][j] =
                        Math.max(
                                alignmentMatrix[i][j - 1] - gapOpenPenalty,
                                deletionScoreMatrix[i][j - 1] - gapExtendPenalty
                        );
                max = Integer.MIN_VALUE;
                if (insertionScoreMatrix[i][j] > max) {
                    max = insertionScoreMatrix[i][j];
                    alignmentMatrix[i][j] = (int) max;
                    tracebackMatrix[i][j] = 'I';
                }
                if (deletionScoreMatrix[i][j] > max) {
                    max = deletionScoreMatrix[i][j];
                    alignmentMatrix[i][j] = (int) max;
                    tracebackMatrix[i][j] = 'D';
                }
                if (matchScoreMatrix[i][j] > max) {
                    max = matchScoreMatrix[i][j];
                    alignmentMatrix[i][j] = (int) max;
                    tracebackMatrix[i][j] = 'M';
                }
            }
        }
        alignmentScore = (int) alignmentMatrix[seq1.length()][seq2.length()];
    /*
    (2) Compute traceback path from global alignment.
    */
        LinkedList<Character> tracebackPath = new LinkedList<>();
        boolean constructTracebackPath = true;
        int i = seq1.length();
        int j = seq2.length();
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
                if (tracebackDirection == '\u0000') // Prevents traceback to walk out of the set alignment band.
                    tracebackDirection = 'M';
            }
        }
        Collections.reverse(tracebackPath);
    /*
    (3) Iterate over traceback path and construct aligned sequences.
    */
        int seq1Index = 0;
        int seq2Index = 0;
        for (Character character : tracebackPath) {
            // Start walking along the traceback path.
            tracebackDirection = character;
            if (tracebackDirection == 'M') {
                // CASE: Match or mismatch of nucleotide and amino-acid sequence.
                seq1Builder.append(seq1Array[seq1Index]);
                seq1Index += 1;
                seq2Builder.append(seq2Array[seq2Index]);
                seq2Index += 1;
            } else if (tracebackDirection == 'D') {
                // CASE: Deletion wrt. to first amino acid sequence.
                seq1Builder.append(Constants.DELETION_OR_GAP_STRING);
                seq2Builder.append(seq2Array[seq2Index]);
                seq2Index += 1;
            } else if (tracebackDirection == 'I') {
                // CASE: Insertion wrt. to first amino acid sequence.
                seq1Builder.append(seq1Array[seq1Index]);
                seq1Index += 1;
                seq2Builder.append(Constants.DELETION_OR_GAP_STRING);
            }
        }
    /*
    (6) Insert results into Triplet.
     */
        return new Triplet<>(alignmentScore, seq1Builder.toString(), seq2Builder.toString());
    }

    /**
     * Computes variants from two aligned sequences.
     * <p>
     * Variants are extracted wrt. the query versus the target sequence and formatted as p:a:r, where p is the variant
     * position, r is the reference content and a is the alternate content.
     *
     * @param targetSequence {@link String} representation of the target sequence (i.e. reference).
     * @param querySequence  {@link String} representation of the query sequence (i.e. the one with variants).
     * @return {@link ArrayList} containing derived variants, c.f. method description for format details.
     */
    @SuppressWarnings("DuplicatedCode")
    public static ArrayList<Triple<Integer, String, String>> getVariantsOfAlignedSequences(String targetSequence,
                                                                                           String querySequence) {
        ArrayList<Triple<Integer, String, String>> variants = new ArrayList<>((int) Math.ceil(targetSequence.length() * 0.5));
        StringBuilder variantBuilder = new StringBuilder();
        StringBuilder referenceBuilder = new StringBuilder();
        char[] targetSequenceArray = targetSequence.toCharArray();
        char[] querySequenceArray = querySequence.toCharArray();
        char targetContent;
        char queryContent;
        int variantStart = 1;
        int lastNonGapPosition = 1;
        int alignmentLength = targetSequence.length();
        boolean isSubstitution = false;
        boolean isInsertion = false;
        boolean isDeletion = false;
        boolean ambiguousSwitch = false;
        for (int i = 0; i < alignmentLength; i++) {
            targetContent = targetSequenceArray[i];
            queryContent = querySequenceArray[i];
            if (targetContent == queryContent) {
                // CASE: Match.
                if (isSubstitution || isInsertion || isDeletion) {
                    variants.add(Triple.of(variantStart, referenceBuilder.toString(), variantBuilder.toString()));
                    referenceBuilder.setLength(0);
                    variantBuilder.setLength(0);
                }
                ambiguousSwitch = false;
                isSubstitution = false;
                isInsertion = false;
                isDeletion = false;
                lastNonGapPosition = i + 1;
            } else if (targetContent == Constants.DELETION_OR_GAP_CHAR) {
                // CASE: Insertion (in variant).
                if (isDeletion) {
                    if (!ambiguousSwitch)
                        Logger.logWarning("Ambiguous Deletion to Insertion switch in aligned sequences ..." + referenceBuilder + targetContent + "...\t..." + variantBuilder + queryContent + "...; Variant will be skipped.");
                    ambiguousSwitch = true;
                } else {
                    if (isSubstitution) {
                        variants.add(Triple.of(variantStart, referenceBuilder.toString(), variantBuilder.toString()));
                        referenceBuilder.setLength(0);
                        variantBuilder.setLength(0);
                    }
                    if (!isInsertion) {
                        variantStart = lastNonGapPosition;
                        referenceBuilder.append(targetSequenceArray[lastNonGapPosition - 1]);
                        variantBuilder.append(targetSequenceArray[lastNonGapPosition - 1]);
                    }
                    referenceBuilder.append(Constants.DELETION_OR_GAP_CHAR);
                    variantBuilder.append(queryContent);
                    isSubstitution = false;
                    isInsertion = true;
                }
            } else if (queryContent == Constants.DELETION_OR_GAP_CHAR) {
                // CASE: Deletion (in variant).
                if (isInsertion) {
                    if (!ambiguousSwitch)
                        Logger.logWarning("Ambiguous Insertion to Deletion switch in aligned sequences ..." + referenceBuilder + targetContent + "...\t..." + variantBuilder + queryContent + "...; Variant will be skipped.");
                    ambiguousSwitch = true;
                } else {
                    if (isSubstitution) {
                        variants.add(Triple.of(variantStart, referenceBuilder.toString(), variantBuilder.toString()));
                        referenceBuilder.setLength(0);
                        variantBuilder.setLength(0);
                    }
                    if (!isDeletion) {
                        variantStart = lastNonGapPosition;
                        referenceBuilder.append(targetSequenceArray[lastNonGapPosition - 1]);
                        variantBuilder.append(targetSequenceArray[lastNonGapPosition - 1]);
                    }
                    referenceBuilder.append(targetContent);
                    variantBuilder.append(Constants.DELETION_OR_GAP_CHAR);
                    isSubstitution = false;
                    isDeletion = true;
                }
            } else {
                // CASE: Mismatch/Substitution.
                if (isSubstitution || isInsertion || isDeletion) {
                    variants.add(Triple.of(variantStart, referenceBuilder.toString(), variantBuilder.toString()));
                    referenceBuilder.setLength(0);
                    variantBuilder.setLength(0);
                }
                ambiguousSwitch = false;
                isSubstitution = true;
                isInsertion = false;
                isDeletion = false;
                lastNonGapPosition = i + 1;
                variantStart = lastNonGapPosition;
                variantBuilder.append(queryContent);
                referenceBuilder.append(targetContent);
            }
        }
        if (isSubstitution || isInsertion || isDeletion) {
            variants.add(Triple.of(variantStart, referenceBuilder.toString(), variantBuilder.toString()));
            variantBuilder.setLength(0);
            referenceBuilder.setLength(0);
        }
        variants.trimToSize();
        return variants;
    }

    /**
     * Enum to store different modes to handle prefix gaps for global sequence alignment.
     */
    public enum GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES {
        /**
         * Gaps at the end are not penalized.
         */
        FREE,
        /**
         * Gaps at the end are penalized normally.
         */
        PENALIZE,
        /**
         * Gaps at the end are not allowed.
         */
        FORBID
    }
}
