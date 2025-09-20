package util;

import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import org.apache.commons.lang3.tuple.Triple;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AmbiguityRNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.transcription.Frame;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.NavigableMap;

/**
 * Utility class for performing various sequence related operations.
 * <p>
 * This class provides static methods for sequence alignment, variant integration, sequence translation, and other related operations. It
 * includes methods for handling nucleotide and protein sequences, as well as utilities for working with gaps and variants.
 */
public final class Bio {

    /**
     * Private constructor to prevent instantiation of this utility class.
     */
    private Bio() {
    }

    /**
     * A transcription engine for translating DNA sequences.
     * <p>
     * See <a href="https://github.com/biojava/biojava-tutorial/blob/master/core/translating.md">https://github
     * .com/biojava/biojava-tutorial/blob/master/core/translating.md</a>
     */
    private static final TranscriptionEngine transcriptionEngine = new TranscriptionEngine.Builder()
            .dnaCompounds(AmbiguityDNACompoundSet.getDNACompoundSet())
            .rnaCompounds(AmbiguityRNACompoundSet.getRNACompoundSet())
            .build();

    /**
     * A cache for storing translated DNA sequences.
     * <p>
     * This static {@link HashMap} is used to store previously translated DNA sequences to improve performance by avoiding redundant
     * translations. The key is the hash code of the translation request (including sequence and direction), and the value is the translated
     * amino acid sequence.
     */
    private static final HashMap<Integer, String> translationCache = new HashMap<>();

    /**
     * Computes optimal pairwise global nucleotide sequence alignment using a gap-affine (Gotoh) banded Needleman-Wunsch algorithm.
     * <p>
     * A simple scoring matrix (match: +1; mismatch: -1) is used.
     *
     * @param sequenceA        The first nucleotide sequence to align.
     * @param sequenceB        The second nucleotide sequence to align.
     * @param gapOpenPenalty   The penalty for opening a gap in the alignment.
     * @param gapExtendPenalty The penalty for extending an existing gap in the alignment.
     * @param noGapPrefix      Prevent (if true) gaps at the beginning of the aligned sequences.
     * @param noGapSuffix      Prevent (if true) gaps at the end of the aligned sequences.
     * @param bandWidth        The width of the band for banded alignment; if below or equal to 0, the full length of sequence B is used.
     * @return A {@link Tuple} containing the aligned sequences.
     */
    public static Tuple<String, String> globalNucleotideSequenceAlignment(String sequenceA, String sequenceB, int gapOpenPenalty,
                                                                          int gapExtendPenalty,
                                                                          boolean noGapPrefix, boolean noGapSuffix, int bandWidth) {
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
        return globalSequenceAlignment(sequenceA, sequenceB, simpleNucleotideScoringMatrixIndexMap, simpleNucleotideScoringMatrix,
                gapOpenPenalty, gapExtendPenalty, noGapPrefix, noGapSuffix, bandWidth);
    }

    /**
     * Computes optimal pairwise global amino acid sequence alignment using a gap-affine (Gotoh) banded Needleman-Wunsch algorithm.
     * <p>
     * This method uses the BLOSUM80 scoring matrix for amino acid matches and mismatches.
     *
     * @param sequenceA        The first protein sequence to align.
     * @param sequenceB        The second protein sequence to align.
     * @param gapOpenPenalty   The penalty for opening a gap in the alignment.
     * @param gapExtendPenalty The penalty for extending an existing gap in the alignment.
     * @param noGapPrefix      Prevent (if true) gaps at the beginning of the aligned sequences.
     * @param noGapSuffix      Prevent (if true) gaps at the end of the aligned sequences.
     * @param bandWidth        The width of the band for banded alignment; if below or equal to 0, the full length of sequence B is used.
     * @return A {@link Tuple} containing the aligned sequences.
     */
    public static Tuple<String, String> globalProteinSequenceAlignment(String sequenceA, String sequenceB, int gapOpenPenalty,
                                                                       int gapExtendPenalty,
                                                                       boolean noGapPrefix, boolean noGapSuffix, int bandWidth) {
        HashMap<Character, Integer> blosum80IndexMap = new HashMap<>() {{
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
            put('X', 23);
            put('*', 24);
        }};
        int[][] blosum80 = {
                {5, -2, -2, -2, -1, -1, -1, 0, -2, -2, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0, -2, -2, -1, -1, -6},
                {-2, 6, -1, -2, -4, 1, -1, -3, 0, -3, -3, 2, -2, -4, -2, -1, -1, -4, -3, -3, -1, -3, 0, -1, -6},
                {-2, -1, 6, 1, -3, 0, -1, -1, 0, -4, -4, 0, -3, -4, -3, 0, 0, -4, -3, -4, 5, -4, 0, -1, -6},
                {-2, -2, 1, 6, -4, -1, 1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4, 5, -5, 1, -1, -6},
                {-1, -4, -3, -4, 9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -4, -2, -4, -1, -6},
                {-1, 1, 0, -1, -4, 6, 2, -2, 1, -3, -3, 1, 0, -4, -2, 0, -1, -3, -2, -3, 0, -3, 4, -1, -6},
                {-1, -1, -1, 1, -5, 2, 6, -3, 0, -4, -4, 1, -2, -4, -2, 0, -1, -4, -3, -3, 1, -4, 5, -1, -6},
                {0, -3, -1, -2, -4, -2, -3, 6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -1, -5, -3, -1, -6},
                {-2, 0, 0, -2, -4, 1, 0, -3, 8, -4, -3, -1, -2, -2, -3, -1, -2, -3, 2, -4, -1, -4, 0, -1, -6},
                {-2, -3, -4, -4, -2, -3, -4, -5, -4, 5, 1, -3, 1, -1, -4, -3, -1, -3, -2, 3, -4, 3, -4, -1, -6},
                {-2, -3, -4, -5, -2, -3, -4, -4, -3, 1, 4, -3, 2, 0, -3, -3, -2, -2, -2, 1, -4, 3, -3, -1, -6},
                {-1, 2, 0, -1, -4, 1, 1, -2, -1, -3, -3, 5, -2, -4, -1, -1, -1, -4, -3, -3, -1, -3, 1, -1, -6},
                {-1, -2, -3, -4, -2, 0, -2, -4, -2, 1, 2, -2, 6, 0, -3, -2, -1, -2, -2, 1, -3, 2, -1, -1, -6},
                {-3, -4, -4, -4, -3, -4, -4, -4, -2, -1, 0, -4, 0, 6, -4, -3, -2, 0, 3, -1, -4, 0, -4, -1, -6},
                {-1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4, 8, -1, -2, -5, -4, -3, -2, -4, -2, -1, -6},
                {1, -1, 0, -1, -2, 0, 0, -1, -1, -3, -3, -1, -2, -3, -1, 5, 1, -4, -2, -2, 0, -3, 0, -1, -6},
                {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2, 1, 5, -4, -2, 0, -1, -1, -1, -1, -6},
                {-3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2, 0, -5, -4, -4, 11, 2, -3, -5, -3, -3, -1, -6},
                {-2, -3, -3, -4, -3, -2, -3, -4, 2, -2, -2, -3, -2, 3, -4, -2, -2, 2, 7, -2, -3, -2, -3, -1, -6},
                {0, -3, -4, -4, -1, -3, -3, -4, -4, 3, 1, -3, 1, -1, -3, -2, 0, -3, -2, 4, -4, 2, -3, -1, -6},
                {-2, -1, 5, 5, -4, 0, 1, -1, -1, -4, -4, -1, -3, -4, -2, 0, -1, -5, -3, -4, 5, -4, 0, -1, -6},
                {-2, -3, -4, -5, -2, -3, -4, -5, -4, 3, 3, -3, 2, 0, -4, -3, -1, -3, -2, 2, -4, 3, -3, -1, -6},
                {-1, 0, 0, 1, -4, 4, 5, -3, 0, -4, -3, 1, -1, -4, -2, 0, -1, -3, -3, -3, 0, -3, 5, -1, -6},
                {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6},
                {-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1}
        };
        return globalSequenceAlignment(sequenceA, sequenceB, blosum80IndexMap, blosum80, gapOpenPenalty, gapExtendPenalty, noGapPrefix,
                noGapSuffix, bandWidth);
    }

    /**
     * Computes optimal pairwise global sequence alignment using a gap-affine (Gotoh) banded Needleman-Wunsch algorithm.
     * <p>
     * Setting the bandWidth to a value less than or equal to 0 will result in a full alignment.
     * <p>
     * Optionally, gaps can be disallowed at the beginning and/or end of the aligned sequences.
     *
     * @param sequenceA         The first sequence to align.
     * @param sequenceB         The second sequence to align.
     * @param symbolsScoreIndex A mapping of characters to their respective indices in the scoring matrix.
     * @param scores            A 2D array representing the scoring matrix for character matches and mismatches.
     * @param gapOpenPenalty    The penalty for opening a gap in the alignment.
     * @param gapExtendPenalty  The penalty for extending an existing gap in the alignment.
     * @param noGapPrefix       Prevent (if true) gaps at the beginning of the aligned sequences.
     * @param noGapSuffix       Prevent (if true) gaps at the end of the aligned sequences.
     * @param bandWidth         The width of the band for banded alignment. Values less than or equal to 0 will result in a full alignment.
     * @return A {@link Tuple} containing the aligned sequences.
     * @throws IllegalArgumentException If the bandWidth is too narrow for the given sequences.
     */
    private static Tuple<String, String> globalSequenceAlignment(String sequenceA, String sequenceB,
                                                                 HashMap<Character, Integer> symbolsScoreIndex,
                                                                 int[][] scores, int gapOpenPenalty,
                                                                 int gapExtendPenalty,
                                                                 boolean noGapPrefix,
                                                                 boolean noGapSuffix,
                                                                 int bandWidth) {
        // Check if either sequence is empty and return accordingly.
        if (sequenceA.isEmpty())
            return new Tuple<>(Constants.GAP.repeat(sequenceB.length()), sequenceB);
        if (sequenceB.isEmpty())
            return new Tuple<>(sequenceA, Constants.GAP.repeat(sequenceA.length()));

        // Define constants for the algorithm.
        final int MIN = Integer.MIN_VALUE / 2;
        int lengthA = sequenceA.length(), lengthB = sequenceB.length(); // Lengths of the sequences to align.
        boolean useBand = bandWidth > 0; // Flag to indicate if banded alignment is used.

        // Validate band-width for banded alignment.
        if (useBand && bandWidth < Math.abs(lengthA - lengthB)) {
            throw new IllegalArgumentException("Banded alignment width (%d) is too narrow for sequences of lengths %d and %d.".formatted(
                    bandWidth, lengthA, lengthB));
        }

        // Initialize sequence components.
        int maxLength = Math.max(lengthA, sequenceB.length()); // This may be exceeded by inserted gaps.
        StringBuilder alignedA = new StringBuilder(maxLength);
        StringBuilder alignedB = new StringBuilder(maxLength);
        char[] symbolsA = sequenceA.toCharArray();
        char[] symbolsB = sequenceB.toCharArray();

        // Initialize DP matrices.
        int[][] D = new int[lengthA + 1][lengthB + 1];
        int[][] P = new int[lengthA + 1][lengthB + 1];
        int[][] Q = new int[lengthA + 1][lengthB + 1];

        // Fill dynamic programming matrices.
        int leftBoundary = 0;
        int rightBoundary = lengthB;
        for (int i = 0; i <= lengthA; i++) {
            // Set the boundaries for banded alignment.
            if (useBand) {
                leftBoundary = Math.max(0, i - bandWidth);
                rightBoundary = Math.min(lengthB, i + bandWidth);
            }
            for (int j = leftBoundary; j <= rightBoundary; j++) {
                // Set cell values outside the band to a worse-than minimum value.
                if (useBand) {
                    if (j > 0 && j == leftBoundary) {
                        D[i][j - 1] = MIN - 1;
                        P[i][j - 1] = MIN - 1;
                        Q[i][j - 1] = MIN - 1;
                    } else if (j < lengthB && j == rightBoundary) {
                        D[i][j + 1] = MIN - 1;
                        P[i][j + 1] = MIN - 1;
                        Q[i][j + 1] = MIN - 1;
                    }
                }
                // Fill the DP matrices based on the current indices.
                if (i == 0 && j == 0) {
                    // The first cell is initialized to zero.
                    D[i][j] = 0;
                    P[i][j] = 0;
                    Q[i][j] = 0;
                } else if (i == 0 && j > 0) {
                    // First row.
                    if (noGapPrefix) D[i][j] = MIN;
                    else D[i][j] = -gapOpenPenalty - j * gapExtendPenalty;
                    P[i][j] = MIN;
                    Q[i][j] = MIN;
                } else if (i > 0 && j == 0) {
                    // First column.
                    if (noGapPrefix) D[i][j] = MIN;
                    else D[i][j] = -gapOpenPenalty - i * gapExtendPenalty;
                    P[i][j] = MIN;
                    Q[i][j] = MIN;
                } else {
                    // Other cells.
                    if (noGapSuffix && (i == lengthA ^ j == lengthB)) {
                        // If no gaps are allowed at the end, set cells in last row or column to minimum value.
                        D[i][j] = MIN;
                        P[i][j] = MIN;
                        Q[i][j] = MIN;
                    } else {
                        // Otherwise calculate the scores based on the previous cells as defined by the Gotoh algorithm.
                        P[i][j] = Math.max(MIN, Math.max(
                                D[i - 1][j] - gapOpenPenalty - gapExtendPenalty,
                                P[i - 1][j] - gapExtendPenalty
                        ));
                        Q[i][j] = Math.max(MIN, Math.max(
                                D[i][j - 1] - gapOpenPenalty - gapExtendPenalty,
                                Q[i][j - 1] - gapExtendPenalty
                        ));
                        D[i][j] = Math.max(MIN, Math.max(
                                D[i - 1][j - 1] + scores[symbolsScoreIndex.get(symbolsA[i - 1])][symbolsScoreIndex.get(symbolsB[j - 1])],
                                Math.max(P[i][j], Q[i][j])
                        ));
                    }
                }
            }
        }

        // Traceback through the matrices to construct aligned sequences.
        int i = sequenceA.length();
        int j = sequenceB.length();
        String matrix = "D"; // Start from D.
        int score;
        while (i > 0 || j > 0) {
            if (i == 0) {
                // Add a gap in sequence A.
                alignedA.append(Constants.GAP);
                alignedB.append(symbolsB[j - 1]);
                j--;
            } else if (j == 0) {
                // Add a gap in sequence B.
                alignedA.append(symbolsA[i - 1]);
                alignedB.append(Constants.GAP);
                i--;
            } else {
                // Determine the current direction and update indices accordingly
                switch (matrix) {
                    case "D" -> {
                        score = D[i - 1][j - 1] + scores[symbolsScoreIndex.get(symbolsA[i - 1])][symbolsScoreIndex.get(symbolsB[j - 1])];
                        if (D[i][j] == score) {
                            alignedA.append(symbolsA[i - 1]);
                            alignedB.append(symbolsB[j - 1]);
                            i--;
                            j--;
                        } else if (D[i][j] == P[i][j]) {
                            matrix = "P";
                        } else if (D[i][j] == Q[i][j]) {
                            matrix = "Q";
                        } else {
                            throw new IllegalStateException("Invalid alignment matrix state at cell D[%d][%d]: %d".formatted(i, j,
                                    D[i][j]));
                        }
                    }
                    case "P" -> {
                        if (P[i][j] == D[i - 1][j] - gapOpenPenalty - gapExtendPenalty) {
                            matrix = "D";
                        }
                        alignedA.append(symbolsA[i - 1]);
                        alignedB.append(Constants.GAP);
                        i--;
                    }
                    case "Q" -> {
                        if (Q[i][j] == D[i][j - 1] - gapOpenPenalty - gapExtendPenalty) {
                            matrix = "D";
                        }
                        alignedA.append(Constants.GAP);
                        alignedB.append(symbolsB[j - 1]);
                        j--;
                    }
                }
            }
        }

        // Reverse the aligned sequences to get the final result and return.
        alignedA.reverse();
        alignedB.reverse();
        return new Tuple<>(alignedA.toString(), alignedB.toString());
    }

    /**
     * Pads a string with gap characters to reach a specified length.
     * <p>
     * This method appends gap characters (defined by {@link Constants#GAP}) to the input string until it reaches the desired length. If the
     * input string is already equal to or longer than the specified length, no padding is added.
     *
     * @param s      The input string to be padded.
     * @param length The desired length of the resulting string.
     * @return The padded string, or the original string if no padding is needed.
     */
    public static String padGaps(String s, int length) {
        return s + Constants.GAP.repeat(Math.max(0, length - s.length()));
    }

    /**
     * Removes all gap characters from the input string.
     * <p>
     * This method replaces all occurrences of the gap character (defined by {@link Constants#GAP}) in the input string with an empty string
     * (defined by {@link Constants#EMPTY}).
     *
     * @param s The input string from which gaps should be removed.
     * @return A new string with all gap characters removed.
     */
    public static String stripGaps(String s) {
        return s.replaceAll(Constants.GAP, Constants.EMPTY);
    }

    /**
     * Integrates variants into a reference sequence.
     * <p>
     * This method modifies the given reference sequence by incorporating the specified variants. Variants are represented as a
     * {@link NavigableMap} where the key is the 0-based position in the reference sequence, and the value is the alternative base
     * sequence.
     * <p>
     * Variants have to be in canonical form, i.e. they must start with a non-gap character that is assumed to be in the coordinate system
     * of the reference sequence. If additional characters follow, these have to be all gaps (indicating a deletion) or all non-gaps
     * (indicating an insertion).
     * <p>
     * The method handles deletions by tracking the number of gaps introduced and ensures that the resulting sequence reflects the
     * integrated variants. Optionally, gaps can be stripped from the final sequence.
     *
     * @param reference The original reference sequence as a {@link String}.
     * @param variants  A {@link NavigableMap} containing the variants to integrate, where the key is the position and the value is the
     *                  alternative base sequence in canonical form.
     * @param stripGaps A {@code boolean} indicating whether to remove gaps from the resulting sequence.
     * @return A {@link String} representing the reference sequence with the integrated variants. If {@code stripGaps} is {@code true}, gaps
     * are removed from the resulting sequence.
     * @throws IllegalArgumentException If the reference sequence is empty.
     */
    public static String integrateVariants(String reference, NavigableMap<Integer, String> variants, boolean stripGaps) {
        // Validate input and return early if no processing is needed.
        if (reference.isEmpty()) throw new IllegalArgumentException("Reference sequence is empty.");
        if (variants.isEmpty()) return reference;

        // Initialize result builder and deletion counter.
        StringBuilder result = new StringBuilder(reference.length());
        int deletions = 0;

        // Process each position in the reference sequence.
        for (int i = 0; i < reference.length(); i++) {
            // Check if a variant exists at the current position.
            if (variants.containsKey(i)) {
                String variant = variants.get(i);
                // Handle deletions by appending gaps if needed.
                if (deletions > 0) {
                    result.append(Constants.GAP);
                    deletions--;
                } else {
                    // Append the first character of the variant.
                    result.append(variant.charAt(0));
                }
                // Handle insertions, deletions, or invalid variants.
                if (isInsertion(variant)) {
                    result.append(variant.substring(1));
                } else if (isDeletion(variant)) {
                    deletions += variant.length() - 1;
                } else if (!isSubstitution(variant)) {
                    throw new IllegalArgumentException("Invalid variant '%s' at position %d.".formatted(variant, i));
                }
            } else {
                // Append gaps or the current character from the reference sequence.
                result.append(deletions > 0 ? Constants.GAP : reference.charAt(i));
                if (deletions > 0) deletions--;
            }
        }

        // Return the final sequence, optionally stripping gaps.
        return stripGaps ? stripGaps(result.toString()) : result.toString();
    }

    /**
     * Translates a DNA sequence into an amino-acid sequence. The translation is always performed in the 1-frame. Utilizes the <a
     * href="https://github.com/biojava/biojava">BioJava library</a> for translation.
     *
     * @param sequence The DNA sequence to translate.
     * @param reverse  Whether to translate the reverse complement of the sequence.
     * @return The translated amino-acid sequence.
     * @throws MusialException If an error occurs during translation.
     */
    public static String translateSequence(String sequence, boolean reverse) throws MusialException {
        if (sequence.isEmpty()) return Constants.EMPTY;
        String cachedTranslationKey = "%s-%s".formatted(reverse ? "rev" : "fwd", sequence);
        // Check if the translation result is already cached.
        if (translationCache.containsKey(cachedTranslationKey.hashCode())) {
            return translationCache.get(cachedTranslationKey.hashCode());
        }
        try {
            // Initialize the DNA sequence.
            Sequence<NucleotideCompound> dna = new DNASequence(sequence);
            String translatedSequence;
            if (reverse)
                translatedSequence =
                        transcriptionEngine.multipleFrameTranslation(dna, Frame.REVERSED_ONE).get(Frame.REVERSED_ONE).getSequenceAsString();
            else
                translatedSequence = transcriptionEngine.multipleFrameTranslation(dna, Frame.ONE).get(Frame.ONE).getSequenceAsString();
            // Add translated stop codon if present at the end of the sequence.
            if (reverse) {
                if (sequence.startsWith("CTA") || sequence.startsWith("TTA") || sequence.startsWith("TCA"))
                    translatedSequence += Constants.TERMINAL_AA;
            } else {
                if (sequence.endsWith("TAG") || sequence.endsWith("TAA") || sequence.endsWith("TGA"))
                    translatedSequence += Constants.TERMINAL_AA;
            }
            // Cache the translation result.
            translationCache.put(cachedTranslationKey.hashCode(), translatedSequence);
            return translatedSequence;
        } catch (Exception e) {
            throw new MusialException("org.biojava.nbio.core.sequence.DNASequence: " + e.getMessage());
        }
    }

    /**
     * Transforms two sequences into canonical VCF variants.
     * <p>
     * The specified reference and alternative are expected to be aligned sequences. Variants are formatted as triples of relative position,
     * reference-, and variant content. The relative position is the 0-based position of the variant in the reference sequence without
     * gaps.
     *
     * @param reference   {@link String} representation of the reference sequence.
     * @param alternative {@link String} representation of the variant/alternative sequence.
     * @return {@link ArrayList} containing derived variants, c.f. method description for format details.
     */
    public static ArrayList<Triple<Integer, String, String>> getCanonicalVariants(String reference, String alternative) {
        if (reference.length() != alternative.length()) {
            throw new IllegalArgumentException("Reference and alternative sequence lengths do not match.");
        }
        ArrayList<Triple<Integer, String, String>> variants = new ArrayList<>();
        StringBuilder referenceBuilder = new StringBuilder();
        StringBuilder alternativeBuilder = new StringBuilder();
        char[] referenceChars = reference.toCharArray();
        char[] alternativeChars = alternative.toCharArray();
        int relativeStartPosition = 0;
        int noInsertions = 0;
        int lastNonGapIndex = 0;
        boolean isSubstitution = false, isInsertion = false, isDeletion = false, ambiguousSwitch = false;
        for (int i = 0; i < reference.length(); i++) {
            if (referenceChars[i] == alternativeChars[i]) { // Match:
                if (isSubstitution || isInsertion || isDeletion) {
                    variants.add(Triple.of(relativeStartPosition, referenceBuilder.toString(), alternativeBuilder.toString()));
                    referenceBuilder.setLength(0);
                    alternativeBuilder.setLength(0);
                }
                isSubstitution = isInsertion = isDeletion = ambiguousSwitch = false;
                lastNonGapIndex = i;
            } else if (referenceChars[i] == Constants.GAP_CHAR) { // Insertion:
                if (isDeletion) {
                    if (!ambiguousSwitch) Logging.logWarning("Skip variant %s > %s due to ambiguous deletion to insertion switch."
                            .formatted(reference, alternative));
                    ambiguousSwitch = true;
                } else {
                    if (!isInsertion && !isSubstitution) {
                        relativeStartPosition = lastNonGapIndex - noInsertions;
                        referenceBuilder.append(referenceChars[lastNonGapIndex]);
                        alternativeBuilder.append(alternativeChars[lastNonGapIndex]);
                    }
                    referenceBuilder.append(Constants.GAP_CHAR);
                    alternativeBuilder.append(alternativeChars[i]);
                    isSubstitution = false;
                    isInsertion = true;
                }
                noInsertions++;
            } else if (alternativeChars[i] == Constants.GAP_CHAR) { // Deletion
                if (isInsertion) {
                    if (!ambiguousSwitch) Logging.logWarning("Skip variant %s > %s due to ambiguous deletion to insertion switch."
                            .formatted(reference, alternative));
                    ambiguousSwitch = true;
                } else {
                    if (!isDeletion && !isSubstitution) {
                        relativeStartPosition = lastNonGapIndex - noInsertions;
                        referenceBuilder.append(referenceChars[lastNonGapIndex]);
                        alternativeBuilder.append(alternativeChars[lastNonGapIndex]);
                    }
                    referenceBuilder.append(referenceChars[i]);
                    alternativeBuilder.append(Constants.GAP_CHAR);
                    isSubstitution = false;
                    isDeletion = true;
                }
            } else { // Substitution
                if (isSubstitution || isInsertion || isDeletion) {
                    variants.add(Triple.of(relativeStartPosition, referenceBuilder.toString(), alternativeBuilder.toString()));
                    referenceBuilder.setLength(0);
                    alternativeBuilder.setLength(0);
                }
                isSubstitution = true;
                isInsertion = isDeletion = ambiguousSwitch = false;
                lastNonGapIndex = i;
                relativeStartPosition = lastNonGapIndex - noInsertions;
                referenceBuilder.append(referenceChars[i]);
                alternativeBuilder.append(alternativeChars[i]);
            }
        }
        if (isSubstitution || isInsertion || isDeletion) {
            variants.add(Triple.of(relativeStartPosition, referenceBuilder.toString(), alternativeBuilder.toString()));
        }
        return variants;
    }

    /**
     * Determines whether a variant is a substitution; i.e., both the reference and alternative base content match a single base of
     * {@link Constants#BASE_SYMBOLS}.
     *
     * @param ref The reference base content.
     * @param alt The alternative base content.
     * @return {@code true} if the variant is a substitution, {@code false} otherwise.
     */
    public static boolean isSubstitution(String ref, String alt) {
        return isSubstitution(ref) && isSubstitution(alt);
    }

    /**
     * Determines whether a given alternative base content represents a substitution.
     * <p>
     * A substitution is defined as a single base from the set of valid nucleotide symbols defined in {@link Constants#BASE_SYMBOLS}.
     *
     * @param alt The alternative base content to check.
     * @return {@code true} if the alternative content represents a substitution, {@code false} otherwise.
     */
    public static boolean isSubstitution(String alt) {
        return alt.length() == 1 && Constants.BASE_SYMBOLS.indexOf(alt.charAt(0)) != -1;
    }

    /**
     * Determines whether a variant is an insertion, i.e.,
     * <ul>
     *     <li>either the alternative base content is a string of any length of {@link Constants#BASE_SYMBOLS}
     *     and the reference base content is a single base of {@link Constants#BASE_SYMBOLS} followed by {@link Constants#GAP}s
     *     matching the alternative content's length (padded canonical),</li>
     *     <li>or the reference base content is a single base of {@link Constants#BASE_SYMBOLS} and the alternative
     *     content is a string of any length of {@link Constants#BASE_SYMBOLS} (un-padded canonical).</li>
     * </ul>
     *
     * @param ref    The reference base content.
     * @param alt    The alternative base content.
     * @param padded Whether the variant is padded by gap symbols.
     * @return {@code true} if the variant is an insertion, {@code false} otherwise.
     */
    public static boolean isInsertion(String ref, String alt, boolean padded) {
        return isDeletion(alt, ref, padded);
    }

    /**
     * Determines whether a variant is an insertion based on its alternative content.
     * <p>
     * This method checks if the alternative base content represents an insertion. An insertion is defined as a string of at least two
     * consecutive bases from the set of valid nucleotide symbols defined in {@link Constants#BASE_SYMBOLS}.
     *
     * @param alt The alternative base content to check.
     * @return {@code true} if the alternative content represents an insertion, {@code false} otherwise.
     */
    public static boolean isInsertion(String alt) {
        return alt.length() > 1 && alt.chars().noneMatch(c -> Constants.BASE_SYMBOLS.indexOf(c) == -1);
    }

    /**
     * Determines whether a variant is a deletion, i.e.,
     * <ul>
     *     <li>either the reference base content is a string of any length of {@link Constants#BASE_SYMBOLS}
     *     and the alternative base content is a single base of {@link Constants#BASE_SYMBOLS} followed by {@link Constants#GAP}s
     *     matching the reference content's length (padded canonical),</li>
     *     <li>or the reference base content is a string of any length of {@link Constants#BASE_SYMBOLS} and the
     *     alternative content is a single base of {@link Constants#BASE_SYMBOLS} (un-padded canonical).</li>
     * </ul>
     *
     * @param ref    The reference base content.
     * @param alt    The alternative base content.
     * @param padded Whether the variant is padded by gap symbols.
     * @return {@code true} if the variant is a deletion, {@code false} otherwise.
     */
    public static boolean isDeletion(String ref, String alt, boolean padded) {
        if (padded) {
            return ref.length() == alt.length()
                    && ref.chars().noneMatch(c -> Constants.BASE_SYMBOLS.indexOf(c) == -1)
                    && isDeletion(alt);
        } else {
            return ref.length() > 1
                    && alt.length() == 1
                    && ref.chars().noneMatch(c -> Constants.BASE_SYMBOLS.indexOf(c) == -1)
                    && Constants.BASE_SYMBOLS.indexOf(alt.charAt(0)) != -1;
        }
    }

    /**
     * Determines whether a variant is a deletion based on its alternative content.
     * <p>
     * This method checks if the alternative base content represents a deletion. A deletion is defined as a string that starts with a valid
     * nucleotide base (from {@link Constants#BASE_SYMBOLS}) followed by one or more gap symbols (defined in {@link Constants#GAP}).
     *
     * @param alt The alternative base content to check.
     * @return {@code true} if the alternative content represents a deletion, {@code false} otherwise.
     */
    public static boolean isDeletion(String alt) {
        return alt.length() > 1 && Constants.BASE_SYMBOLS.indexOf(alt.charAt(0)) != -1 && alt.chars().skip(1).noneMatch(c -> c != Constants.GAP_CHAR);
    }

    /**
     * Determines whether a variant is canonical.
     * <p>
     * A variant is canonical if it is:
     * <ul>
     *     <li>a single nucleotide variant (SNV) ({@link #isSubstitution}),</li>
     *     <li>an un-padded canonical insertion ({@link #isInsertion}), or</li>
     *     <li>an un-padded canonical deletion ({@link #isDeletion}).</li>
     * </ul>
     *
     * @param ref The reference base content.
     * @param alt The alternative base content.
     * @return {@code true} if the variant is canonical, {@code false} otherwise.
     */
    public static boolean isCanonical(String ref, String alt) {
        return isSubstitution(ref, alt)
                || isInsertion(ref, alt, false)
                || isDeletion(ref, alt, false);
    }

    /**
     * Determines whether a variant is padded canonical.
     * <p>
     * A variant is padded canonical if it is:
     * <ul>
     *     <li>a single nucleotide variant (SNV) ({@link #isSubstitution}),</li>
     *     <li>a padded canonical insertion ({@link #isInsertion}), or</li>
     *     <li>a padded canonical deletion ({@link #isDeletion}).</li>
     * </ul>
     *
     * @param ref The reference base content.
     * @param alt The alternative base content.
     * @return {@code true} if the variant is padded canonical, {@code false} otherwise.
     */
    public static boolean isPaddedCanonical(String ref, String alt) {
        return isSubstitution(ref, alt)
                || isInsertion(ref, alt, true)
                || isDeletion(ref, alt, true);
    }
}
