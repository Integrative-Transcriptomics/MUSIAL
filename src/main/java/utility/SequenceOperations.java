package utility;

import datastructure.Contig;
import datastructure.Feature;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import org.apache.commons.lang3.tuple.Triple;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AmbiguityRNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.transcription.Frame;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;

import java.io.IOException;
import java.util.*;

/**
 * Utility class for performing various sequence operations.
 * <p>
 * This class provides static methods for sequence alignment, variant integration,
 * sequence translation, and other related operations. It includes methods for
 * handling nucleotide and protein sequences, as well as utilities for working
 * with gaps and variants.
 */
public final class SequenceOperations {

    /**
     * A cache for storing translated DNA sequences.
     * <p>
     * This static {@link HashMap} is used to store previously translated DNA sequences
     * to improve performance by avoiding redundant translations. The key is the hash code
     * of the translation request (including sequence and direction), and the value is the
     * translated amino acid sequence.
     */
    private final static HashMap<Integer, String> translationCache = new HashMap<>();

    /**
     * Performs global nucleotide sequence alignment using a simple scoring matrix.
     * <p>
     * This method aligns two nucleotide sequences using a gap-affine Needleman-Wunsch algorithm.
     * It utilizes a predefined scoring matrix for nucleotide matches, mismatches, and gaps.
     * <p>
     * The scoring matrix is defined as follows:
     * <ul>
     *   <li>Match: +1</li>
     *   <li>Mismatch: -1</li>
     *   <li>Gap: -1</li>
     * </ul>
     *
     * @param sequenceA        The first nucleotide sequence to align.
     * @param sequenceB        The second nucleotide sequence to align.
     * @param gapOpenPenalty   The penalty for opening a gap in the alignment.
     * @param gapExtendPenalty The penalty for extending an existing gap in the alignment.
     * @param allowGapPrefix   Whether to allow gaps at the beginning of the aligned sequences.
     * @param allowGapSuffix   Whether to allow gaps at the end of the aligned sequences.
     * @param banded           Whether to use banded alignment (restricted by the sequence's length difference).
     * @return A {@link Tuple} containing the aligned sequences.
     */
    public static Tuple<String, String> globalNucleotideSequenceAlignment(String sequenceA, String sequenceB, int gapOpenPenalty, int gapExtendPenalty,
                                                                          boolean allowGapPrefix, boolean allowGapSuffix, boolean banded) {
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
                gapOpenPenalty, gapExtendPenalty, allowGapPrefix, allowGapSuffix, banded);
    }

    /**
     * Performs global protein sequence alignment using the BLOSUM80 scoring matrix.
     * <p>
     * This method aligns two protein sequences using a gap-affine Needleman-Wunsch algorithm.
     * It utilizes the BLOSUM80 scoring matrix for amino acid matches, mismatches, and gaps.
     * <p>
     * The scoring matrix is defined as follows:
     * <ul>
     *   <li>Match: Based on BLOSUM80 values.</li>
     *   <li>Mismatch: Based on BLOSUM80 values.</li>
     *   <li>Gap penalties: Defined by the gap open and gap extend penalties.</li>
     * </ul>
     *
     * @param sequenceA        The first protein sequence to align.
     * @param sequenceB        The second protein sequence to align.
     * @param gapOpenPenalty   The penalty for opening a gap in the alignment.
     * @param gapExtendPenalty The penalty for extending an existing gap in the alignment.
     * @param allowGapPrefix   Whether to allow gaps at the beginning of the aligned sequences.
     * @param allowGapSuffix   Whether to allow gaps at the end of the aligned sequences.
     * @param banded           Whether to use banded alignment (restricted by the sequence's length difference).
     * @return A {@link Tuple} containing the aligned sequences.
     */
    public static Tuple<String, String> globalProteinSequenceAlignment(String sequenceA, String sequenceB, int gapOpenPenalty, int gapExtendPenalty,
                                                                       boolean allowGapPrefix, boolean allowGapSuffix, boolean banded) {
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
        return globalSequenceAlignment(sequenceA, sequenceB, blosum80IndexMap, blosum80, gapOpenPenalty, gapExtendPenalty, allowGapPrefix,
                allowGapSuffix, banded);
    }

    /**
     * Computes optimal pairwise global sequence alignment using a gap-affine banded Needleman-Wunsch algorithm.
     * <p>
     * This method aligns two sequences using a scoring matrix and gap penalties. It supports banded alignment
     * for performance optimization and handles marginal gaps based on the specified modes.
     * <p>
     * The alignment matrices are structured as follows:
     * <ul>
     *   <li>The y-axis corresponds to the first sequence (sequenceY).</li>
     *   <li>The x-axis corresponds to the second sequence (sequenceX).</li>
     *   <li>Insertions are traced back by walking vertically, and deletions are traced back by walking horizontally.</li>
     * </ul>
     *
     * @param sequenceA        The first sequence to align.
     * @param sequenceB        The second sequence to align.
     * @param symbolScoreIndex A mapping of characters to their respective indices in the scoring matrix.
     * @param scores           A 2D array representing the scoring matrix for character matches and mismatches.
     * @param gapOpenPenalty   The penalty for opening a gap in the alignment.
     * @param gapExtendPenalty The penalty for extending an existing gap in the alignment.
     * @param allowGapPrefix   Whether to allow gaps at the beginning of the aligned sequences.
     * @param allowGapSuffix   Whether to allow gaps at the end of the aligned sequences.
     * @param banded           Whether to use banded alignment (restricted by the sequence's length difference).
     * @return A {@link Tuple} containing the aligned sequences.
     * @throws IllegalArgumentException If the bandWidth is too narrow for the given sequences.
     */
    private static Tuple<String, String> globalSequenceAlignment(String sequenceA, String sequenceB,
                                                                 HashMap<Character, Integer> symbolScoreIndex,
                                                                 int[][] scores, int gapOpenPenalty,
                                                                 int gapExtendPenalty,
                                                                 boolean allowGapPrefix,
                                                                 boolean allowGapSuffix,
                                                                 boolean banded) {
        // Set band's width by sequence's length, if specified.
        Integer width = null;
        if (banded) {
            width = Math.abs(sequenceA.length() - sequenceB.length());
        }

        int capacity = Math.max(sequenceA.length(), sequenceB.length());
        StringBuilder alignedA = new StringBuilder(capacity);
        StringBuilder alignedB = new StringBuilder(capacity);
        char[] symbolsA = sequenceA.toCharArray();
        char[] symbolsB = sequenceB.toCharArray();

        // Initialize DP matrices.
        double[][] D = initializeMatrix(sequenceA.length() + 1, sequenceB.length() + 1);
        double[][] P = initializeMatrix(sequenceA.length() + 1, sequenceB.length() + 1);
        double[][] Q = initializeMatrix(sequenceA.length() + 1, sequenceB.length() + 1);

        // Fill dynamic programming matrices.
        P[0][0] = 0; // Insertion (in sequence A wrt. sequence B) matrix initialization.
        Q[0][0] = 0; // Deletion (in sequence A wrt. sequence B) matrix initialization.
        D[0][0] = 0; // Substitution matrix initialization.
        int gapCost;
        if (allowGapPrefix) {
            for (int i = 1; i < (Objects.isNull(width) ? sequenceA.length() + 1 : width + 1); i++) {
                gapCost = -gapOpenPenalty - i * gapExtendPenalty;
                D[i][0] = gapCost;
            }
            for (int j = 1; j < (Objects.isNull(width) ? sequenceB.length() + 1 : width + 1); j++) {
                gapCost = -gapOpenPenalty - j * gapExtendPenalty;
                D[0][j] = gapCost;
            }
        }
        for (int i = 1; i <= sequenceA.length(); i++) {
            int jLeftBound = Math.max(1, (Objects.isNull(width) ? 1 : i - width));
            int jRightBound = Math.min(sequenceB.length(), (Objects.isNull(width) ? sequenceB.length() : i + width));
            for (int j = jLeftBound; j <= jRightBound; j++) {
                P[i][j] = Math.max(
                        D[i - 1][j] - gapOpenPenalty - gapExtendPenalty,
                        P[i - 1][j] - gapExtendPenalty
                );
                Q[i][j] = Math.max(
                        D[i][j - 1] - gapOpenPenalty - gapExtendPenalty,
                        Q[i][j - 1] - gapExtendPenalty
                );
                D[i][j] = Math.max(
                        D[i - 1][j - 1] + scores[symbolScoreIndex.get(symbolsA[i - 1])][symbolScoreIndex.get(symbolsB[j - 1])],
                        Math.max(P[i][j], Q[i][j])
                );
            }
        }

        // Traceback through the matrices to construct aligned sequences.
        int i = sequenceA.length();
        int j = sequenceB.length();
        // Start from the scoreMatrix.
        String currentMatrix = "D";
        double score;
        while (i > 0 || j > 0) {

            // If the alignment is not allowed to end with a gap, enforce the first traceback step to align the last symbols.
            if (!allowGapSuffix) {
                alignedA.append(symbolsA[i - 1]);
                alignedB.append(symbolsB[j - 1]);
                i--;
                j--;
                allowGapSuffix = true;
            }

            // If banded alignment is used, ensure that the current indices are within the band's width.
            if (Objects.nonNull(width) && Math.abs(i - j) >= width) {
                if (i > 0) {
                    alignedA.append(symbolsA[i - 1]);
                    i--;
                } else {
                    alignedA.append(Constants.gapString);
                }
                if (j > 0) {
                    alignedB.append(symbolsB[j - 1]);
                    j--;
                } else {
                    alignedB.append(Constants.gapString);
                }
            }

            // Default.
            switch (currentMatrix) {
                case "D" -> {
                    score = D[i][j];
                    if (i > 0 && j > 0 && score == D[i - 1][j - 1] + scores[symbolScoreIndex.get(symbolsA[i - 1])][symbolScoreIndex.get(symbolsB[j - 1])]) {
                        alignedA.append(symbolsA[i - 1]);
                        alignedB.append(symbolsB[j - 1]);
                        i--;
                        j--;
                    } else if (i > 0 && score == P[i][j]) {
                        currentMatrix = "P";
                    } else if (i == 0) {
                        alignedA.append(Constants.gapString);
                        alignedB.append(symbolsB[j - 1]);
                        j--;
                    } else if (j > 0 && score == Q[i][j]) {
                        currentMatrix = "Q";
                    } else if (j == 0) {
                        alignedA.append(symbolsA[i - 1]);
                        alignedB.append(Constants.gapString);
                        i--;
                    }
                }
                case "P" -> {
                    score = P[i][j];
                    alignedA.append(symbolsA[i - 1]);
                    alignedB.append(Constants.gapString);
                    i--;
                    if (i > 0 && score == P[i][j] - gapExtendPenalty) {
                        currentMatrix = "P";
                    } else {
                        currentMatrix = "D";
                    }
                }
                case "Q" -> {
                    score = Q[i][j];
                    alignedA.append(Constants.gapString);
                    alignedB.append(symbolsB[j - 1]);
                    j--;
                    if (j > 0 && score == Q[i][j] - gapExtendPenalty) {
                        currentMatrix = "Q";
                    } else {
                        currentMatrix = "D";
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
     * Initializes a 2D matrix with the specified number of rows and columns.
     * <p>
     * This method creates a 2D array of doubles with the given dimensions and fills
     * each cell with {@link Double#NEGATIVE_INFINITY}.
     *
     * @param rows The number of rows in the matrix.
     * @param cols The number of columns in the matrix.
     * @return A 2D array of doubles initialized with {@link Double#NEGATIVE_INFINITY}.
     */
    private static double[][] initializeMatrix(int rows, int cols) {
        double[][] matrix = new double[rows][cols];
        Arrays.stream(matrix).forEach(row -> Arrays.fill(row, Double.NEGATIVE_INFINITY));
        return matrix;
    }

    /**
     * Pads a string with gap characters to reach a specified length.
     * <p>
     * This method appends gap characters (defined by {@link Constants#gapString})
     * to the input string until it reaches the desired length. If the input string
     * is already equal to or longer than the specified length, no padding is added.
     *
     * @param s      The input string to be padded.
     * @param length The desired length of the resulting string.
     * @return The padded string, or the original string if no padding is needed.
     */
    public static String padGaps(String s, int length) {
        return s + Constants.gapString.repeat(Math.max(0, length - s.length()));
    }

    /**
     * Removes all gap characters from the input string.
     * <p>
     * This method replaces all occurrences of the gap character (defined by {@link Constants#gapString})
     * in the input string with an empty string (defined by {@link Constants#EMPTY}).
     *
     * @param s The input string from which gaps should be removed.
     * @return A new string with all gap characters removed.
     */
    public static String stripGaps(String s) {
        return s.replaceAll(Constants.gapString, Constants.EMPTY);
    }

    /**
     * Integrates variants into a reference sequence for a given feature.
     * <p>
     * This method processes a reference sequence from a specified contig and feature, integrating
     * variants provided in a map. Variants can include single nucleotide variants (SNVs),
     * insertions, and deletions. The resulting sequence can optionally have gaps stripped.
     * <p>
     * Upstream deletions are handled by skipping affected positions and logging a warning.
     *
     * @param contig    The {@link Contig} object containing the reference sequence.
     * @param feature   The {@link Feature} object specifying the region of interest.
     * @param variants  A {@link NavigableMap} mapping positions to variant strings.
     * @param stripGaps A boolean indicating whether to remove gaps from the resulting sequence.
     * @return A {@link String} representing the updated sequence with integrated variants.
     * @throws IOException              If an error occurs while accessing the contig sequence.
     * @throws IllegalArgumentException If the contig does not have a sequence or the feature is incompatible.
     */
    public static String integrateVariants(Contig contig, Feature feature, NavigableMap<Integer, String> variants, boolean stripGaps) throws IOException {
        // Validate contig and feature compatibility.
        if (!contig.hasSequence()) {
            throw new IllegalArgumentException("Contig %s does not have a sequence.".formatted(contig.name));
        }
        if (!Objects.equals(feature.contig, contig.name)) {
            throw new IllegalArgumentException("Contig %s is not the parent of feature %s.".formatted(feature.name, contig.name));
        }

        // Initialize variables for processing.
        char[] referenceChars = contig.getSubsequence(feature.start, feature.end).toCharArray();
        StringBuilder result = new StringBuilder(referenceChars.length);
        int deletionCount = 0;

        // Iterate through the reference sequence positions.
        for (int pos = feature.start, idx = 0; idx < referenceChars.length; pos++, idx++) {
            if (variants.containsKey(pos)) {
                String variant = variants.get(pos);

                // Handle upstream deletions.
                if (deletionCount > 0) {
                    result.append(Constants.gapString);
                    deletionCount--;
                    Logging.logWarning("Skip variant %s at position %d due to upstream deletion.".formatted(variant, pos));
                    continue;
                }

                // Process variant types.
                switch (contig.getVariantInformation(pos, variant).type) {
                    case SNV, INSERTION -> result.append(variant);
                    case DELETION -> {
                        result.append(variant.charAt(0));
                        deletionCount += variant.length() - 1;
                    }
                }
            } else {
                // Handle gaps from deletions or append reference character.
                if (deletionCount > 0) {
                    result.append(Constants.gapString);
                    deletionCount--;
                } else {
                    result.append(referenceChars[idx]);
                }
            }
        }

        // Return the final sequence, optionally stripping gaps.
        return stripGaps ? stripGaps(result.toString()) : result.toString();
    }

    /**
     * Translates a DNA sequence into an amino-acid sequence. The translation is always performed in the 1-frame.
     * Utilizes the <a href="https://github.com/biojava/biojava">BioJava library</a> for translation.
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
            // Define ambiguity compound sets. See: https://github.com/biojava/biojava-tutorial/blob/master/core/translating.md
            AmbiguityDNACompoundSet ambiguityDNACompoundSet = AmbiguityDNACompoundSet.getDNACompoundSet();
            CompoundSet<NucleotideCompound> nucleotideCompoundSet = AmbiguityRNACompoundSet.getRNACompoundSet();
            // Initialize the transcription engine. See: https://github.com/biojava/biojava-tutorial/blob/master/core/translating.md
            TranscriptionEngine engine = new
                    TranscriptionEngine.Builder().dnaCompounds(ambiguityDNACompoundSet).rnaCompounds(nucleotideCompoundSet).build();
            // Initialize the DNA sequence.
            Sequence<NucleotideCompound> dna = new DNASequence(sequence);
            String translatedSequence;
            if (reverse)
                translatedSequence = engine.multipleFrameTranslation(dna, Frame.REVERSED_ONE).get(Frame.REVERSED_ONE).getSequenceAsString();
            else
                translatedSequence = engine.multipleFrameTranslation(dna, Frame.ONE).get(Frame.ONE).getSequenceAsString();
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
     * reference-, and variant content. The relative position is the 0-based position of the variant in the reference sequence without gaps.
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
            } else if (referenceChars[i] == Constants.gapChar) { // Insertion:
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
                    referenceBuilder.append(Constants.gapChar);
                    alternativeBuilder.append(alternativeChars[i]);
                    isSubstitution = false;
                    isInsertion = true;
                }
                noInsertions++;
            } else if (alternativeChars[i] == Constants.gapChar) { // Deletion
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
                    alternativeBuilder.append(Constants.gapChar);
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
}
