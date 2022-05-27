package utility;

import com.google.common.base.Splitter;
import datastructure.FeatureEntry;
import datastructure.VariantsDictionary;
import exceptions.MusialBioException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.javatuples.Triplet;

/**
 * Comprises static methods used to solve specifically bioinformatics problems.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class Bio {

  /**
   * Enum to store different modes to handle prefix gaps for global sequence alignment.
   */
  public enum GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES {
    FREE, PENALIZE, FORBID
  }

  /**
   * Character to indicate a gap/missing information regarding a nucleotide and amino acid sequence.
   */
  public static final char NO_MATCH_CHAR = '#';
  /**
   * Character used as a separator in padded sequences.
   */
  public static final char SEPARATOR_CHAR = '/';
  /**
   * One-letter code character used to indicate a translated stop codon.
   */
  public static final char TERMINATION_AA1 = '*';
  /**
   * Three-letter code String used to indicate a translated stop codon.
   */
  public static final String TERMINATION_AA3 = "TER";
  /**
   * One-letter code character used to indicate the translation of an incomplete codon.
   */
  public static final char ANY_AA1 = 'X';
  /**
   * Three-letter code string used to indicate the translation of an incomplete codon.
   */
  public static final String ANY_AA3 = "ANY";
  /**
   * TODO
   */
  public static final char DELETION_AA1 = '-';
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
    put(String.valueOf(ANY_AA1), ANY_AA3);
    put(String.valueOf(TERMINATION_AA1), TERMINATION_AA3);
  }};
  /**
   * Hash map mapping amino-acid one letter symbols to indices used to access substitution matrix values.
   */
  public static final HashMap<Character, Integer> AA1_PAM120_INDEX = new HashMap<>() {{
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
    put(ANY_AA1, 20);
    put(TERMINATION_AA1, 21);
  }};
  /**
   * The PAM120 matrix used for alignment computation.
   * <p>
   * Includes values for 'X'; any amino-acid and '*'; termination.
   */
  private static final int[][] PAM120 = {
      {1, -1, 0, 0, -1, 0, 0, 0, -1, 0, -1, -1, -1, -1, 0, 0, 0, -2, -1, 0, -1, -4},
      {-1, 2, 0, -1, -1, 0, -1, -1, 0, -1, -1, 1, 0, -2, 0, 0, -1, 0, -2, -1, -1, -4},
      {0, 0, 1, 1, -2, 0, 0, 0, 1, -1, -1, 0, -1, -1, -1, 0, 0, -2, -1, -1, -1, -4},
      {0, -1, 1, 2, -2, 0, 1, 0, 0, -1, -2, 0, -1, -2, -1, 0, 0, -3, -2, -1, -1, -4},
      {-1, -1, -2, -2, 3, -2, -2, -2, -1, -1, -3, -2, -2, -2, -1, 0, -1, -3, 0, -1, -1, -4},
      {0, 0, 0, 0, -2, 2, 1, -1, 1, -1, -1, 0, 0, -2, 0, -1, -1, -2, -2, -1, -1, -4},
      {0, -1, 0, 1, -2, 1, 2, 0, 0, -1, -2, 0, -1, -2, -1, 0, -1, -3, -2, -1, -1, -4},
      {0, -1, 0, 0, -2, -1, 0, 2, -1, -1, -2, -1, -1, -2, -1, 0, 0, -3, -2, -1, -1, -4},
      {-1, 0, 1, 0, -1, 1, 0, -1, 2, -1, -1, -1, -1, -1, 0, -1, -1, -1, 0, -1, -1, -4},
      {0, -1, -1, -1, -1, -1, -1, -1, -1, 2, 0, -1, 1, 0, -1, -1, 0, -2, -1, 1, -1, -4},
      {-1, -1, -1, -2, -3, -1, -2, -2, -1, 0, 2, -1, 1, 0, -1, -1, -1, -1, -1, 0, -1, -4},
      {-1, 1, 0, 0, -2, 0, 0, -1, -1, -1, -1, 2, 0, -2, -1, 0, 0, -2, -2, -1, -1, -4},
      {-1, 0, -1, -1, -2, 0, -1, -1, -1, 1, 1, 0, 3, 0, -1, -1, 0, -2, -1, 0, -1, -4},
      {-1, -2, -1, -2, -2, -2, -2, -2, -1, 0, 0, -2, 0, 3, -2, -1, -1, 0, 2, -1, -1, -4},
      {0, 0, -1, -1, -1, 0, -1, -1, 0, -1, -1, -1, -1, -2, 2, 0, 0, -2, -2, -1, -1, -4},
      {0, 0, 0, 0, 0, -1, 0, 0, -1, -1, -1, 0, -1, -1, 0, 1, 1, -1, -1, -1, -1, -4},
      {0, -1, 0, 0, -1, -1, -1, 0, -1, 0, -1, 0, 0, -1, 0, 1, 1, -2, -1, 0, -1, -4},
      {-2, 0, -2, -3, -3, -2, -3, -3, -1, -2, -1, -2, -2, 0, -2, -1, -2, 4, -1, -3, -1, -4},
      {-1, -2, -1, -2, 0, -2, -2, -2, 0, -1, -1, -2, -1, 2, -2, -1, -1, -1, 3, -1, -1, -4},
      {0, -1, -1, -1, -1, -1, -1, -1, -1, 1, 0, -1, 0, -1, -1, -1, 0, -3, -1, 2, -1, -4},
      {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -4},
      {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 0}
  };

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
  public static String translateCodon(String codon, boolean asAA3, boolean includeTermination,
                                      boolean includeIncomplete)
      throws MusialBioException {
    if (codon.length() != 3) {
      if (!includeIncomplete) {
        throw new MusialBioException("Unable to translate codon " + codon + " with length different from three.");
      } else {
        if (asAA3) {
          return ANY_AA3;
        } else {
          return String.valueOf(ANY_AA1);
        }
      }
    } else if (codon.contains("N")) { // TODO: Add constant definition for nucleotide symbols.
      if (asAA3) {
        return ANY_AA3;
      } else {
        return String.valueOf(ANY_AA1);
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
   * TODO
   *
   * @param s1
   * @param s2
   * @return
   */
  public static double computeJaccardIndex(String s1, String s2) {
    return computeJaccardIndex(s1, s2, 21);
  }

  /**
   * TODO
   *
   * @param s1
   * @param s2
   * @param k
   * @return
   */
  public static double computeJaccardIndex(String s1, String s2, int k) {
    if (s1.length() < k || s2.length() < k) {
      return 0.0F;
    }
    HashSet<String> s1Kmers = new HashSet<>();
    HashSet<String> s2Kmers = new HashSet<>();
    for (int i = 0; i + k <= s1.length(); i++) {
      s1Kmers.add(s1.substring(i, i + k));
    }
    for (int i = 0; i + k <= s2.length(); i++) {
      s2Kmers.add(s2.substring(i, i + k));
    }
    HashSet<String> intersection = new HashSet<>(s1Kmers);
    intersection.retainAll(s2Kmers);
    HashSet<String> union = new HashSet<>(s1Kmers);
    union.addAll(s2Kmers);
    return (double) intersection.size() / (double) union.size();
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
  public static String translateNucSequence(Iterable<String> splitNucSequence, boolean includeTermination,
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
   * @throws MusialBioException If any codon with a length different than three is detected.
   */
  public static String translateNucSequence(String nucSequence, boolean includeTermination,
                                            boolean includeIncomplete, boolean asSense) throws MusialBioException {
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
   * @throws MusialBioException If any base occurs in the sequence that can not be reversed, i.e. for which no
   *                            complement base is known.
   */
  public static String reverseComplement(String sequence) throws MusialBioException {
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
   * @throws MusialBioException If no complement base for the specified base is known.
   */
  public static Character invertBase(char base) throws MusialBioException {
    return switch (base) {
      case 'A' -> 'T';
      case 'C' -> 'G';
      case 'G' -> 'C';
      case 'T' -> 'A';
      default -> base;
    };
  }

  /**
   * Removes all gap symbols ("-") from the passed sequence.
   *
   * @param sequence The sequence from which gaps should be removed.
   * @return {@link String}, the passed sequence without any gap symbols.
   */
  public static String removeGaps(String sequence) {
    return sequence.replace("-", "");
  }

  /**
   * TODO
   *
   * @param aaSeq1
   * @param aaSeq2
   * @return
   */
  public static Triplet<Integer, String, String> globalNucleotideSequenceAlignment(String nucSeq1, String nucSeq2,
                                                                                   GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES left_mode,
                                                                                   GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES right_mode) {
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
        {-1, -1, -1, -1, 1},
    };
    return globalSequenceAlignment(nucSeq1, nucSeq2, simpleNucleotideScoringMatrixIndexMap,
        simpleNucleotideScoringMatrix,
        2, 1, left_mode, right_mode);
  }

  /**
   * TODO
   *
   * @param nucSeq1
   * @param nucSeq2
   * @return
   */
  public static Triplet<Integer, String, String> globalAminoAcidSequenceAlignment(String aaSeq1, String aaSeq2,
                                                                                  GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES left_mode,
                                                                                  GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES right_mode) {
    return globalSequenceAlignment(aaSeq1, aaSeq2, AA1_PAM120_INDEX, PAM120, 5, 4, left_mode, right_mode);
  }

  /**
   * FIXME
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
  public static Triplet<Integer, String, String> globalSequenceAlignment(String seq1, String seq2,
                                                                         HashMap<Character, Integer> scoringMatrixIndexMap,
                                                                         int[][] scoringMatrix, int gapOpenPenalty,
                                                                         int gapExtendPenalty,
                                                                         GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES left_mode,
                                                                         GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES right_mode) {
    /*
    Notes on the alignment matrices:
    - The y-axis will yield the proteins amino acid sequence.
    - The x-axis will yield the translated reference gene amino acid sequence.
    - Indels are wrt. the proteins amino acid sequence; thus insertions are traced back by walking vertically and
    deletions are traced back by walking horizontally in the matrix.
     */
    int[][] alignmentMatrix = new int[seq1.length() + 1][seq2.length() + 1];
    int[][] matchScoreMatrix = new int[seq1.length() + 1][seq2.length() + 1];
    int[][] insertionScoreMatrix = new int[seq1.length() + 1][seq2.length() + 1];
    int[][] deletionScoreMatrix = new int[seq1.length() + 1][seq2.length() + 1];
    char[][] tracebackMatrix = new char[seq1.length() + 1][seq2.length() + 1];
    char[] aaSeq1Array = seq1.toCharArray();
    char[] aaSeq2Array = seq2.toCharArray();
    int alignmentScore;
    StringBuilder aaSeq1Builder = new StringBuilder();
    StringBuilder aaSeq2Builder = new StringBuilder();
    /*
    (1) Compute global sequence alignment.
     */
    alignmentMatrix[0][0] = 0;
    matchScoreMatrix[0][0] = 0;
    insertionScoreMatrix[0][0] = 0;
    deletionScoreMatrix[0][0] = 0;
    int gapCost;
    // i -> PREFIX
    for (int i = 0; i < seq1.length() + 1; i++) {
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
    for (int j = 0; j < seq2.length() + 1; j++) {
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
      for (int j = 1; j < seq2.length() + 1; j++) {
        matchScoreMatrix[i][j] =
            alignmentMatrix[i - 1][j - 1]
                + scoringMatrix[scoringMatrixIndexMap.get(aaSeq1Array[i - 1])][scoringMatrixIndexMap
                .get(aaSeq2Array[j - 1])];
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
    alignmentScore = alignmentMatrix[seq1.length()][seq2.length()];
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
      }
    }
    Collections.reverse(tracebackPath);
    /*
    (3) Iterate over traceback path and construct aligned nucleotide and amino-acid sequence.
    */
    int aaSeq1Index = 0;
    int aaSeq2Index = 0;
    for (Character character : tracebackPath) {
      // Start walking along the traceback path.
      tracebackDirection = character;
      if (tracebackDirection == 'M') {
        // CASE: Match or mismatch of nucleotide and amino-acid sequence.
        aaSeq1Builder.append(aaSeq1Array[aaSeq1Index]);
        aaSeq1Index += 1;
        aaSeq2Builder.append(aaSeq2Array[aaSeq2Index]);
        aaSeq2Index += 1;
      } else if (tracebackDirection == 'D') {
        // CASE: Deletion wrt. to first amino acid sequence.
        aaSeq1Builder.append(DELETION_AA1);
        aaSeq2Builder.append(aaSeq2Array[aaSeq2Index]);
        aaSeq2Index += 1;
      } else if (tracebackDirection == 'I') {
        // CASE: Insertion wrt. to first amino acid sequence.
        aaSeq1Builder.append(aaSeq1Array[aaSeq1Index]);
        aaSeq1Index += 1;
        aaSeq2Builder.append(DELETION_AA1);
      }
    }
    /*
    (6) Insert results into Triplet.
     */
    return new Triplet<>(alignmentScore, aaSeq1Builder.toString(), aaSeq2Builder.toString());
  }

  public static ArrayList<String> getVariantsOfAlignedSequences(String referenceSequence,
                                                                String variantSequence,
                                                                String variantStringPrefix) {
    ArrayList<String> variants = new ArrayList<>();
    StringBuilder variantBuilder = new StringBuilder();
    StringBuilder referenceBuilder = new StringBuilder();
    char[] referenceSequenceArray = referenceSequence.toCharArray();
    char[] variantSequenceArray = variantSequence.toCharArray();
    char referenceContent;
    char variantContent;
    int variantStart = 1;
    int alignmentLength = referenceSequence.length();
    boolean isSubstitution = false;
    boolean isInsertion = false;
    boolean isDeletion = false;
    boolean referenceHasPrefixGap = true;
    for (int i = 0; i < alignmentLength; i++) {
      referenceContent = referenceSequenceArray[i];
      variantContent = variantSequenceArray[i];
      if (referenceHasPrefixGap && (referenceContent != Bio.DELETION_AA1)) {
        referenceHasPrefixGap = false;
      }
      if (referenceContent == variantContent) {
        // CASE: Match.
        if (isSubstitution || isInsertion || isDeletion) {
          variants.add(variantStringPrefix + "@" + variantStart + "@" + variantBuilder + "@" + referenceBuilder);
          variantBuilder.setLength(0);
          referenceBuilder.setLength(0);
        }
        isSubstitution = false;
        isInsertion = false;
        isDeletion = false;
      } else if (referenceContent == Bio.DELETION_AA1) {
        // CASE: Insertion (in variant).
        if (isSubstitution || isDeletion) {
          variants.add(variantStringPrefix + "@" + variantStart + "@" + variantBuilder + "@" + referenceBuilder);
          variantBuilder.setLength(0);
          referenceBuilder.setLength(0);
        }
        if (!isInsertion) {
          if (referenceHasPrefixGap) {
            variantStart = 0;
          } else {
            variantBuilder.append(referenceSequenceArray[i - 1]);
            variantStart = i;
            referenceBuilder.append(referenceSequenceArray[i - 1]);
          }
        }
        variantBuilder.append(variantContent);
        isSubstitution = false;
        isInsertion = true;
        isDeletion = false;
      } else if (variantContent == Bio.DELETION_AA1) {
        // CASE: Deletion (in variant).
        if (isSubstitution || isInsertion) {
          variants.add(variantStringPrefix + "@" + variantStart + "@" + variantBuilder + "@" + referenceBuilder);
          variantBuilder.setLength(0);
          referenceBuilder.setLength(0);
        }
        if (!isDeletion) {
          // Alignment is forbidden to start with a gap.
          variantBuilder.append(referenceSequenceArray[i - 1]);
          variantStart = i;
          referenceBuilder.append(referenceSequenceArray[i - 1]);
        }
        referenceBuilder.append(referenceContent);
        variantBuilder.append(variantContent);
        isSubstitution = false;
        isInsertion = false;
        isDeletion = true;
      } else {
        // CASE: Mismatch/Substitution.
        if (isDeletion || isInsertion) {
          variants.add(variantStringPrefix + "@" + variantStart + "@" + variantBuilder + "@" + referenceBuilder);
          variantBuilder.setLength(0);
          referenceBuilder.setLength(0);
        }
        if (!isSubstitution) {
          variantStart = i + 1;
        }
        variantBuilder.append(variantContent);
        referenceBuilder.append(referenceContent);
        isSubstitution = true;
        isInsertion = false;
        isDeletion = false;
      }
    }
    if (isSubstitution || isInsertion || isDeletion) {
      variants.add(variantStringPrefix + "@" + variantStart + "@" + variantBuilder + "@" + referenceBuilder);
      variantBuilder.setLength(0);
      referenceBuilder.setLength(0);
    }
    return variants;
  }


  /**
   * @param referenceAllele
   * @param alternateAllele
   * @return
   */
  public static ArrayList<String> resolveAmbiguousVariant(String referenceAllele,
                                                          String alternateAllele) {
    Triplet<Integer, String, String> alignedAlleles =
        globalNucleotideSequenceAlignment(referenceAllele, alternateAllele,
            GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FORBID, GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE);
    return getVariantsOfAlignedSequences(alignedAlleles.getValue1(), alignedAlleles.getValue2(),
        "");
  }

  /**
   * TODO
   *
   * @param sequence
   * @return
   */
  public static String getSequenceIdentifier(String sequence) {
    HashMap<String, String> m = new HashMap<>() {{
      put("-", "X");
      put("0", "A");
      put("1", "D");
      put("2", "E");
      put("3", "F");
      put("4", "I");
      put("5", "G");
      put("6", "O");
      put("7", "U");
      put("8", "J");
      put("9", "K");
    }};
    return Stream.of(String.valueOf(sequence.hashCode())).map(m::get).collect(Collectors.joining());
  }

  /**
   * TODO
   *
   * @param sequence
   * @return
   */
  public static String getSequenceIdentifier(String sequence, String reference) {
    String sequenceIdentifier = getSequenceIdentifier(sequence);
    double jaccardIndex = computeJaccardIndex(sequence, reference);
    return sequenceIdentifier + "x" + jaccardIndex;
  }

  public static ArrayList<Triplet<String, String, ArrayList<String>>> inferProteoform(
      VariantsDictionary variantsDictionary, String fId, String sId)
      throws MusialBioException {
    if (variantsDictionary.features.get(fId).allocatedProtein != null) {
      String sampleVSwab = variantsDictionary.samples.get(sId).annotations.get(fId + "vSwab");
      if (sampleVSwab != null) {
        int variantPosition;
        int relativeVariantPosition;
        String variantContent;
        String relativeVariantContent;
        String effectAnnotation;
        FeatureEntry featureEntry = variantsDictionary.features.get(fId);
        String referenceSequence = featureEntry.nucleotideSequence;
        if (!featureEntry.isSense) {
          referenceSequence = Bio.reverseComplement(referenceSequence);
        }
        /*
        1. Part: Infer relative variant position and content, i.e. wrt. feature start and coding direction.
         */
        String[] sampleVariants = sampleVSwab.split("\\|");
        HashMap<Integer, String> sampleVariantsMap = new HashMap<>();
        for (String sV : sampleVariants) {
          variantContent = sV.split("@")[0];
          variantPosition = Integer.parseInt(sV.split("@")[1]);
          relativeVariantPosition = featureEntry.isSense ? variantPosition - featureEntry.start + 1 :
              featureEntry.end - variantPosition + 1;
          effectAnnotation =
              variantsDictionary.variants.get(variantPosition).get(variantContent).annotations
                  .get(featureEntry.annotations.get("geneName") +
                      "_EFF");
          if (effectAnnotation == null) {
            Logging.logWarning("Skip variant " + sV + " due to missing effect annotation.");
            continue;
          }
          if (!featureEntry.isSense) {
            // In the case of anti-sense coding the variant is converted/shifted relative to coding direction.
            if (variantContent.length() == 1) {
              // CASE: Single residue substitution.
              relativeVariantContent = String.valueOf(Bio.invertBase(variantContent.charAt(0)));
            } else {
              if (variantContent.contains("-")) {
                // CASE: Deletion.
                relativeVariantContent = variantContent.substring(1);
                relativeVariantPosition -= relativeVariantContent.length() + 1;
                relativeVariantContent =
                    referenceSequence.charAt(relativeVariantPosition - 1) + relativeVariantContent;
              } else {
                // CASE: Insertion.
                relativeVariantContent = Bio.reverseComplement(variantContent.substring(1));
                relativeVariantPosition -= relativeVariantContent.length() + 1;
                relativeVariantContent =
                    referenceSequence.charAt(relativeVariantPosition - 1) + relativeVariantContent;
              }
            }
          } else {
            relativeVariantContent = variantContent;
          }
          sampleVariantsMap.put(
              relativeVariantPosition,
              relativeVariantContent + "|" + variantContent + "@" +
                  variantPosition + "|" + effectAnnotation.split("\\|")[0]
          );
        }
        /*
        2. Part: Iterate (codon wise) over reference sequence and infer variable segments.
         */
        int frameImbalance = 0;
        int skipPositions = 0;
        int variableSegmentStart = 0;
        int insertedPositions = 0;
        boolean isVariableSegment = false;
        String referenceAminoacidSequence = Bio.translateNucSequence(referenceSequence, true, true, true);
        String variableSegmentSequence = "";
        StringBuilder variableSegmentNucleotides = new StringBuilder();
        StringBuilder variableSegmentAminoacids = new StringBuilder();
        ArrayList<Triplet<String, String, ArrayList<String>>> proteoformVariableSegments = new ArrayList<>();
        ArrayList<String> affectingVariants = new ArrayList<>();
        for (int cdn = 0; cdn < referenceSequence.length() / 3; cdn++) {
          for (int rpos = 3 * (cdn + 1) - 2; rpos <= 3 * (cdn + 1); rpos++) {
            if (skipPositions > 0) {
              skipPositions -= 1;
              continue;
            }
            if (sampleVariantsMap.containsKey(rpos)) {
              isVariableSegment = isVariableSegment || !sampleVariantsMap.get(rpos).split("\\|")[2].contains(
                  "synonymous");
              variableSegmentStart = cdn + 1;
              relativeVariantContent = sampleVariantsMap.get(rpos).split("\\|")[0];
              variantContent = sampleVariantsMap.get(rpos).split("\\|")[1];
              if (relativeVariantContent.contains("-")) {
                frameImbalance -= relativeVariantContent.chars().filter(c -> c == '-').count();
                skipPositions += relativeVariantContent.chars().filter(c -> c == '-').count();
              } else {
                frameImbalance += relativeVariantContent.length() - 1;
              }
              relativeVariantContent = relativeVariantContent.replace("-", "");
              variableSegmentNucleotides.append(relativeVariantContent);
              affectingVariants.add(variantContent);
            } else {
              variableSegmentNucleotides.append(referenceSequence.charAt(rpos - 1));
            }
          }
          if (isVariableSegment) {
            if (variableSegmentNucleotides.length() % 3 == 0 && frameImbalance % 3 == 0) {
              // CASE: Frameshift resolved or no frameshift.
              variableSegmentAminoacids
                  .append(Bio.translateNucSequence(variableSegmentNucleotides.toString(), true, true,
                      true));
              // Determine position suffix, i.e. number of added or deleted positions.
              String positionSuffix = "";
              if (frameImbalance < 0) {
                try {
                  int excessLength =
                      (variableSegmentStart + variableSegmentAminoacids.length() - 1) -
                          referenceAminoacidSequence.length();
                  if (excessLength > 0) {
                    positionSuffix = "+" + excessLength;
                    variableSegmentSequence = variableSegmentAminoacids.toString();
                  } else {
                    frameImbalance = frameImbalance * -1;
                    positionSuffix = "-" + (frameImbalance / 3);
                    variableSegmentSequence = Bio.globalAminoAcidSequenceAlignment(
                        referenceAminoacidSequence.substring(
                            variableSegmentStart - 1,
                            variableSegmentStart + variableSegmentAminoacids.length() + (frameImbalance / 3) - 1
                        ),
                        variableSegmentAminoacids.toString(),
                        GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE,
                        GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE
                    ).getValue2();
                  }
                } catch (Exception e) {
                  System.out.println(referenceAminoacidSequence + IO.LINE_SEPARATOR);
                  System.out
                      .println("|".repeat(variableSegmentStart - 1) + variableSegmentAminoacids + IO.LINE_SEPARATOR);
                  System.out.println(variableSegmentStart);
                  System.out.println(frameImbalance);
                  throw e;
                }
              } else {
                positionSuffix = "+" + (frameImbalance / 3);
                insertedPositions += (frameImbalance / 3);
                variableSegmentSequence = variableSegmentAminoacids.toString();
              }
              proteoformVariableSegments.add(new Triplet<>(
                  variableSegmentStart + positionSuffix,
                  variableSegmentSequence,
                  affectingVariants
              ));
              // Reset all variables.
              frameImbalance = 0;
              variableSegmentStart = 0;
              isVariableSegment = false;
              variableSegmentNucleotides.setLength(0);
              variableSegmentAminoacids.setLength(0);
              affectingVariants = new ArrayList<>();
            }
          } else {
            variableSegmentNucleotides.setLength(0);
          }
        }
        // If after iterating over all positions the variable segment has not ended.
        if (variableSegmentNucleotides.length() != 0) {
          // Pad with Ns in order to complete last codon.
          int padNs = 0;
          if (variableSegmentNucleotides.length() % 3 != 0) {
            padNs = (3 - variableSegmentNucleotides.length() % 3);
            frameImbalance += padNs;
            variableSegmentNucleotides.append("N".repeat(padNs));
          }
          variableSegmentAminoacids
              .append(Bio.translateNucSequence(variableSegmentNucleotides.toString(), true, true,
                  true));
          // Determine position suffix, i.e. number of added or deleted positions.
          String positionSuffix = "";
          if (frameImbalance < 0) {
            // FIXME: Inv. bug with end segment being too long while frame imbalance is negative.
            int excessLength =
                (variableSegmentStart + variableSegmentAminoacids.length() - 1) - referenceAminoacidSequence.length();
            if (excessLength > 0) {
              positionSuffix = "+" + excessLength;
              variableSegmentSequence = variableSegmentAminoacids.toString();
            } else {
              frameImbalance = frameImbalance * -1;
              positionSuffix = "-" + (frameImbalance / 3);
              variableSegmentSequence = Bio.globalAminoAcidSequenceAlignment(
                  referenceAminoacidSequence.substring(
                      variableSegmentStart - 1,
                      variableSegmentStart + variableSegmentAminoacids.length() + (frameImbalance / 3) - 2
                  ),
                  variableSegmentAminoacids.toString(),
                  GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE,
                  GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE
              ).getValue2();
            }
          } else {
            positionSuffix = "+" + (frameImbalance / 3);
            variableSegmentSequence = variableSegmentAminoacids.toString();
          }
          proteoformVariableSegments.add(new Triplet<>(
              variableSegmentStart + positionSuffix,
              variableSegmentSequence,
              affectingVariants
          ));
        }
        return proteoformVariableSegments;
      }
    }
    return new ArrayList<>();
  }
}
