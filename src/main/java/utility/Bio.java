package utility;

import com.google.common.base.Splitter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;

/**
 * Comprises static methods that implement algorithms to solve specifically bioinformatics problems.
 * <p>
 * Currently only methods for sequence alignment exist which are used to align a nucleotide with a amino acid
 * sequence. This is done in a special way by first translating the nucleotide sequence into an amino acid sequence.
 * The alignment information is used to transform the translated nucleotide sequence back into a nucleotide sequence.
 * This process is applied enable a visualization of the respective sequences and should not be used for general
 * local sequence alignment.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class Bio {

  /**
   * Character to indicate a gap/missing information regarding a nucleotide and amino acid sequence.
   */
  public static final char NO_MATCH_CHAR = '?';
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
   * The BLOSUM62 substitution matrix.
   */
  private static final SubstitutionMatrix<AminoAcidCompound> BLOSUM62 = SimpleSubstitutionMatrix.getBlosum62();

  /**
   * Matches a amino-acid and nucleotide sequence.
   * <p>
   * This is done by translating the nucleotide sequence, for every possible frameshift if applicable, into a amino
   * acid sequence and aligning the translated nucleotide sequence with the amino-acid sequence. The best alignment,
   * here the one with the least number of gaps, is used as the best matching sequences.
   *
   * @param aaSequence  {@link String} representing a amino-acid sequence.
   * @param nucSequence {@link String} representing a nucleotide sequence.
   * @return {@link ArrayList<String>} comprising the padded amino acid (at index 0) and nucleotide sequence (at
   * index 1) that were derived as the best alignment.
   * @throws CompoundNotFoundException If a compound, i.e. amino acid, could not be translated.
   */
  public static ArrayList<String> matchSequences(String aaSequence, String nucSequence)
      throws CompoundNotFoundException {
    // Initialize hash map to store matched sequences.
    HashMap<Integer, ArrayList<String>> alignedSequences = new HashMap<>();
    // Remove, but store index of, all gap symbols in the nucleotide sequence (this would rise errors in the sequence
    // alignment process).
    HashSet<Integer> gapsFromSamples = getGapsFromSequence(nucSequence);
    nucSequence = removeGapsFromSequence(nucSequence);
    // Check if nucleotide sequence is fully coding sequence, i.e. divisible by three.
    int residual = nucSequence.length() % 3;
    // Based on the residual assign frame shift and split nucleotide sequence.
    String translatedNucSequence;
    ArrayList<String> splitNucSequence = new ArrayList<>();
    if (residual == 0) {
      // Truncate no characters.
      for (String s : Splitter.fixedLength(3).split(nucSequence)) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, "", "", alignedSequences, gapsFromSamples);
    } else if (residual == 1) {
      // Truncate character at start.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(1))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, nucSequence.substring(0, 1), "",
          alignedSequences, gapsFromSamples);
      // Truncate character at end.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(0, nucSequence.length() - 1))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, "",
          nucSequence.substring(nucSequence.length() - 1), alignedSequences, gapsFromSamples);
    } else {
      // Truncate two characters at start.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(2))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, nucSequence.substring(0, 2), "",
          alignedSequences, gapsFromSamples);
      // Truncate two characters at end.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(0, nucSequence.length() - 2))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, "",
          nucSequence.substring(nucSequence.length() - 2), alignedSequences, gapsFromSamples);
      // Truncate one character at start and end.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(1, nucSequence.length() - 1))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, nucSequence.substring(0, 1),
          nucSequence.substring(nucSequence.length() - 1), alignedSequences, gapsFromSamples);
    }
    // The sequences with the fewest number of mismatches characters, i.e. alignment gaps, are chosen as the correct
    // ones.
    int minScore = Integer.MAX_VALUE;
    ArrayList<String> bestAlignment = new ArrayList<>();
    for (Integer score : alignedSequences.keySet()) {
      if (score < minScore) {
        bestAlignment = alignedSequences.get(score);
      }
      minScore = score;
    }
    return bestAlignment;
  }

  /**
   * Translates a nucleotide sequence, split into codons of length three, into a single-letter amino acid sequence.
   *
   * @param splitNucSequence {@link ArrayList<String>} representing a nucleotide sequence that was split into codons
   *                         of length 3.
   * @return {@link String} representing the translated nucleotide sequence.
   */
  private static String translateNucSequence(ArrayList<String> splitNucSequence) {
    StringBuilder translatedNucSequenceBuilder = new StringBuilder();
    for (String s : splitNucSequence) {
      translatedNucSequenceBuilder.append(CODON_MAP.get(s));
    }
    return translatedNucSequenceBuilder.toString();
  }

  /**
   * Conducts a local alignment of amino acid sequences using the BLOSUM62 scoring matrix.
   * <p>
   * Thereby, one amino acid sequence is derived from a nucleotide sequence.
   * <p>
   * Alignment gaps are replaced by special characters to separate alignment gaps from such gaps that are
   * introduces to the passed nucleotide sequence by insertions of samples. In addition, the non-nucleotide-derived
   * amino acid sequence is padded with symbols to represent a truncation of the original nucleotide sequence that
   * may be present from a possible frameshift.
   * <p>
   * The aligned translated nucleotide sequence is mapped back to a nucleotide sequence and, both, the nucleotide
   * sequence and amino acid sequence are stored in a map accessible via the number of gaps in the alignment.
   *
   * @param aaSequence              {@link String} representing a non-nucleotide-derived amino acid sequence.
   * @param translatedAASequence    {@link String} representing a amino acid sequence that originates from translating a
   *                                nucleotide sequence.
   * @param splitNucleotideSequence {@link String} the original nucleotide sequence split into chunks of length three, i.e.
   *                                codons.
   * @param truncatedStart          {@link String} nucleotides that were truncated from the start of the nucleotide sequence
   *                                due to a frameshift.
   * @param truncatedEnd            {@link String} nucleotides that were truncated from the end of the nucleotide sequence due
   *                                to a frameshift.
   * @param resultsMap              {@link HashMap} mapping the alignment score, here the number of gaps, to the aligned sequences.
   * @param gapPositions            {@link HashSet<Integer>} containing positions at which a gap was removed from the
   *                                original nucleotide sequence.
   * @throws CompoundNotFoundException If a compound, i.e. amino acid, could not be translated.
   */
  private static void alignSequences(String aaSequence, String translatedAASequence,
                                     ArrayList<String> splitNucleotideSequence, String truncatedStart,
                                     String truncatedEnd, HashMap<Integer, ArrayList<String>> resultsMap,
                                     HashSet<Integer> gapPositions)
      throws CompoundNotFoundException {
    StringBuilder paddedNucSequenceBuilder = new StringBuilder();
    StringBuilder paddedAASequenceBuilder = new StringBuilder();
    StringBuilder paddedAASequenceAnnotations = new StringBuilder();
    char[] alignedAASequence;
    char[] alignedTranslatedAASequence;
    SequencePair<ProteinSequence, AminoAcidCompound> alignedSequences;
    if (aaSequence.length() > translatedAASequence.length()) {
      // Align translated amino acid sequence against non-translated amino acid sequence.
      alignedSequences = Alignments.getPairwiseAlignment(
          new ProteinSequence(translatedAASequence),
          new ProteinSequence(aaSequence),
          Alignments.PairwiseSequenceAlignerType.GLOBAL,
          new SimpleGapPenalty(6, 5), // GOP and GEP were adjusted to score one less than the worst amino-acid
          // substitution.
          BLOSUM62
      );
      alignedTranslatedAASequence = alignedSequences.getAlignedSequences().get(0).toString().toCharArray();
      alignedAASequence = alignedSequences.getAlignedSequences().get(1).toString().toCharArray();
    } else {
      // Align non-translated amino acid sequence with translated amino acid sequence.
      alignedSequences = Alignments.getPairwiseAlignment(
          new ProteinSequence(aaSequence),
          new ProteinSequence(translatedAASequence),
          Alignments.PairwiseSequenceAlignerType.GLOBAL,
          new SimpleGapPenalty(6, 5), // GOP and GEP were adjusted to score one less than the worst amino-acid
          // substitution.
          BLOSUM62
      );
      alignedTranslatedAASequence = alignedSequences.getAlignedSequences().get(1).toString().toCharArray();
      alignedAASequence = alignedSequences.getAlignedSequences().get(0).toString().toCharArray();
    }
    /*
    (1) Characters truncated at the start are re-appended.
     */
    paddedNucSequenceBuilder.append(truncatedStart);
    paddedAASequenceBuilder.append(String.valueOf(NO_MATCH_CHAR).repeat(truncatedStart.length()));
    paddedAASequenceAnnotations.append("T".repeat(truncatedStart.length()));
    /*
    (2) Iterate over split nucleotide sequence.
    */
    int alignmentIndex = 0;
    String currentCodon;
    for (String s : splitNucleotideSequence) {
      currentCodon = s;
      if (CODON_MAP.get(currentCodon).equals("")) {
        // CASE: Stop codon is stored in split nucleotide sequence.
        // Add coding nucleotides to nucleotide sequence.
        paddedNucSequenceBuilder.append(currentCodon);
        // Add NO_MATCH_CHAR to amino acid sequence.
        paddedAASequenceBuilder.append(String.valueOf(NO_MATCH_CHAR).repeat(3));
        // Add annotation to amino acid sequence annotations.
        paddedAASequenceAnnotations.append("@".repeat(3));
      } else {
        // CASE: AA coding codon is stored in split nucleotide sequence.
        // Check if any gaps are inserted into the translated amino-acid sequence before the next matching alignment
        // position.
        char alignedTranslatedAAChar = alignedTranslatedAASequence[alignmentIndex];
        char alignedAAChar = alignedAASequence[alignmentIndex];
        if (alignedTranslatedAAChar == '-') {
          while (alignedTranslatedAAChar == '-') {
            // Add NO_MATCH_CHAR to nucleotide sequence.
            paddedNucSequenceBuilder.append(String.valueOf(NO_MATCH_CHAR).repeat(3));
            // Add content to amino acid sequence.
            paddedAASequenceBuilder.append(String.valueOf(alignedAAChar).repeat(3));
            // Add annotation to amino acid sequence annotations.
            paddedAASequenceAnnotations.append("X".repeat(3));
            alignmentIndex += 1;
            alignedTranslatedAAChar = alignedTranslatedAASequence[alignmentIndex];
            alignedAAChar = alignedAASequence[alignmentIndex];
          }
        }
        // Add content to nucleotide sequence.
        paddedNucSequenceBuilder.append(currentCodon);
        // Add content to amino acid sequence.
        paddedAASequenceBuilder
            .append(String.valueOf((alignedAAChar == '-') ? NO_MATCH_CHAR : alignedAAChar).repeat(3));
        // Add annotation to amino acid sequence annotations.
        if (alignedAAChar != alignedTranslatedAAChar) {
          paddedAASequenceAnnotations.append("X".repeat(3));
        } else {
          paddedAASequenceAnnotations.append("O".repeat(3));
        }
        alignmentIndex += 1;
      }
    }
    // If alignmentIndex does not match alignment length remaining symbols exist.
    while (alignmentIndex < (alignedAASequence.length - 1)) {
      char alignedTranslatedAAChar = alignedTranslatedAASequence[alignmentIndex];
      char alignedAAChar = alignedAASequence[alignmentIndex];
      // Add content to nucleotide sequence.
      paddedNucSequenceBuilder
          .append(String.valueOf((alignedTranslatedAAChar == '-') ? NO_MATCH_CHAR : alignedTranslatedAAChar).repeat(3));
      // Add content to amino acid sequence.
      paddedAASequenceBuilder.append(String.valueOf((alignedAAChar == '-') ? NO_MATCH_CHAR : alignedAAChar).repeat(3));
      // Add annotation to amino acid sequence annotations.
      if (alignedAAChar != alignedTranslatedAAChar) {
        paddedAASequenceAnnotations.append("X".repeat(3));
      } else {
        paddedAASequenceAnnotations.append("O".repeat(3));
      }
      alignmentIndex += 1;
    }
    /*
    (3) Characters truncated at the end are re-appended.
     */
    paddedNucSequenceBuilder.append(truncatedEnd);
    paddedAASequenceBuilder.append(String.valueOf(NO_MATCH_CHAR).repeat(truncatedEnd.length()));
    paddedAASequenceAnnotations.append("T".repeat(truncatedEnd.length()));
    /*
    (4) Insert gaps in both sequences that are caused by sample insertions.
     */
    int insertedGapsShift = 0;
    for (Integer gapPosition : gapPositions) {
      char[] currentPaddedNucChars = paddedNucSequenceBuilder.toString().toCharArray();
      int noMatchShift = 0;
      for (int i = 0; i < currentPaddedNucChars.length; i++) {
        if (i < gapPosition + noMatchShift) {
          if (currentPaddedNucChars[i] == NO_MATCH_CHAR) {
            noMatchShift += 1;
          }
        }
      }
      paddedNucSequenceBuilder.insert(gapPosition + noMatchShift + insertedGapsShift, "-");
      paddedAASequenceBuilder.insert(gapPosition + noMatchShift + insertedGapsShift, "-");
      paddedAASequenceAnnotations.insert(gapPosition + noMatchShift + insertedGapsShift, "-");
      insertedGapsShift += 1;
    }
    /*
    (5) Insert results into results map.
     */
    String paddedNucSequence = paddedNucSequenceBuilder.toString();
    int paddingNucSequence = (int) paddedNucSequence.codePoints().filter(c -> c == NO_MATCH_CHAR).count();
    String paddedAASequence = paddedAASequenceBuilder.toString();
    int paddingAASequence = (int) paddedAASequence.codePoints().filter(c -> c == NO_MATCH_CHAR).count();
    ArrayList<String> resultsList = new ArrayList<>();
    resultsList.add(paddedAASequence);
    resultsList.add(paddedNucSequence);
    resultsList.add(paddedAASequenceAnnotations.toString());
    resultsMap.put((paddingAASequence + paddingNucSequence), resultsList);
  }

  /**
   * Returns the reverse complement of the passed nucleotide sequence.
   *
   * @param sequence {@link String} representing a nucleotide sequence.
   * @return {@link String} representing the reverse complement of the passed sequence.
   */
  public static String getReverseComplement(String sequence) {
    StringBuilder reverseComplementBuilder = new StringBuilder();
    char[] sequenceCharArray = sequence.toCharArray();
    for (int i = sequenceCharArray.length - 1; i >= 0; i--) {
      char sequenceAtI = sequenceCharArray[i];
      switch (sequenceAtI) {
        case 'A':
          reverseComplementBuilder.append('T');
          continue;
        case 'C':
          reverseComplementBuilder.append('G');
          continue;
        case 'G':
          reverseComplementBuilder.append('C');
          continue;
        case 'T':
          reverseComplementBuilder.append('A');
      }
    }
    return reverseComplementBuilder.toString();
  }

  /**
   * Returns the complement of a nucleotide base, i.e. switches A to T, C to G and vice versa.
   *
   * @param base {@link Character} base to invert.
   * @return {@link Character} the inverted base.
   */
  public static Character invertBase(char base) {
    switch (base) {
      case 'A':
        return 'T';
      case 'C':
        return 'G';
      case 'G':
        return 'C';
      case 'T':
        return 'A';
      default:
        return '?';
    }
  }

  /**
   * Returns the position q of a position p on the reverse complement of a genomic feature.
   * TODO: End should be the full genome length, not only the end position of the feature.
   *
   * @param position {@link Integer} the position that should be converted.
   * @param start    {@link Integer} the start position of the genomic feature on to the full genome.
   * @param end      {@link Integer} the end position of the genomic feature on to the full genome.
   * @return {@link Integer} representing the position on the reverse complement strand.
   */
  public static int getPositionOnReverseComplement(int position, int start, int end) {
    return end - (position - start);
  }

  /**
   * Computes the positions (0-based indexed) at which a gap symbol ("-") is present in a {@link String}.
   *
   * @param sequence The sequence for which gap indices should be computed.
   * @return {@link HashSet<Integer>} of indices at which a gap symbol is present in the passed sequence.
   */
  public static HashSet<Integer> getGapsFromSequence(String sequence) {
    HashSet<Integer> gapIndices = new HashSet<>();
    char[] sequenceCharArray = sequence.toCharArray();
    for (int i = 0; i < sequenceCharArray.length; i++) {
      if (sequenceCharArray[i] == '-') {
        gapIndices.add(i);
      }
    }
    return gapIndices;
  }

  /**
   * Removes all gap symbols ("-") from the passed sequence.
   *
   * @param sequence The sequence from which gaps should be removed.
   * @return {@link String}, the passed sequence without any gap symbols.
   */
  public static String removeGapsFromSequence(String sequence) {
    return sequence.replace("-", "");
  }

}
