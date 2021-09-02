package utility;

import com.google.common.base.Splitter;
import java.util.ArrayList;
import java.util.HashMap;
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
  private static final char PADDING_CHAR = '$';
  /**
   * Character to pad an amino acid sequence if the corresponding nucleotide sequence was truncated due to a frameshift.
   */
  public static final char TRUNCATED_CHAR = '&';
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
   */
  public static ArrayList<String> matchSequences(String aaSequence, String nucSequence)
      throws CompoundNotFoundException {
    // Initialize hash map to store matched sequences.
    HashMap<Integer, ArrayList<String>> paddedSequences = new HashMap<>();
    // Check if nucleotide sequence is fully coding sequence, i.e. divisible by three.
    int residual = nucSequence.length() % 3;
    // Check length difference of sequences.
    int lengthDifference = aaSequence.length() - (nucSequence.length() / 3);
    // Based on the residual assign frame shift and split nucleotide sequence.
    String translatedNucSequence;
    ArrayList<String> splitNucSequence = new ArrayList<>();
    if (residual == 0) {
      // Truncate no characters.
      for (String s : Splitter.fixedLength(3).split(nucSequence)) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, lengthDifference, "", "", paddedSequences);
    } else if (residual == 1) {
      // Truncate character at start.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(1))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, lengthDifference, nucSequence.substring(0, 1)
          , "", paddedSequences);
      // Truncate character at end.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(0, nucSequence.length() - 1))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, lengthDifference, "",
          nucSequence.substring(nucSequence.length() - 1), paddedSequences);
    } else {
      // Truncate two characters at start.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(2))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, lengthDifference, nucSequence.substring(0, 2)
          , "", paddedSequences);
      // Truncate two characters at end.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(0, nucSequence.length() - 2))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, lengthDifference, "",
          nucSequence.substring(nucSequence.length() - 2), paddedSequences);
      // Truncate one character at start and end.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(1, nucSequence.length() - 1))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      alignSequences(aaSequence, translatedNucSequence, splitNucSequence, lengthDifference, nucSequence.substring(0, 1)
          , nucSequence.substring(nucSequence.length() - 1), paddedSequences);
    }
    // Next we choose the sequences as correct with the fewest number of padding characters, i.e. alignment gaps.
    int minScore = Integer.MAX_VALUE;
    ArrayList<String> bestAlignment = new ArrayList<>();
    for (Integer score : paddedSequences.keySet()) {
      if (score < minScore) {
        bestAlignment = paddedSequences.get(score);
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
   * @param aaSequence           {@link String} representing a non-nucleotide-derived amino acid sequence.
   * @param translatedAASequence {@link String} representing a amino acid sequence that originates from translating a
   *                             nucleotide sequence.
   * @param splitNucSequence     {@link String} the original nucleotide sequence split into chunks of length three, i.e.
   *                             codons.
   * @param lengthDifference     {@link Integer} the length difference of the two sequences.
   * @param truncatedStart       {@link String} nucleotides that were truncated from the start of the nucleotide sequence
   *                             due to a frameshift.
   * @param truncatedEnd         {@link String} nucleotides that were truncated from the end of the nucleotide sequence due
   *                             to a frameshift.
   * @param resultsMap           {@link HashMap} mapping the alignment score, here the number of gaps, to the aligned sequences.
   * @throws CompoundNotFoundException If a compound, i.e. amino acid, could not be translated.
   */
  private static void alignSequences(String aaSequence, String translatedAASequence, ArrayList<String> splitNucSequence,
                                     int lengthDifference, String truncatedStart, String truncatedEnd, HashMap<Integer,
      ArrayList<String>> resultsMap)
      throws CompoundNotFoundException {
    StringBuilder paddedNucSequenceBuilder = new StringBuilder();
    StringBuilder paddedAASequenceBuilder = new StringBuilder();
    String alignedAASequence;
    String alignedTranslatedAASequence;
    if (lengthDifference > 0) {
      // CASE: Amino acid sequence is longer.
      SequencePair<ProteinSequence, AminoAcidCompound> alignedSequences = Alignments.getPairwiseAlignment(
          new ProteinSequence(translatedAASequence),
          new ProteinSequence(aaSequence),
          Alignments.PairwiseSequenceAlignerType.LOCAL,
          new SimpleGapPenalty(6, 5), // GOP and GEP were adjusted to score one less than the worst amino-acid
          // substitution.
          BLOSUM62
      );
      alignedTranslatedAASequence = alignedSequences.getQuery().getSequenceAsString();
      alignedAASequence = alignedSequences.getTarget().getSequenceAsString();
    } else {
      // CASE: Translated amino acid sequence is longer or sequences match in their length.
      SequencePair<ProteinSequence, AminoAcidCompound> alignedSequences = Alignments.getPairwiseAlignment(
          new ProteinSequence(aaSequence),
          new ProteinSequence(translatedAASequence),
          Alignments.PairwiseSequenceAlignerType.LOCAL,
          new SimpleGapPenalty(6, 5),
          BLOSUM62
      );
      alignedTranslatedAASequence = alignedSequences.getTarget().getSequenceAsString();
      alignedAASequence = alignedSequences.getQuery().getSequenceAsString();
    }
    // Replace gap characters of aligned amino acid sequence with padding characters.
    paddedAASequenceBuilder.append(String.valueOf(TRUNCATED_CHAR).repeat(truncatedStart.length()));
    paddedAASequenceBuilder.append(alignedAASequence.replace('-', PADDING_CHAR));
    paddedAASequenceBuilder.append(String.valueOf(TRUNCATED_CHAR).repeat(truncatedEnd.length()));
    String paddedAASequence = paddedAASequenceBuilder.toString();
    int paddingAASequence = (int) paddedAASequence.codePoints().filter(c -> c == PADDING_CHAR).count();
    // Translate translated aligned nucleotide sequence back into padded nucleotide sequence.
    paddedNucSequenceBuilder.append(truncatedStart);
    char[] alignedTranslatedAASequenceChars = alignedTranslatedAASequence.toCharArray();
    int gapShift = 0;
    for (int i = 0; i < alignedTranslatedAASequenceChars.length; i++) {
      char c = alignedTranslatedAASequenceChars[i];
      if (c == '-') {
        paddedNucSequenceBuilder.append(String.valueOf(PADDING_CHAR).repeat(3));
        gapShift += 1;
      } else {
        paddedNucSequenceBuilder.append(splitNucSequence.get(i - gapShift));
      }
    }
    paddedNucSequenceBuilder.append(truncatedEnd);
    String paddedNucSequence = paddedNucSequenceBuilder.toString();
    int paddingNucSequence = (int) paddedNucSequence.codePoints().filter(c -> c == PADDING_CHAR).count();
    // Insert results into results map.
    ArrayList<String> resultsList = new ArrayList<>();
    resultsList.add(paddedAASequence);
    resultsList.add(paddedNucSequence);
    resultsMap.put((paddingAASequence + paddingNucSequence), resultsList);
  }

}
