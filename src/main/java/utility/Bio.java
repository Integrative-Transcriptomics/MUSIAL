package utility;

import com.google.common.base.Splitter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;

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
  public static final char NO_MATCH_CHAR = '_';
  /**
   * Character to indicate a truncation of a sequence.
   */
  public static final char TRUNCATED_CHAR = '#';
  /**
   * Character to indicate stop codon.
   */
  public static final char STOP_CODON_CHAR = '@';
  /**
   * Character used as a separator in padded sequences.
   */
  public static final char SEPARATOR_CHAR = '&';
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
    put("B", "ASX");
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
  }};
  /**
   * TODO: Use other or own substitution matrix, especially to consider internal stop codons, if present.
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
    // Remove all gap symbols in the nucleotide sequence (this would rise errors in the sequence
    // alignment process).
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
      padSequences(aaSequence, translatedNucSequence, splitNucSequence, "", "", alignedSequences);
    } else if (residual == 1) {
      // Truncate character at start.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(1))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      padSequences(aaSequence, translatedNucSequence, splitNucSequence, nucSequence.substring(0, 1), "",
          alignedSequences);
      // Truncate character at end.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(0, nucSequence.length() - 1))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      padSequences(aaSequence, translatedNucSequence, splitNucSequence, "",
          nucSequence.substring(nucSequence.length() - 1), alignedSequences);
    } else {
      // Truncate two characters at start.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(2))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      padSequences(aaSequence, translatedNucSequence, splitNucSequence, nucSequence.substring(0, 2), "",
          alignedSequences);
      // Truncate two characters at end.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(0, nucSequence.length() - 2))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      padSequences(aaSequence, translatedNucSequence, splitNucSequence, "",
          nucSequence.substring(nucSequence.length() - 2), alignedSequences);
      // Truncate one character at start and end.
      for (String s : Splitter.fixedLength(3).split(nucSequence.substring(1, nucSequence.length() - 1))) {
        splitNucSequence.add(s);
      }
      translatedNucSequence = translateNucSequence(splitNucSequence);
      padSequences(aaSequence, translatedNucSequence, splitNucSequence, nucSequence.substring(0, 1),
          nucSequence.substring(nucSequence.length() - 1), alignedSequences);
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
   * <p>
   * Insertions are currently handled separate and not included into the padded sequences directly.
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
   * @throws CompoundNotFoundException If a compound, i.e. amino acid, could not be translated.
   */
  private static void padSequences(String aaSequence, String translatedAASequence,
                                   ArrayList<String> splitNucleotideSequence, String truncatedStart,
                                   String truncatedEnd, HashMap<Integer, ArrayList<String>> resultsMap)
      throws CompoundNotFoundException {
    StringBuilder paddedNucSequenceBuilder = new StringBuilder();
    StringBuilder paddedAASequenceBuilder = new StringBuilder();
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
          // substitution in the BLOSUM62 matrix.
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
          // substitution in the BLOSUM62 matrix.
          BLOSUM62
      );
      alignedTranslatedAASequence = alignedSequences.getAlignedSequences().get(1).toString().toCharArray();
      alignedAASequence = alignedSequences.getAlignedSequences().get(0).toString().toCharArray();
    }
    /*
    (1) Characters truncated at the start are re-appended.
     */
    if ( truncatedStart.length() != 0 ) {
      paddedNucSequenceBuilder.append(truncatedStart).append(SEPARATOR_CHAR);
      paddedAASequenceBuilder.append(TRUNCATED_CHAR).append(SEPARATOR_CHAR);
    }
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
        paddedNucSequenceBuilder.append(currentCodon).append(SEPARATOR_CHAR);
        // Add NO_MATCH_CHAR to amino acid sequence.
        paddedAASequenceBuilder.append(STOP_CODON_CHAR).append(SEPARATOR_CHAR);
      } else {
        // CASE: AA coding codon is stored in split nucleotide sequence.
        // Check if any gaps are inserted into the translated amino-acid sequence before the next matching alignment
        // position.
        char alignedTranslatedAAChar = alignedTranslatedAASequence[alignmentIndex];
        char alignedAAChar = alignedAASequence[alignmentIndex];
        if (alignedTranslatedAAChar == '-') {
          while (alignedTranslatedAAChar == '-') {
            // Add NO_MATCH_CHAR to nucleotide sequence.
            paddedNucSequenceBuilder.append(String.valueOf(NO_MATCH_CHAR).repeat(3)).append(SEPARATOR_CHAR);
            // Add content to amino acid sequence.
            paddedAASequenceBuilder.append(AA1TO3.get(String.valueOf(alignedAAChar))).append(SEPARATOR_CHAR);
            alignmentIndex += 1;
            alignedTranslatedAAChar = alignedTranslatedAASequence[alignmentIndex];
            alignedAAChar = alignedAASequence[alignmentIndex];
          }
        }
        // Add content to nucleotide sequence.
        paddedNucSequenceBuilder.append(currentCodon).append(SEPARATOR_CHAR);
        // Add content to amino acid sequence.
        paddedAASequenceBuilder
            .append((alignedAAChar == '-') ? String.valueOf(NO_MATCH_CHAR).repeat(3) :
                AA1TO3.get(String.valueOf(alignedAAChar)))
            .append(SEPARATOR_CHAR);
        alignmentIndex += 1;
      }
    }
    // If alignmentIndex does not match alignment length we have a trailing "gap" in the translated nucleotide sequence.
    while (alignmentIndex < (alignedAASequence.length - 1)) {
      char alignedTranslatedAAChar = alignedTranslatedAASequence[alignmentIndex];
      char alignedAAChar = alignedAASequence[alignmentIndex];
      // Add content to nucleotide sequence.
      paddedNucSequenceBuilder
          .append((alignedTranslatedAAChar == '-') ? String.valueOf(NO_MATCH_CHAR).repeat(3) :
              "?".repeat(3)).append(SEPARATOR_CHAR);
      // Add content to amino acid sequence.
      paddedAASequenceBuilder.append((alignedAAChar == '-') ? String.valueOf(NO_MATCH_CHAR).repeat(3) :
          AA1TO3.get(String.valueOf(alignedAAChar))).append(SEPARATOR_CHAR);
      alignmentIndex += 1;
    }
    /*
    (3) Characters truncated at the end are re-appended.
     */
    if ( truncatedEnd.length() != 0 ) {
      paddedNucSequenceBuilder.append(truncatedEnd).append(SEPARATOR_CHAR);
      paddedAASequenceBuilder.append(TRUNCATED_CHAR).append(SEPARATOR_CHAR);
    }
    /*
    DEPRECATED:
    (4) Insert gaps in both sequences that are caused by sample insertions.
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
      insertedGapsShift += 1;
    }
     */
    /*
    (5) Insert results into results map.
     */
    paddedNucSequenceBuilder.setLength( paddedNucSequenceBuilder.length() - 1 );
    String paddedNucSequence = paddedNucSequenceBuilder.toString();
    int paddingNucSequence = (int) paddedNucSequence.codePoints().filter(c -> c == NO_MATCH_CHAR).count();
    paddedAASequenceBuilder.setLength( paddedAASequenceBuilder.length() - 1 );
    String paddedAASequence = paddedAASequenceBuilder.toString();
    int paddingAASequence = (int) paddedAASequence.codePoints().filter(c -> c == NO_MATCH_CHAR).count();
    ArrayList<String> resultsList = new ArrayList<>();
    resultsList.add(paddedAASequence);
    resultsList.add(paddedNucSequence);
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
      reverseComplementBuilder.append(invertBase(sequenceAtI));
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
    return switch (base) {
      case 'A' -> 'T';
      case 'a' -> 't';
      case 'C' -> 'G';
      case 'c' -> 'g';
      case 'G' -> 'C';
      case 'g' -> 'c';
      case 'T' -> 'A';
      case 't' -> 'a';
      case '.' -> '.';
      case ':' -> ':';
      case 'N' -> 'N';
      case '-' -> '-';
      default -> '?';
    };
  }

  /**
   * Returns the position q of a position p on the reverse complement of a genomic feature.
   *
   * @param position       {@link Integer} the position that should be converted.
   * @param sequenceLength {@link Integer} the length of the parent sequence (of the genomic
   *                       feature) or full contig/chromosome/plasmid.
   * @return {@link Integer} representing the position on the reverse complement strand.
   */
  public static int getPositionOnReverseComplement(int position, int sequenceLength) {
    return sequenceLength + (1 - position);
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
