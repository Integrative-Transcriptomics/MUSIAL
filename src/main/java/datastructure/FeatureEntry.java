package datastructure;

import exceptions.MusialBioException;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import org.apache.commons.io.FilenameUtils;
import org.javatuples.Triplet;
import utility.Bio;
import utility.IO;
import utility.Logging;

/**
 * Internal representation of a reference sequence location that is subject to analysis. May represent the full
 * genome, a single gene, contigs or plasmids and chromosomes.
 * <p>
 * Pairs of instances of this class and {@link SampleEntry} instances are used internally as so called 'Run
 * entries' to specify a set of analysis tasks, i.e. which sample is analyzed with respect to which specified
 * reference feature.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class FeatureEntry {

  /**
   * The internal name of the entry.
   */
  public final String name;
  /**
   * The reference DNA sequence of the entry.
   */
  public String nucleotideSequence;
  /**
   * The location of the entry on the reference, i.e. for genes the contig or chromosome the feature is located on.
   */
  public final String chromosome;
  /**
   * The 1-based indexed starting position of the feature.
   */
  public final int start;
  /**
   * The 1-based indexed end position of the feature.
   */
  public final int end;
  /**
   * Indicates if the feature is located on the sense strand.
   */
  public final boolean isSense;
  /**
   * TODO
   */
  public AllocatedProteinEntry allocatedProtein;
  /**
   * TODO
   */
  public transient File pdbFile;
  /**
   * TODO
   */
  public final HashMap<String, String> annotations = new HashMap<>();

  /**
   * Constructor of {@link FeatureEntry}.
   *
   * @param entryName     {@link String} representing the internal name of the reference feature to analyze.
   * @param isGene        {@link Boolean} whether the feature represents a single gene or not.
   * @param entryLocation {@link String} the name of the reference location (contig, chromosome, plasmid) the feature
   *                      is located on.
   * @param entryStart    {@link Integer} The 1-based indexed starting position of the feature on the reference.
   * @param entryEnd      {@link Integer} The 1-based indexed end position of the feature on the reference.
   * @throws MusialBioException If the specified locus is ambiguous.
   */
  public FeatureEntry(String entryName, String entryLocation,
                      int entryStart,
                      int entryEnd) throws MusialBioException {
    this.name = entryName;
    this.chromosome = entryLocation;
    if (entryStart >= 0 && entryEnd > 0 && (entryEnd > entryStart)) {
      // CASE: Feature is on sense strand.
      this.isSense = true;
      this.start = entryStart + 1;
      this.end = entryEnd;
    } else if (entryStart < 0 && entryEnd <= 0 && (entryEnd < entryStart)) {
      // CASE: Feature is on anti-sense strand.
      this.isSense = false;
      this.start = -entryStart + 1;
      this.end = -entryEnd;
    } else {
      throw new MusialBioException(
          "Failed to add gene feature due to faulty position data (start, end)\t" + entryStart + ", " + entryEnd);
    }
  }

  public static void imputeSequence(FeatureEntry featureEntry, FastaContainer referenceSequence) {
    featureEntry.nucleotideSequence = referenceSequence.getSequence(featureEntry.start, featureEntry.end);
  }

  public static void imputeProtein(FeatureEntry featureEntry) throws IOException, MusialBioException {
    HashMap<String, String> proteinSequences = IO.getSequencesFromPdbStructure(IO.readStructure(featureEntry.pdbFile));
    String chainId;
    String chainSeq;
    char[] alignedChainSeq;
    ArrayList<String> pdbLines = IO.readLinesFromFile(featureEntry.pdbFile.getAbsolutePath());
    String translFeatureSequence = Bio.translateNucSequence(featureEntry.nucleotideSequence, true, true,
        featureEntry.isSense);
    char[] alignedTranslFeatureSequence;
    StringBuilder paddedChainSeqBuilder = new StringBuilder();
    List<String> diffSegments;
    int prefixPeptideLength = 0;
    String[] splitPDBLine;
    for (Map.Entry<String, String> proteinSequenceEntry : proteinSequences.entrySet()) {
      paddedChainSeqBuilder.setLength(0);
      chainId = proteinSequenceEntry.getKey();
      chainSeq = proteinSequenceEntry.getValue();
      // Align protein and transl. reference sequence to infer correct starting position and truncated prefix peptides.
      Triplet<Integer, String, String> alignedSequences =
          Bio.globalAminoAcidSequenceAlignment(chainSeq, translFeatureSequence,
              Bio.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FREE, Bio.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FREE);
      alignedChainSeq = alignedSequences.getValue1().toCharArray();
      alignedTranslFeatureSequence = alignedSequences.getValue2().toCharArray();
      for (int pos = 0; pos < alignedChainSeq.length; pos++) {
        if (alignedChainSeq[pos] != alignedTranslFeatureSequence[pos]) {
          if (alignedChainSeq[pos] == Bio.DELETION_AA1) {
            paddedChainSeqBuilder.append(Character.toLowerCase(alignedTranslFeatureSequence[pos]));
          } else if (alignedTranslFeatureSequence[pos] == Bio.DELETION_AA1) {
            throw new MusialBioException("Failed to allocate protein " + featureEntry.pdbFile.getName() + " to " +
                "feature " + featureEntry.name + " due to gaps in the aligned translated feature nucleotide sequence.");
          }
        } else {
          paddedChainSeqBuilder.append(alignedChainSeq[pos]);
        }
      }
      diffSegments = Arrays.stream(
          paddedChainSeqBuilder.toString().split("(?=\\p{Upper})")
      ).filter(s -> s.length() > 1).collect(Collectors.toList());
      prefixPeptideLength = Character.isUpperCase(diffSegments.get(0).charAt(0)) ? 0 : diffSegments.get(0).length();
      if (diffSegments.size() > 2) {
        Logging.logWarning("Feature " + featureEntry.name + " disaccords in " + diffSegments.size() + " segments " +
            "with allocated protein " + featureEntry.pdbFile.getName() + ": Structure may be inappropriate.");
      }
      proteinSequences.put(chainId, paddedChainSeqBuilder.toString());
    }
    featureEntry.allocatedProtein = new AllocatedProteinEntry(
        FilenameUtils.removeExtension(featureEntry.pdbFile.getName()),
        String.join(IO.LINE_SEPARATOR, pdbLines),
        proteinSequences
    );
  }

}