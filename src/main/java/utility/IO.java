package utility;

import datastructure.FastaEntry;
import datastructure.FeatureAnalysisEntry;
import datastructure.VariablePosition;
import datastructure.VariablePositionsTable;
import exceptions.MusialIOException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.TreeSet;
import main.Musial;
import org.apache.commons.io.FilenameUtils;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import tools.ArgumentsParser;

/**
 * This class comprises static methods used for reading and writing files.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class IO {
  /**
   * OS dependent line separator
   */
  public static final String LINE_SEPARATOR = System.getProperty("line.separator");

  /**
   * Reads a file and returns its content line-wise as list.
   *
   * @param filePath {@link String} representing a file path.
   * @return {@link ArrayList} of {@link String} objects, each representing one line of the file accessed via the
   * passed file path.
   * @throws FileNotFoundException If the specified file path does not lead to any file.
   */
  public static List<String> getLinesFromFile(String filePath) throws FileNotFoundException {
    ArrayList<String> lines = new ArrayList<>();
    Scanner scanner = new Scanner(new File(filePath));
    while (scanner.hasNextLine()) {
      String nextLine = scanner.nextLine();
      if (!nextLine.trim().isEmpty()) {
        lines.add(nextLine.trim());
      }
    }
    scanner.close();
    return lines;
  }

  /**
   * Reads a .pdb file an returns the parsed file as {@link Structure}.
   *
   * @param file {@link File} object pointing to the .pdb file to parse.
   * @return A {@link Structure} instance representation of the specified .pdb file.
   * @throws MusialIOException If the .pdb file can not be parsed.
   */
  public static Structure getStructureFromPdb(File file) throws MusialIOException {
    // Initialize pdbReader.
    PDBFileReader pdbReader = new PDBFileReader();
    FileParsingParameters params = new FileParsingParameters();
    params.setParseSecStruc(true);
    params.setAlignSeqRes(true);
    pdbReader.setFileParsingParameters(params);
    // Read the `.pdb` file.
    Structure pdbStructure;
    try {
      pdbStructure = pdbReader.getStructure(file);
    } catch (IOException e) {
      throw new MusialIOException("Failed to read file:\t" + file.getAbsolutePath());
    }
    return pdbStructure;
  }

  /**
   * Reads a .gff file from the path specified by a {@link File} object and returns a
   * {@link FeatureList} object.
   *
   * @param file {@link File} specifying a .gff file.
   * @return {@link FeatureList} comprising the entries from the read .gff file.
   * @throws IOException If the {@link GFF3Reader} throws any {@link IOException}.
   */
  public static FeatureList readGFF(File file) throws IOException {
    System.setOut(Musial.EMPTY_STREAM);
    System.setErr(Musial.EMPTY_STREAM);
    FeatureList features = GFF3Reader.read(file.getCanonicalPath());
    System.setOut(Musial.ORIGINAL_OUT_STREAM);
    System.setErr(Musial.ORIGINAL_ERR_STREAM);
    return features;
  }

  /**
   * Reads a .fasta file from a {@link File} object into a {@link HashSet} of {@link FastaEntry} instances. Each
   * element contains the header and sequence of the respective entry from the .fasta file.
   *
   * @param file {@link File} specifying a .fasta file.
   * @return {@link HashSet} of {@link FastaEntry} instances.
   * @throws IOException If the respective file can not be found or accessed.
   */
  public static HashSet<FastaEntry> readFasta(File file) throws IOException {
    HashSet<FastaEntry> fastaEntries = new HashSet<>();
    @SuppressWarnings("resource")
    BufferedReader br = new BufferedReader(new FileReader(file));
    String currLine;
    String currHeader = "";
    StringBuilder currSequence = new StringBuilder();
    while ((currLine = br.readLine()) != null) {
      if (currLine.startsWith(">")) {
        if (currHeader.length() > 0) {
          fastaEntries.add(new FastaEntry(currHeader, currSequence.toString()));
          currSequence = new StringBuilder();
        }
        currHeader = currLine.substring(1);
      } else if (!currLine.startsWith(";")) {
        currSequence.append(currLine.trim());
      }
    }
    fastaEntries.add(new FastaEntry(currHeader, currSequence.toString()));
    return fastaEntries;
  }

  /**
   * Accepts a {@link File} object representing a directory that is subject to deletion.
   *
   * @param file {@link File} specifying the directory to delete.
   * @throws IOException If any file deletion procedure failed, for example the specified directory does not exist.
   */
  public static void deleteDirectory(File file) throws IOException {
    if (file.isDirectory()) {
      //noinspection ResultOfMethodCallIgnored
      Files.walk(file.toPath())
          .sorted(Comparator.reverseOrder())
          .map(Path::toFile)
          .forEach(File::delete);
    }
  }

  /**
   * Tries to generate the directory specified by the passed {@link File} object.
   *
   * @param file {@link File} object representing a directory.
   * @throws MusialIOException If the directory could not be generated.
   */
  public static void generateDirectory(File file) throws MusialIOException {
    if (!file.mkdirs()) {
      throw new MusialIOException("Failed to generate output directory:\t" + file.getAbsolutePath());
    }
  }

  /**
   * Tries to generate the file specified by the passed {@link File} object.
   *
   * @param file {@link File} object representing a file.
   * @throws MusialIOException If the file could not be generated.
   */
  public static void generateFile(File file) throws MusialIOException {
    try {
      if (!file.createNewFile()) {
        throw new MusialIOException("Failed to generate output directory:\t" + file.getAbsolutePath());
      }
    } catch (IOException e) {
      throw new MusialIOException("Failed to generate output directory:\t" + file.getAbsolutePath());
    }

  }

  /**
   * Copies the file pointed to with file to the file pointed to with target.
   *
   * @param file   {@link File} object, the source file.
   * @param target {@link File} object, the target file.
   * @throws MusialIOException If the copy procedure fails.
   */
  public static void copyFile(File file, File target) throws MusialIOException {
    try {
      Files.copy(file.getAbsoluteFile().toPath(), target.getAbsoluteFile().toPath());
    } catch (IOException e) {
      throw new MusialIOException(
          "Failed to copy file " + file.getAbsolutePath() + " to target " + target.getAbsolutePath());
    }
  }

  /**
   * Writes output `.tsv` file containing nucleotide variants in a tab delimited format.
   * <p>
   * The expected output format is:
   * Position \t Reference (optional) \t [SAMPLE_1_NAME] \t ...
   * x (Integer) \t [Nucleotide at position x in reference] (optional) \t [Nucleotide at position x in sample 1] \t ...
   *
   * @param outFile                {@link File} instance pointing to the output file.
   * @param referenceAnalysisId    {@link String} internal reference feature identifier of which the output is written.
   * @param variablePositionsTable {@link VariablePositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeSnvTable(File outFile, String referenceAnalysisId,
                                   VariablePositionsTable variablePositionsTable)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter snvTableWriter = new FileWriter(outFile.getAbsolutePath());
      snvTableWriter.write(
          variablePositionsTable.getSampleNamesHeaderTabDelimited(referenceAnalysisId) + LINE_SEPARATOR
      );
      Iterator<String> variantPositions = variablePositionsTable.getVariantPositions(referenceAnalysisId);
      // Initialize the shift caused by any insertions with 0.
      int shift = 0;
      while (variantPositions.hasNext()) {
        String variantPosition = variantPositions.next();
        // If the current position contains an "I" the shift is increased by one. As the positions are sorted with
        // respect to the number of "I" characters, e.g. 1,2,2I,2II,2III,3 for an insertion of length 3 at position
        // 2, and positions are iterated over in ascending order, the number of "I"s is not relevant.
        if (variantPosition.contains("I")) {
          shift += 1;
        }
        snvTableWriter.write(
            variablePositionsTable
                .getPositionContentTabDelimited(referenceAnalysisId, variantPosition, shift) +
                LINE_SEPARATOR
        );
      }
      snvTableWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Writes output `.tsv` file containing annotations in a tab delimited format.
   * <p>
   * The expected output format is:
   * Position \t [SAMPLE_1_NAME] \t ...
   * x (Integer) \t [Quality,Coverage,...] \t ...
   *
   * @param outFile                {@link File} instance pointing to the output file.
   * @param referenceAnalysisId    {@link String} internal reference feature identifier of which the output is written.
   * @param variablePositionsTable {@link VariablePositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeSnvAnnotations(File outFile, String referenceAnalysisId,
                                         VariablePositionsTable variablePositionsTable) throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter snvAnnotationsWriter = new FileWriter(outFile.getAbsolutePath());
      snvAnnotationsWriter.write(
          variablePositionsTable.getSampleNamesHeaderTabDelimited(referenceAnalysisId) + LINE_SEPARATOR
      );
      Iterator<String> variantPositions = variablePositionsTable.getVariantPositions(referenceAnalysisId);
      int shift = 0;
      while (variantPositions.hasNext()) {
        String variantPosition = variantPositions.next();
        if (variantPosition.contains("I")) {
          shift += 1;
        }
        snvAnnotationsWriter.write(
            variablePositionsTable.getPositionAnnotationTabDelimited(referenceAnalysisId, variantPosition, shift) +
                LINE_SEPARATOR
        );
      }
      snvAnnotationsWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Writes output `.fasta` file containing the sequence of variants of all samples and the reference (optional).
   *
   * @param outFile                {@link File} instance pointing to the output file.
   * @param referenceAnalysisId    {@link String} internal reference feature identifier of which the output is written.
   * @param variablePositionsTable {@link VariablePositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeVariantAlignment(File outFile, String referenceAnalysisId,
                                           VariablePositionsTable variablePositionsTable)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter variantAlignmentWriter = new FileWriter(outFile.getAbsolutePath());
      variantAlignmentWriter.write(variablePositionsTable.getSampleVariantsSequence(referenceAnalysisId,
          "Reference", true));
      for (String sampleName : variablePositionsTable.getSampleNamesOf(referenceAnalysisId)) {
        variantAlignmentWriter.write(variablePositionsTable.getSampleVariantsSequence(referenceAnalysisId, sampleName
            , true));
      }
      variantAlignmentWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Writes output `.fasta` file containing the full length sequences of all samples and the reference (optional).
   *
   * @param outFile                {@link File} instance pointing to the output file.
   * @param referenceAnalysisId    {@link String} internal reference feature identifier of which the output is written.
   * @param variablePositionsTable {@link VariablePositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeFullAlignment(File outFile, String referenceAnalysisId,
                                        VariablePositionsTable variablePositionsTable)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter sequenceAlignmentWriter = new FileWriter(outFile.getAbsolutePath());
      sequenceAlignmentWriter.write(variablePositionsTable.getSampleFullSequence(referenceAnalysisId,
          "Reference", true));
      for (String sampleName : variablePositionsTable.getSampleNamesOf(referenceAnalysisId)) {
        sequenceAlignmentWriter.write(variablePositionsTable.getSampleFullSequence(referenceAnalysisId, sampleName,
            true));
      }
      sequenceAlignmentWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Writes a `.tsv` file containing per sample statistics.
   *
   * @param outFile                {@link File} instance pointing to the output file.
   * @param referenceAnalysisId    {@link String} internal reference feature identifier of which the output is written.
   * @param variablePositionsTable {@link VariablePositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeSampleStatistics(File outFile, String referenceAnalysisId,
                                           VariablePositionsTable variablePositionsTable) throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter sampleStatisticsWriter = new FileWriter(outFile.getAbsolutePath());
      sampleStatisticsWriter.write(
          variablePositionsTable.getStatisticsHeaderTabDelimited("Sample") + LINE_SEPARATOR
      );
      List<String> sampleNames = variablePositionsTable.getSampleNamesOf(referenceAnalysisId);
      for (String sampleName : sampleNames) {
        sampleStatisticsWriter.write(
            variablePositionsTable.getSampleStatisticsTabDelimited(referenceAnalysisId, sampleName) +
                LINE_SEPARATOR
        );
      }
      sampleStatisticsWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Writes a `.tsv` file containing per position statistics.
   * <p>
   * One row contains information about the counts of each possible variant content (see {@link datastructure.VariablePosition})
   *
   * @param outFile                {@link File} instance pointing to the output file.
   * @param referenceAnalysisId    {@link String} internal reference feature identifier of which the output is written.
   * @param variablePositionsTable {@link VariablePositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writePositionStatistics(File outFile, String referenceAnalysisId,
                                             VariablePositionsTable variablePositionsTable) throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter positionStatisticsWriter = new FileWriter(outFile.getAbsolutePath());
      positionStatisticsWriter.write(
          variablePositionsTable.getStatisticsHeaderTabDelimited("Position") + LINE_SEPARATOR
      );
      Iterator<String> variantPositions = variablePositionsTable.getVariantPositions(referenceAnalysisId);
      int shift = 0;
      while (variantPositions.hasNext()) {
        String variantPosition = variantPositions.next();
        if (variantPosition.contains("I")) {
          shift += 1;
        }
        positionStatisticsWriter.write(
            variablePositionsTable.getPositionStatisticsTabDelimited(referenceAnalysisId, variantPosition, shift) +
                LINE_SEPARATOR
        );
      }
      positionStatisticsWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Writes a nucleotide to amino-acid sequence alignment for a specified feature (referenceAnalysisId) and protein
   * structure (.pdb file accessed via referenceAnalysisId).
   * <p>
   * One alignment and one output file is generated per chain of the protein structure.
   *
   * @param outDirectory           {@link File} object specifying the output directory.
   * @param featureAnalysisEntry   {@link FeatureAnalysisEntry} specifying the reference feature to analyze.
   * @param variablePositionsTable {@link VariablePositionsTable} storing the corresponding data.
   * @param arguments              {@link ArgumentsParser} arguments parsed from the command line.
   * @throws MusialIOException If any output file does not exist or could not be accessed.
   */
  public static void writeProteinAlignment(File outDirectory, FeatureAnalysisEntry featureAnalysisEntry,
                                           VariablePositionsTable variablePositionsTable, ArgumentsParser arguments)
      throws MusialIOException {
    if (!outDirectory.exists()) {
      throw new MusialIOException("The specified output directory does not exist:\t" + outDirectory.getAbsolutePath());
    }
    File proteinAlignmentOutFile = null;
    try {
      String referenceAnalysisId = featureAnalysisEntry.identifier;
      String referenceSequence = variablePositionsTable.getSampleFullSequence(referenceAnalysisId,
          "Reference", false);
      if (!featureAnalysisEntry.isSense) {
        referenceSequence = Bio.getReverseComplement(referenceSequence);
      }
      File pdbFile = arguments.getPdInputFiles().get(referenceAnalysisId);
      String pdbName = FilenameUtils.removeExtension(pdbFile.getName());
      Structure pdbStructure = getStructureFromPdb(pdbFile);
      HashMap<String, String> proteinSequences = Bio.getSequencesFromPdb(pdbStructure);
      FileWriter proteinAlignmentWriter;
      StringBuilder lineStringBuilder = new StringBuilder();
      for (Chain chain : pdbStructure.getChains()) {
        String chainId = chain.getName();
        String chainAsymId = chain.getId();
        proteinAlignmentOutFile = new File(outDirectory.getAbsolutePath() + "/proteinAlignment_" + chainId +
            ".tsv");
        IO.generateFile(proteinAlignmentOutFile);
        proteinAlignmentWriter = new FileWriter(proteinAlignmentOutFile);
        // (1) Write file header.
        proteinAlignmentWriter.write("#fileformat\tASFv1.0" + LINE_SEPARATOR);
        proteinAlignmentWriter
            .write("#date\t" + DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss").format(LocalDateTime.now()) +
                LINE_SEPARATOR);
        proteinAlignmentWriter.write("#source\t" + Musial.CLASS_NAME + Musial.VERSION + LINE_SEPARATOR);
        proteinAlignmentWriter
            .write("#input\tPROTEIN=" + pdbName + "\tNUCLEOTIDES=" + referenceAnalysisId + LINE_SEPARATOR);
        proteinAlignmentWriter.write("Type\tGenomePositions\tGenomeContent\tEncodedAminoAcid\tProteinPosition" +
            "\tProteinContent\tMissenseMutations" + LINE_SEPARATOR);
        // Match amino-acid and nucleotide sequence.
        String aminoAcidSequence = proteinSequences.get(chainId);
        ArrayList<String> matchedSequences = Bio.matchSequences(aminoAcidSequence, referenceSequence);
        String alignedAminoAcidSequence = matchedSequences.get(0);
        Iterator<String> alignedAminoAcidSequenceIterator =
            Arrays.stream(alignedAminoAcidSequence.split(String.valueOf(Bio.SEPARATOR_CHAR))).iterator();
        String alignedReferenceSequence = matchedSequences.get(1);
        Iterator<String> alignedReferenceSequenceIterator =
            Arrays.stream(alignedReferenceSequence.split(String.valueOf(Bio.SEPARATOR_CHAR))).iterator();
        Iterator<String> insertionPositionsIterator =
            variablePositionsTable.getInsertionPositions(referenceAnalysisId);
        int aminoAcidGapCount = 0;
        int aminoAcidPosition = Bio.getFirstResidueNumberFromProteinStructure(pdbStructure, chainAsymId);
        // TODO: What happens if features positions start with an insertion from any sample? Is this even possible
        //  from .vcf files?
        TreeSet<String> featurePositions = variablePositionsTable.getAllPositions(referenceAnalysisId);
        int referenceGapCount = 0;
        int referencePosition;
        if (featureAnalysisEntry.isSense) {
          referencePosition = Integer.parseInt(featurePositions.first());
        } else {
          referencePosition = -Integer.parseInt(featurePositions.last());
        }
        // (2) Write alignment information to file.
        String alignedAminoAcidSequenceEntry;
        String alignedReferenceSequenceEntry;
        while (alignedAminoAcidSequenceIterator.hasNext() || alignedReferenceSequenceIterator.hasNext()) {
          alignedAminoAcidSequenceEntry = alignedAminoAcidSequenceIterator.next();
          alignedReferenceSequenceEntry = alignedReferenceSequenceIterator.next();
          // (2i) Write truncated sequence information.
          if (alignedAminoAcidSequenceEntry.equals(String.valueOf(Bio.TRUNCATED_CHAR))) {
            lineStringBuilder.append("TRC\t");
            int truncationLength = alignedReferenceSequenceEntry.length();
            for (int i = 0; i < truncationLength; i++) {
              if (i == truncationLength - 1) {
                lineStringBuilder.append(referencePosition).append("\t");
              } else {
                lineStringBuilder.append(referencePosition).append(",");
              }
              referencePosition += 1;
            }
            lineStringBuilder.append(alignedReferenceSequenceEntry).append("\t\t\t\t").append(LINE_SEPARATOR);
            proteinAlignmentWriter.write(lineStringBuilder.toString());
            lineStringBuilder = new StringBuilder();
          } else {
            // Write alignment information.
            lineStringBuilder.append("ALN\t");
            // Whether to check for mis-sense mutations.
            boolean referenceIsAligned = false;
            // Write information from nucleotide/reference sequence.
            if (alignedReferenceSequenceEntry.equals(String.valueOf(Bio.NO_MATCH_CHAR).repeat(3))) {
              // CASE: An alignment gap is present in the nucleotide sequence.
              lineStringBuilder.append("_\t_\t_\t");
            } else {
              // CASE: Content from nucleotide sequence is aligned.
              // Add information about genome position.
              int referenceEntryLength = alignedReferenceSequenceEntry.length();
              for (int i = 0; i < referenceEntryLength; i++) {
                if (i == referenceEntryLength - 1) {
                  lineStringBuilder.append(referencePosition).append("\t");
                } else {
                  lineStringBuilder.append(referencePosition).append(",");
                }
                referencePosition += 1;
              }
              // Add information about genome content.
              lineStringBuilder.append(alignedReferenceSequenceEntry).append("\t");
              // Add information about amino-acid encoded by genome content.
              String codon = Bio.CODON_MAP.get(alignedReferenceSequenceEntry);
              if (codon.equals("")) {
                // Stop codon.
                lineStringBuilder.append("STOP").append("\t");
              } else {
                // Encoded amino acid.
                lineStringBuilder.append(Bio.AA1TO3.get(codon)).append("\t");
              }
              referenceIsAligned = true;
            }
            // Write information from amino-acid sequence.
            if (alignedAminoAcidSequenceEntry.equals(String.valueOf(Bio.NO_MATCH_CHAR).repeat(3)) ||
                alignedAminoAcidSequenceEntry.equals(String.valueOf(Bio.STOP_CODON_CHAR))) {
              // CASE: An alignment gap is present in the amino-acid sequence.
              lineStringBuilder.append("_\t_\t");
            } else {
              // CASE: Content from amino-acid sequence is aligned.
              lineStringBuilder.append(aminoAcidPosition).append("\t").append(alignedAminoAcidSequenceEntry).append(
                  "\t");
              aminoAcidPosition += 1;
            }
            // Write information about possible mis-sense mutations from samples.
            boolean hasMissense = false;
            if (referenceIsAligned) {
              // CASE: Search for mis-sense mutations in samples.
              for (String sampleName : variablePositionsTable.getSampleNamesOf(referenceAnalysisId)) {
                char[] codon = alignedReferenceSequenceEntry.toCharArray();
                char[] alternateCodon = Arrays.copyOf(codon, 3);
                String codonString = new String(codon);
                String alternateCodonString;
                String referenceAminoAcid = Bio.CODON_MAP.get(codonString);
                String alternativeAminoAcid;
                VariablePosition variablePosition;
                // Check each reference position of the current codon if any alternate calls for the sample are present.
                // Given that the position is i, then the positions i-3, i-2 and i-1 have to be looked at as the
                // increment of the position in the steps before is already at the next codons position.
                for (int i = 3; i > 0; i--) {
                  if (featureAnalysisEntry.isSense) {
                    variablePosition = variablePositionsTable
                        .getVariablePosition(referenceAnalysisId, sampleName, String.valueOf(referencePosition - i));
                  } else {
                    variablePosition = variablePositionsTable
                        .getVariablePosition(referenceAnalysisId, sampleName, String.valueOf(-referencePosition + i));
                  }
                  if (variablePosition != null &&
                      variablePosition.content != VariablePosition.NO_CALL &&
                      variablePosition.content != VariablePosition.REFERENCE &&
                      variablePosition.content != VariablePosition.REFERENCE_DISCARDED) {
                    if (featureAnalysisEntry.isSense) {
                      alternateCodon[3 - i] = variablePosition.content;
                    } else {
                      alternateCodon[3 - i] = Bio.invertBase(variablePosition.content);
                    }
                  } else {
                    alternateCodon[3 - i] = codon[3 - i];
                  }
                }
                alternateCodonString = new String(alternateCodon);
                if (codonString.contains("-")) {
                  // CASE: Deletion in sample leads to frame shift.
                  hasMissense = true;
                  lineStringBuilder.append(sampleName).append("|").append(alternateCodonString).append("|DEL")
                      .append(",");
                } else {
                  // CASE: Possible alternative amino acid due to mutation.
                  alternativeAminoAcid = Bio.CODON_MAP.get(alternateCodonString.toUpperCase());
                  if (!referenceAminoAcid.equals(alternativeAminoAcid)) {
                    hasMissense = true;
                    if (alternativeAminoAcid.equals("")) {
                      // Alternate codon is a stop codon-
                      lineStringBuilder.append(sampleName).append("|").append(alternateCodonString).append("|")
                          .append("STOP").append(",");
                    } else {
                      // Alternate codon encodes a non-reference amino-acid.
                      lineStringBuilder.append(sampleName).append("|").append(alternateCodonString).append("|")
                          .append(Bio.AA1TO3.get(alternativeAminoAcid)).append(",");
                    }
                  }
                }
              }
              // Remove last appended , character.
              if (hasMissense) {
                lineStringBuilder.setLength(lineStringBuilder.length() - 1);
              }
              lineStringBuilder.append(LINE_SEPARATOR);
            } else {
              // CASE: No reference sequence position is aligned.
              lineStringBuilder.append("_").append(LINE_SEPARATOR);
            }
          }
          // Write line and empty line string builder.
          proteinAlignmentWriter.write(lineStringBuilder.toString());
          lineStringBuilder = new StringBuilder();
        }
        // (4) Write insertion positions information.
        String insertionPosition;
        VariablePosition variablePosition;
        while (insertionPositionsIterator.hasNext()) {
          insertionPosition = insertionPositionsIterator.next();
          lineStringBuilder.append("INS\t").append(insertionPosition).append("\t\t\t\t");
          // Check samples with insertions at the specified positions.
          for (String sampleName : variablePositionsTable.getSampleNamesOf(referenceAnalysisId)) {
            variablePosition = variablePositionsTable.getVariablePosition(referenceAnalysisId, sampleName,
                insertionPosition);
            if (variablePosition != null) {
              lineStringBuilder.append(sampleName).append("|").append(variablePosition.content).append(",");
            }
          }
          // Remove last appended , character.
          lineStringBuilder.setLength(lineStringBuilder.length() - 1);
          lineStringBuilder.append(LINE_SEPARATOR);
          proteinAlignmentWriter.write(lineStringBuilder.toString());
          lineStringBuilder = new StringBuilder();
        }
        // Close writer.
        proteinAlignmentWriter.close();
      }
    } catch (FileNotFoundException e) {
      throw new MusialIOException(e.getMessage());
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + proteinAlignmentOutFile.getAbsolutePath());
    } catch (CompoundNotFoundException e) {
      e.printStackTrace();
    }
  }

  /**
   * Computes a residue-residue contact map for a specified .pdb file (accessed via a referenceAnalysisId).
   *
   * @param outFile             {@link File} object specifying the output file.
   * @param referenceAnalysisId {@link String} used as a key in a {@link VariablePositionsTable} to specify the
   * @param arguments           {@link ArgumentsParser} arguments parsed from the command line.
   * @throws MusialIOException If any output file does not exist or could not be accessed.
   */
  public static void writeProteinContactMap(File outFile, String referenceAnalysisId, ArgumentsParser arguments)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      // Initialize output file.
      FileWriter proteinContactWriter = new FileWriter(outFile);
      File pdbFile = arguments.getPdInputFiles().get(referenceAnalysisId);
      Structure pdbStructure = getStructureFromPdb(pdbFile);
      StringBuilder lineStringBuilder = new StringBuilder();
      // Extract list of residues from .pdb file:
      HashMap<String, Group> residues = Bio.getResiduesFromStructure(pdbStructure);
      // Write header of output file.
      lineStringBuilder.append("\t").append(String.join("\t", residues.keySet())).append(LINE_SEPARATOR);
      proteinContactWriter.write(lineStringBuilder.toString());
      lineStringBuilder = new StringBuilder();
      // Compute list of residues sidechain centers of mass.
      HashMap<String, ArrayList<Double>> residuesSidechainCenterOfMass = Bio.getSidechainCentersOfMass(residues);
      // For each residue we have to iterate over each other residue in order to compute the distance in 3D space.
      for (String residueIdentifier : residues.keySet()) {
        lineStringBuilder.append(residueIdentifier).append("\t");
        for (String rI2 : residues.keySet()) {
          if (residueIdentifier.equals(rI2)) {
            lineStringBuilder.append(0.0).append("\t");
          } else {
            double distance3D = Math.sqrt(
                Math.pow(residuesSidechainCenterOfMass.get(residueIdentifier).get(0) -
                        residuesSidechainCenterOfMass.get(rI2).get(0)
                    , 2) +
                    Math.pow(residuesSidechainCenterOfMass.get(residueIdentifier).get(1) -
                        residuesSidechainCenterOfMass.get(rI2).get(1), 2) +
                    Math.pow(residuesSidechainCenterOfMass.get(residueIdentifier).get(2) -
                        residuesSidechainCenterOfMass.get(rI2).get(2), 2));
            lineStringBuilder.append(distance3D).append("\t");
          }
        }
        lineStringBuilder.setLength(lineStringBuilder.length() - 1);
        lineStringBuilder.append(LINE_SEPARATOR);
        proteinContactWriter.write(lineStringBuilder.toString());
        lineStringBuilder = new StringBuilder();
      }
      proteinContactWriter.close();
    } catch (FileNotFoundException e) {
      throw new MusialIOException(e.getMessage());
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }
}