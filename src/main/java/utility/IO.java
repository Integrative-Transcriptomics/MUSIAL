package utility;

import com.google.common.collect.Lists;
import datastructure.FastaEntry;
import datastructure.GeneFeature;
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
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.Scanner;
import main.Musial;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
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
   * Directory suffix of a EAGER output folder pointing to the used .vcf file
   */
  private static final String EAGER_VCF_SUFFIX = "/10-GATKGenotype";
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
   * Reads a `.pdb` file and returns the therein stored proteins amino-acid sequence.
   *
   * @param file {@link File} object pointing to a `.pdb` file.
   * @return {@link String} representing the amino acid sequence of the protein stored in the passed `.pdb` file.
   * @throws MusialIOException If the specified `.pdb` file does not exist or could not be accessed.
   */
  public static String getSequenceFromPdb(File file) throws MusialIOException {
    // Initialize pdbReader.
    PDBFileReader pdbReader = new PDBFileReader();
    FileParsingParameters params = new FileParsingParameters();
    params.setParseSecStruc(true);
    params.setAlignSeqRes(true);
    pdbReader.setFileParsingParameters(params);
    // Initialize string builders for amino-acid sequence.
    StringBuilder aaSequenceStringBuilder = new StringBuilder();
    // Read the `.pdb` file.
    Structure pdbStructure;
    try {
      pdbStructure = pdbReader.getStructure(file);
    } catch (IOException e) {
      throw new MusialIOException("Failed to read file:\t" + file.getAbsolutePath());
    }
    // Parse nucleotide information from chains.
    int modelNr = 0;
    List<Chain> chains = pdbStructure.getChains(modelNr);
    for (Chain chain : chains) {
      aaSequenceStringBuilder.append(chain.getAtomSequence());
    }
    return aaSequenceStringBuilder.toString();
  }

  /**
   * Reads a `.pdb` file and returns the file content as a string.
   *
   * @param file {@link File} object pointing to a `.pdb` file.
   * @return {@link String} representing the file content of the specified `.pdb` file.
   * @throws FileNotFoundException If the specified `.pdb` file does not exist.
   */
  public static String getPdbAsString(File file) throws FileNotFoundException {
    StringBuilder pdbStringBuilder = new StringBuilder();
    List<String> pdbLines = IO.getLinesFromFile(file.getAbsolutePath());
    for (int i = 0; i < pdbLines.size(); i++) {
      pdbStringBuilder.append(pdbLines.get(i));
      if (i != (pdbLines.size() - 1)) {
        pdbStringBuilder.append(IO.LINE_SEPARATOR);
      }
    }
    return pdbStringBuilder.toString();
  }

  /**
   * Scans an EAGER output directory for .vcf files and returns, for each sample, the sample name and .vcf files path.
   *
   * @param outputDirectory {@link String} specifying an EAGER output directory.
   * @return A size two array containing two {@link ArrayList} objects, the first comprising the .vcf file paths and
   * the second comprising the sample names.
   */
  public static ArrayList<ArrayList<String>> getInputFromEagerOutput(String outputDirectory) {
    ArrayList<String> vcfFilePaths = new ArrayList<>();
    ArrayList<String> sampleNameList = new ArrayList<>();
    File eagerOutput = new File(outputDirectory);
    for (File sampleDir : Objects.requireNonNull(eagerOutput.listFiles())) {
      if (sampleDir.isDirectory()) {
        File vcfDirectory = new File(sampleDir.getAbsolutePath() + EAGER_VCF_SUFFIX);
        sampleNameList.add(sampleDir.getName());
        for (File vcfFile : Objects.requireNonNull(vcfDirectory.listFiles())) {
          if (vcfFile.getName().endsWith(".vcf") || vcfFile.getName().endsWith(".vcf.gz")) {
            vcfFilePaths.add(vcfFile.getName());
            break;
          }
        }
      }
    }
    ArrayList<ArrayList<String>> result = new ArrayList<>();
    result.add(vcfFilePaths);
    result.add(sampleNameList);
    return result;
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
   * Writes output `.tsv` file containing nucleotide variants in a tab delimited format.
   * <p>
   * The expected output format is:
   * Position \t Reference (optional) \t [SAMPLE_1_NAME] \t ...
   * x (Integer) \t [Nucleotide at position x in reference] (optional) \t [Nucleotide at position x in sample 1] \t ...
   *
   * @param outFile                {@link File} instance pointing to the output file.
   * @param referenceAnalysisId    {@link String} internal reference feature identifier of which the output is written.
   * @param variablePositionsTable {@link VariablePositionsTable} instance containing variant information.
   * @param addReference           {@link Boolean} whether reference information should be added.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeSnvTable(File outFile, String referenceAnalysisId,
                                   VariablePositionsTable variablePositionsTable, boolean addReference)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter snvTableWriter = new FileWriter(outFile.getAbsolutePath());
      snvTableWriter.write(
          variablePositionsTable.getSampleNamesHeaderTabDelimited(referenceAnalysisId, addReference) + LINE_SEPARATOR
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
                .getPositionContentTabDelimited(referenceAnalysisId, variantPosition, shift, addReference) +
                LINE_SEPARATOR
        );
      }
      snvTableWriter.close();
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
   * @param addReference           {@link Boolean} whether reference information should be added.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeSequenceAlignment(File outFile, String referenceAnalysisId,
                                            VariablePositionsTable variablePositionsTable, boolean addReference)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter sequenceAlignmentWriter = new FileWriter(outFile.getAbsolutePath());
      if (addReference) {
        sequenceAlignmentWriter.write(variablePositionsTable.getSampleFullSequence(referenceAnalysisId,
            "Reference", true));
      }
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
   * Writes output `.fasta` file containing the sequence of variants of all samples and the reference (optional).
   *
   * @param outFile                {@link File} instance pointing to the output file.
   * @param referenceAnalysisId    {@link String} internal reference feature identifier of which the output is written.
   * @param variablePositionsTable {@link VariablePositionsTable} instance containing variant information.
   * @param addReference           {@link Boolean} whether reference information should be added.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeVariantAlignment(File outFile, String referenceAnalysisId,
                                           VariablePositionsTable variablePositionsTable, boolean addReference)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter variantAlignmentWriter = new FileWriter(outFile.getAbsolutePath());
      if (addReference) {
        variantAlignmentWriter.write(variablePositionsTable.getSampleVariantsSequence(referenceAnalysisId,
            "Reference", true));
      }
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
   * Writes a `.tsv` file containing per sample statistics.
   *
   * @param outFile                {@link File} instance pointing to the output file.
   * @param referenceAnalysisId    {@link String} internal reference feature identifier of which the output is written.
   * @param variablePositionsTable {@link VariablePositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  @SuppressWarnings("unused")
  public static void writeSampleStatistics(File outFile, String referenceAnalysisId,
                                           VariablePositionsTable variablePositionsTable) throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter sampleStatisticsWriter = new FileWriter(outFile.getAbsolutePath());
      /*
        TODO
       */
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
          variablePositionsTable.getPositionStatisticsHeaderTabDelimited() + LINE_SEPARATOR
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
          variablePositionsTable.getSampleNamesHeaderTabDelimited(referenceAnalysisId, false) + LINE_SEPARATOR
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
   * Writes a `.json` file storing the data for visualization with MUSIAL IVE.
   * <p>
   * In the current version 2.0 the file contains the following data:
   * - The name of the reference feature (e.g. a gene name).
   * - The reference nucleotide sequence.
   * - Per sample per (variant) position variant bases.
   * - Per sample per (variant) position annotations (i.e. coverage and quality).
   * - (Optional) The name of the protein corresponding to the reference feature, currently the name of the `.pdb` file.
   * - (Optional) The amino acid sequence of the protein corresponding to the reference feature.
   *
   * @param outFile                {@link File} object specifying the output file.
   * @param referenceAnalysisId    {@link String} used as a key in a {@link VariablePositionsTable} to specify the
   *                               reference feature for which the config file is generated.
   * @param variablePositionsTable {@link VariablePositionsTable} storing the corresponding data.
   * @param arguments              {@link ArgumentsParser} arguments parsed from the command line.
   * @throws MusialIOException If any output file does not exist or could not be accessed.
   */
  @SuppressWarnings("unchecked")
  public static void writeIveConfig(File outFile, String referenceAnalysisId,
                                    VariablePositionsTable variablePositionsTable, ArgumentsParser arguments)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      // Initialize JSON object.
      JSONObject mainObj = new JSONObject();
      // If protein information is not added, simply add reference sequence information.
      String referenceSequence = variablePositionsTable.getSampleFullSequence(referenceAnalysisId,
          "Reference", false);
      if (!arguments.isIncludeStructureInformation()) {
        mainObj.put("ReferenceSequence", referenceSequence);
        mainObj.put("ReferenceName", referenceAnalysisId);
      } else {
        // Otherwise, first match protein and reference sequence before adding them to JSON.
        File pdbFile = arguments.getPdInputFiles().get(referenceAnalysisId);
        String proteinStructure = getPdbAsString(pdbFile);
        String proteinSequence = getSequenceFromPdb(pdbFile);
        ArrayList<String> matchedSequences = Bio.matchSequences(proteinSequence, referenceSequence);
        proteinSequence = matchedSequences.get(0);
        referenceSequence = matchedSequences.get(1);
        mainObj.put("ReferenceSequence", referenceSequence);
        mainObj.put("ReferenceSequenceName", referenceAnalysisId);
        mainObj.put("ProteinSequence", proteinSequence);
        mainObj.put("ProteinStructure", proteinStructure);
        mainObj.put("ProteinName",
            FilenameUtils.removeExtension(FilenameUtils.getBaseName(pdbFile.getAbsolutePath())));
      }
      // Add sample names.
      JSONArray sampleNames = new JSONArray();
      sampleNames.addAll(variablePositionsTable.getSampleNamesOf(referenceAnalysisId));
      mainObj.put("SampleNames", sampleNames);
      // For each variant position, add sample variants, annotation and overall counts.
      JSONObject perPositionVariants = new JSONObject();
      JSONObject perPositionAnnotations = new JSONObject();
      JSONObject perPositionCounts = new JSONObject();
      int shift = 0;
      // TODO: Temporary solution to get correct reference position offset/start.
      int locationOffset = 0;
      if (arguments.getIncludedGeneFeatures().size() > 0) {
        for (GeneFeature includedGeneFeature : arguments.getIncludedGeneFeatures()) {
          if (includedGeneFeature.featureName.equals(referenceAnalysisId)) {
            locationOffset = includedGeneFeature.startPosition;
          }
        }
      }
      for (Iterator<String> it = variablePositionsTable.getVariantPositions(referenceAnalysisId); it.hasNext(); ) {
        String variantPosition = it.next();
        int position = Integer.parseInt(variantPosition) - locationOffset + 1; // Add 1 due to 1-based indexing.
        if (variantPosition.contains("I")) {
          shift += 1;
        }
        JSONArray positionVariants = new JSONArray();
        JSONArray positionAnnotations = new JSONArray();
        JSONArray positionCounts = new JSONArray();
        ArrayList<String> positionVariantsList =
            Lists.newArrayList(variablePositionsTable
                .getPositionContentTabDelimited(referenceAnalysisId, variantPosition, shift, false).split("\t"));
        positionVariantsList.remove(0);
        ArrayList<String> positionAnnotationsList =
            Lists.newArrayList(variablePositionsTable.getPositionAnnotationTabDelimited(referenceAnalysisId,
                variantPosition, shift).split("\t"));
        positionAnnotationsList.remove(0);
        ArrayList<String> positionCountsList =
            Lists.newArrayList(variablePositionsTable.getPositionStatisticsTabDelimited(referenceAnalysisId,
                variantPosition, shift).split("\t"));
        positionCountsList.remove(0);
        positionVariants.addAll(positionVariantsList);
        positionAnnotations.addAll(positionAnnotationsList);
        positionCounts.addAll(positionCountsList);
        perPositionVariants.put(position, positionVariants);
        perPositionAnnotations.put(position, positionAnnotations);
        perPositionCounts.put(position, positionCounts);
      }
      mainObj.put("PerPositionVariants", perPositionVariants);
      mainObj.put("PerPositionAnnotations", perPositionAnnotations);
      mainObj.put("PerPositionCounts", perPositionCounts);
      // Write JSON to file.
      FileWriter iveConfigWriter = new FileWriter(outFile.getAbsolutePath());
      iveConfigWriter.write(mainObj.toString());
      iveConfigWriter.close();
    } catch (FileNotFoundException e) {
      throw new MusialIOException(e.getMessage());
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    } catch (CompoundNotFoundException e) {
      e.printStackTrace();
    }
  }


}
