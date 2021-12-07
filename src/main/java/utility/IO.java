package utility;

import static org.forester.util.ForesterUtil.round;

import com.google.common.base.Splitter;
import components.SnpEffAnnotator;
import datastructure.FastaEntry;
import datastructure.FeatureAnalysisEntry;
import datastructure.PositionStatistics;
import datastructure.SampleStatistics;
import datastructure.VariantContent;
import datastructure.VariantPositionsTable;
import exceptions.MusialBioException;
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
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.concurrent.ConcurrentSkipListSet;
import main.Musial;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.checkerframework.checker.units.qual.A;
import org.json.simple.JSONObject;
import components.ArgumentsParser;

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
        currHeader = currLine.replace(">", "");
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
   *
   * @param outFile               {@link File} instance pointing to the output file.
   * @param referenceAnalysisId   {@link String} internal reference feature identifier of which the output is written.
   * @param variantPositionsTable {@link VariantPositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeVariantsContentTable(File outFile, String referenceAnalysisId,
                                               VariantPositionsTable variantPositionsTable)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter snvTableWriter = new FileWriter(outFile.getAbsolutePath());
      StringBuilder lineBuilder = new StringBuilder();
      // Write header with meta-information.
      lineBuilder.append("#date=")
          .append(DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss").format(LocalDateTime.now())).append(";source=")
          .append(Musial.CLASS_NAME).append(Musial.VERSION).append(";reference=").append(referenceAnalysisId)
          .append(";samples=").append(variantPositionsTable.getSampleNames(referenceAnalysisId).size()).append(";" +
          "variantPositions=").append(variantPositionsTable.getVariantPositionsSet(referenceAnalysisId).size())
          .append(LINE_SEPARATOR);
      snvTableWriter.write(lineBuilder.toString());
      lineBuilder = new StringBuilder();
      // Write column headers.
      lineBuilder.append("Position\t").append("Reference\t").append(String.join("\t",
          variantPositionsTable.getSampleNames(referenceAnalysisId))).append(LINE_SEPARATOR);
      snvTableWriter.write(lineBuilder.toString());
      // Iterate over variant positions and write tab delimited content.
      Iterator<String> variantPositions = variantPositionsTable.getVariantPositions(referenceAnalysisId);
      HashSet<String> sampleNames = variantPositionsTable.getSampleNames(referenceAnalysisId);
      while (variantPositions.hasNext()) {
        String variantPosition = variantPositions.next();
        lineBuilder = new StringBuilder();
        lineBuilder.append(variantPosition).append("\t");
        // Append content of reference.
        if (variantPosition.contains("+")) {
          lineBuilder.append(VariantContent.INSERTION_DUMMY).append("\t");
        } else {
          Iterator<VariantContent> referenceContentIterator = variantPositionsTable.getContent(referenceAnalysisId,
              "Reference", variantPosition);
          if (referenceContentIterator != null) {
            lineBuilder.append(referenceContentIterator.next().content).append("\t");
          }
        }
        // Append content of each sample.
        for (String sampleName : sampleNames) {
          Iterator<VariantContent> sampleContentIterator =
              variantPositionsTable.getContent(referenceAnalysisId, sampleName, variantPosition);
          if (sampleContentIterator == null) {
            if (variantPosition.contains("I")) {
              lineBuilder.append(VariantContent.INSERTION_DUMMY).append("\t");
            } else {
              lineBuilder.append(VariantContent.REFERENCE).append("\t");
            }
          } else {
            StringBuilder sampleContentBuilder = new StringBuilder();
            while (sampleContentIterator.hasNext()) {
              VariantContent sampleContent = sampleContentIterator.next();
              sampleContentBuilder.append(sampleContent.content);
              sampleContentBuilder.append(VariantContent.SEPARATOR);
            }
            sampleContentBuilder.setLength(sampleContentBuilder.length() - 1);
            lineBuilder.append(sampleContentBuilder).append("\t");
          }
        }
        lineBuilder.setLength(lineBuilder.length() - 1);
        lineBuilder.append(LINE_SEPARATOR);
        snvTableWriter.write(lineBuilder.toString());
      }
      snvTableWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Writes output `.tsv` file containing annotations in a tab delimited format.
   *
   * @param outFile               {@link File} instance pointing to the output file.
   * @param referenceAnalysisId   {@link String} internal reference feature identifier of which the output is written.
   * @param variantPositionsTable {@link VariantPositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeVariantsAnnotationTable(File outFile, String referenceAnalysisId,
                                                  VariantPositionsTable variantPositionsTable)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter snvAnnotationsWriter = new FileWriter(outFile.getAbsolutePath());
      StringBuilder lineBuilder = new StringBuilder();
      // Write header with meta-information.
      lineBuilder.append("#date=")
          .append(DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss").format(LocalDateTime.now())).append(";source=")
          .append(Musial.CLASS_NAME).append(Musial.VERSION).append(";reference=").append(referenceAnalysisId)
          .append(";samples=").append(variantPositionsTable.getSampleNames(referenceAnalysisId).size()).append(";" +
          "variantPositions=").append(variantPositionsTable.getVariantPositionsSet(referenceAnalysisId).size())
          .append(LINE_SEPARATOR);
      snvAnnotationsWriter.write(lineBuilder.toString());
      lineBuilder = new StringBuilder();
      // Write column headers.
      lineBuilder.append("Position\t").append(String.join("\t",
          variantPositionsTable.getSampleNames(referenceAnalysisId))).append(LINE_SEPARATOR);
      snvAnnotationsWriter.write(lineBuilder.toString());
      // Iterate over variant positions and write tab delimited content.
      Iterator<String> variantPositions = variantPositionsTable.getVariantPositions(referenceAnalysisId);
      HashSet<String> sampleNames = variantPositionsTable.getSampleNames(referenceAnalysisId);
      while (variantPositions.hasNext()) {
        String variantPosition = variantPositions.next();
        lineBuilder = new StringBuilder();
        lineBuilder.append(variantPosition).append("\t");
        // Append annotation of each sample.
        for (String sampleName : sampleNames) {
          Iterator<VariantContent> sampleContentIterator =
              variantPositionsTable.getContent(referenceAnalysisId, sampleName, variantPosition);
          if (sampleContentIterator != null) {
            StringBuilder sampleContentBuilder = new StringBuilder();
            while (sampleContentIterator.hasNext()) {
              VariantContent sampleContent = sampleContentIterator.next();
              sampleContentBuilder.append("QUALITY=");
              sampleContentBuilder.append(round(sampleContent.quality, 2));
              sampleContentBuilder.append(";");
              sampleContentBuilder.append("COVERAGE=");
              sampleContentBuilder.append(round(sampleContent.coverage, 2));
              sampleContentBuilder.append(";");
              sampleContentBuilder.append("FREQUENCY=");
              sampleContentBuilder.append(round(sampleContent.frequency, 2));
              sampleContentBuilder.append(";");
              if (sampleContent.annotations.size() > 0) {
                for (Map.Entry<String, String> entry : sampleContent.annotations.entrySet()) {
                  sampleContentBuilder.append(entry.getKey()).append("=").append(entry.getValue()).append(";");
                }
              }
              sampleContentBuilder.setLength(sampleContentBuilder.length() - 1);
              sampleContentBuilder.append(VariantContent.SEPARATOR);
            }
            sampleContentBuilder.setLength(sampleContentBuilder.length() - 1);
            lineBuilder.append(sampleContentBuilder).append("\t");
          } else {
            lineBuilder.append("\t");
          }
        }
        lineBuilder.setLength(lineBuilder.length() - 1);
        lineBuilder.append(LINE_SEPARATOR);
        snvAnnotationsWriter.write(lineBuilder.toString());
      }
      snvAnnotationsWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Writes output `.fasta` file containing the full length sequences of all samples and the reference.
   *
   * @param outFile               {@link File} instance pointing to the output file.
   * @param referenceAnalysisId   {@link String} internal reference feature identifier of which the output is written.
   * @param variantPositionsTable {@link VariantPositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeSequences(File outFile, String referenceAnalysisId,
                                    VariantPositionsTable variantPositionsTable)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter sequenceWriter = new FileWriter(outFile.getAbsolutePath());
      sequenceWriter.write(">Reference|" + referenceAnalysisId + LINE_SEPARATOR);
      sequenceWriter.write(String.join(LINE_SEPARATOR,
          Splitter.fixedLength(80).split(Bio.getSequenceFromTable(variantPositionsTable, referenceAnalysisId,
              "Reference",
              variantPositionsTable.getAllPositionsIterator(referenceAnalysisId), false, true))));
      sequenceWriter.write(LINE_SEPARATOR);
      for (String sampleName : variantPositionsTable.getSampleNames(referenceAnalysisId)) {
        sequenceWriter.write(">Sample|" + sampleName + LINE_SEPARATOR);
        sequenceWriter.write(String.join(LINE_SEPARATOR,
            Splitter.fixedLength(80).split(Bio.getSequenceFromTable(variantPositionsTable, referenceAnalysisId,
                sampleName, variantPositionsTable.getAllPositionsIterator(referenceAnalysisId), false, true))));
        sequenceWriter.write(LINE_SEPARATOR);
      }
      sequenceWriter.close();
    } catch (IOException e) {
      throw new MusialIOException(
          "Failed to write to output file:\t" + outFile.getAbsolutePath() + ", Caused by: " + e.getMessage());
    }
  }

  /**
   * Writes a `.tsv` file containing per sample statistics.
   *
   * @param outFile               {@link File} instance pointing to the output file.
   * @param referenceAnalysisId   {@link String} internal reference feature identifier of which the output is written.
   * @param variantPositionsTable {@link VariantPositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writeSampleStatistics(File outFile, String referenceAnalysisId,
                                           VariantPositionsTable variantPositionsTable) throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter sampleStatisticsWriter = new FileWriter(outFile.getAbsolutePath());
      StringBuilder lineBuilder = new StringBuilder();
      // Write header with meta-information.
      lineBuilder.append("#date=")
          .append(DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss").format(LocalDateTime.now())).append(";source=")
          .append(Musial.CLASS_NAME).append(Musial.VERSION).append(";reference=").append(referenceAnalysisId)
          .append(";samples=").append(variantPositionsTable.getSampleNames(referenceAnalysisId).size())
          .append(LINE_SEPARATOR);
      sampleStatisticsWriter.write(lineBuilder.toString());
      lineBuilder = new StringBuilder();
      // Write column headers.
      lineBuilder.append(SampleStatistics.getHeaderString());
      sampleStatisticsWriter.write(lineBuilder.toString());
      lineBuilder = new StringBuilder();
      // Write sample statistics.
      HashSet<String> sampleNames = variantPositionsTable.getSampleNames(referenceAnalysisId);
      for (String sampleName : sampleNames) {
        lineBuilder.append(sampleName).append("\t")
            .append(variantPositionsTable.getSampleStatistics(referenceAnalysisId, sampleName).getContentString());
        sampleStatisticsWriter.write(lineBuilder.toString());
        lineBuilder = new StringBuilder();
      }
      sampleStatisticsWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Writes a `.tsv` file containing per position statistics.
   * <p>
   * One row contains information about the counts of each possible variant content (see {@link VariantContent})
   *
   * @param outFile               {@link File} instance pointing to the output file.
   * @param referenceAnalysisId   {@link String} internal reference feature identifier of which the output is written.
   * @param variantPositionsTable {@link VariantPositionsTable} instance containing variant information.
   * @throws MusialIOException If the output file could not be generated or not written to.
   */
  public static void writePositionStatistics(File outFile, String referenceAnalysisId,
                                             VariantPositionsTable variantPositionsTable) throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter positionStatisticsWriter = new FileWriter(outFile.getAbsolutePath());
      StringBuilder lineBuilder = new StringBuilder();
      // Write header with meta-information.
      lineBuilder.append("#date=")
          .append(DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss").format(LocalDateTime.now())).append(";source=")
          .append(Musial.CLASS_NAME).append(Musial.VERSION).append(";reference=").append(referenceAnalysisId)
          .append(";samples=").append(variantPositionsTable.getSampleNames(referenceAnalysisId).size())
          .append(LINE_SEPARATOR);
      positionStatisticsWriter.write(lineBuilder.toString());
      lineBuilder = new StringBuilder();
      // Write column headers.
      lineBuilder.append(PositionStatistics.getHeaderString());
      positionStatisticsWriter.write(lineBuilder.toString());
      lineBuilder = new StringBuilder();
      // Write sample statistics.
      Iterator<String> variantPositions = variantPositionsTable.getVariantPositions(referenceAnalysisId);
      while (variantPositions.hasNext()) {
        String variantPosition = variantPositions.next();
        String printablePosition;
        if (variantPosition.contains("I")) {
          printablePosition =
              variantPosition.replace("I", "") + "+" + variantPosition.chars().filter(c -> c == 'I').count();
        } else {
          printablePosition = variantPosition;
        }
        lineBuilder.append(printablePosition).append("\t");
        lineBuilder.append(
            variantPositionsTable.getPositionStatistics(referenceAnalysisId, variantPosition).getContentString());
        positionStatisticsWriter.write(lineBuilder.toString());
        lineBuilder = new StringBuilder();
      }
      positionStatisticsWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Writes an annotated summary of all detected variants into an .tsv file.
   *
   * @param outFile               {@link File} instance pointing to the output file.
   * @param referenceAnalysisId   {@link String} internal reference feature identifier of which the output is written.
   * @param variantPositionsTable {@link VariantPositionsTable} instance containing variant information.
   * @throws MusialIOException  If the output file could not be generated or not written to.
   * @throws MusialBioException If any position with no reference content entry is accessed.
   */
  public static void writeSnpEffSummary(File outFile, String referenceAnalysisId,
                                        VariantPositionsTable variantPositionsTable)
      throws MusialIOException, MusialBioException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      FileWriter snpEffSummaryWriter = new FileWriter(outFile.getAbsolutePath());
      StringBuilder lineBuilder = new StringBuilder();
      // Write header with meta-information.
      lineBuilder.append("#date=")
          .append(DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss").format(LocalDateTime.now())).append(";source=")
          .append(Musial.CLASS_NAME).append(Musial.VERSION).append(";reference=").append(referenceAnalysisId)
          .append(";samples=").append(variantPositionsTable.getSampleNames(referenceAnalysisId).size())
          .append(LINE_SEPARATOR);
      snpEffSummaryWriter.write(lineBuilder.toString());
      lineBuilder = new StringBuilder();
      // Write column headers.
      lineBuilder.append(SnpEffAnnotator.snpEffSummaryHeader);
      snpEffSummaryWriter.write(lineBuilder.toString());
      lineBuilder = new StringBuilder();
      // Collect and write annotation summary.
      ConcurrentSkipListSet<String> variantPositions =
          variantPositionsTable.getVariantPositionsSet(referenceAnalysisId);
      HashMap<String, ArrayList<String>> perAnnotationSamples;
      ArrayList<String> annotationSampleNames;
      HashSet<String> sampleNames = variantPositionsTable.getSampleNames(referenceAnalysisId);
      Iterator<VariantContent> variantContentIterator;
      VariantContent variantContent;
      String snpEffAnnotation;
      Iterator<VariantContent> referenceContentIterator;
      char referenceContent;
      for (String variantPosition : variantPositions) {
        perAnnotationSamples = new HashMap<>();
        // Skip any insertion positions that exceed the +1 suffix to avoid duplications.
        if (variantPosition.contains("+") && !variantPosition.split("\\+")[1].equals("1")) {
          continue;
        }
        for (String sampleName : sampleNames) {
          variantContentIterator = variantPositionsTable.getContent(referenceAnalysisId, sampleName, variantPosition);
          if (variantContentIterator != null && variantContentIterator.hasNext()) {
            while (variantContentIterator.hasNext()) {
              variantContent = variantContentIterator.next();
              if (variantContent.annotations.containsKey(SnpEffAnnotator.snpEffAnnotationAttributeKey)) {
                snpEffAnnotation = variantContent.annotations.get(SnpEffAnnotator.snpEffAnnotationAttributeKey);
                if (perAnnotationSamples.containsKey(snpEffAnnotation)) {
                  perAnnotationSamples.get(snpEffAnnotation).add(sampleName);
                } else {
                  annotationSampleNames = new ArrayList<>();
                  annotationSampleNames.add(sampleName);
                  perAnnotationSamples.put(snpEffAnnotation, annotationSampleNames);
                }
              }
            }
          }
        }
        for (Map.Entry<String, ArrayList<String>> entry : perAnnotationSamples.entrySet()) {
          // If the position reflects an insertion only the prefix of the position string is used.
          if (variantPosition.contains("+")) {
            variantPosition = variantPosition.split("\\+")[0];
          }
          referenceContentIterator = variantPositionsTable.getContent(referenceAnalysisId, "Reference",
              variantPosition);
          if (referenceContentIterator != null && referenceContentIterator.hasNext()) {
            referenceContent =
                referenceContentIterator.next().content;
          } else {
            throw new MusialBioException("Tried to access reference content at position " + variantPosition + " of " +
                "feature " + referenceAnalysisId + ", but no entry exits.");
          }
          lineBuilder.append(variantPosition).append("\t").append(referenceContent).append("\t")
              .append(String.join(",", entry.getValue())).append("\t").append(String.join("\t", entry.getKey().split(
              "\\|"))).append(IO.LINE_SEPARATOR);
          snpEffSummaryWriter.write(lineBuilder.toString());
          lineBuilder = new StringBuilder();
        }
      }
      snpEffSummaryWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

  /**
   * Allocates sample and reference data to provided protein data and writes the allocation into a JSON file.
   *
   * @param outFile               {@link File} instance pointing to the output file.
   * @param featureAnalysisEntry  {@link FeatureAnalysisEntry} specifying the reference feature to analyze.
   * @param variantPositionsTable {@link VariantPositionsTable} storing the corresponding data.
   * @param arguments             {@link ArgumentsParser} arguments parsed from the command line.
   * @throws MusialIOException If any output file does not exist or could not be accessed.
   */
  @SuppressWarnings("unchecked")
  public static void writeStructureMatching(File outFile, FeatureAnalysisEntry featureAnalysisEntry,
                                            VariantPositionsTable variantPositionsTable, ArgumentsParser arguments)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      // Access the reference nucleotide sequence; the reverse complement is built, if the reference is located on
      // the ani-sense strand.
      String referenceAnalysisId = featureAnalysisEntry.name;
      String referenceSequence = Bio.getSequenceFromTable(variantPositionsTable, referenceAnalysisId, "Reference",
          variantPositionsTable.getAllPositionsIterator(referenceAnalysisId), false, false);
      if (!featureAnalysisEntry.isSense) {
        referenceSequence = Bio.getReverseComplement(referenceSequence);
      }
      // Access the amino-acid sequence from the .pdb file specified to hold the structure information related
      // to the reference feature.
      File pdbFile = arguments.getPdInputFiles().get(referenceAnalysisId);
      Structure pdbStructure = getStructureFromPdb(pdbFile);
      // (1) Generate JSONObject and add meta-information.
      JSONObject allocationJson = new JSONObject();
      allocationJson.put("Date", DateTimeFormatter.ofPattern("yyyy/MM/dd").format(LocalDateTime.now()));
      allocationJson.put("Software", Musial.CLASS_NAME + Musial.VERSION);
      allocationJson.put("GeneReferenceName", featureAnalysisEntry.name);
      allocationJson.put("GeneReferenceIsSense", featureAnalysisEntry.isSense);
      allocationJson.put("GeneReference5'Position", featureAnalysisEntry.locationStart);
      allocationJson.put("GeneReference3'Position", featureAnalysisEntry.locationEnd);
      allocationJson.put("GeneReferenceLocus", featureAnalysisEntry.referenceSequenceLocation);
      allocationJson.put("GeneReferenceLength", referenceSequence.length());
      allocationJson.put("ProteinReferenceFile", pdbFile.getName());

      HashMap<String, JSONObject> perChainAllocation = new HashMap<>();
      // Each chain of the protein is considered separately.
      for (Chain chain : pdbStructure.getChains()) {
        String chainId = chain.getName();
        String chainAsymId = chain.getId();
        // (1) Generate JSONObject with allocation data.
        JSONObject chainAllocationData = Bio.matchDataToStructure(variantPositionsTable, referenceAnalysisId,
            pdbStructure, chainId, chainAsymId, featureAnalysisEntry.isSense);
        // (2) Store chain allocation data in map.
        perChainAllocation.put(chainId, chainAllocationData);
      }
      for (Map.Entry<String, JSONObject> entry : perChainAllocation.entrySet()) {
        allocationJson.put("DataChain" + entry.getKey(), entry.getValue());
      }
      FileWriter structureAllocationWriter = new FileWriter(outFile.getAbsolutePath());
      structureAllocationWriter.write(allocationJson.toJSONString());
      structureAllocationWriter.close();
    } catch (IOException | MusialBioException e) {
      throw new MusialIOException(
          "Failed to write to output file:\t" +
              outFile.getAbsolutePath() + ", Caused by: " + e.getMessage());
    }
  }

  /**
   * Computes a residue-residue distance map for a specified .pdb file (accessed via a referenceAnalysisId).
   *
   * @param outFile             {@link File} object specifying the output file.
   * @param referenceAnalysisId {@link String} internal reference feature identifier of which the output is written.
   * @param arguments           {@link ArgumentsParser} arguments parsed from the command line.
   * @throws MusialIOException If any output file does not exist or could not be accessed.
   */
  public static void writeProteinDistanceMap(File outFile, String referenceAnalysisId, ArgumentsParser arguments)
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
      lineStringBuilder.append("#date=")
          .append(DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss").format(LocalDateTime.now())).append(";");
      lineStringBuilder.append("source=").append(Musial.CLASS_NAME).append(Musial.VERSION).append(";");
      lineStringBuilder.append("protein=").append(pdbFile.getName()).append(LINE_SEPARATOR);
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

  /**
   * Writes run information to a text file.
   *
   * @param outFile             {@link File} object specifying the output file.
   * @param arguments           {@link ArgumentsParser} arguments parsed from the command line.
   * @param referenceAnalysisId {@link String} internal reference feature identifier of which the output is written.
   * @throws MusialIOException If any output file does not exist or could not be accessed.
   */
  public static void writeRunInfo(File outFile, ArgumentsParser arguments, String referenceAnalysisId)
      throws MusialIOException {
    if (!outFile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outFile.getAbsolutePath());
    }
    try {
      // Initialize output file.
      FileWriter runInfoWriter = new FileWriter(outFile);
      runInfoWriter.write("METAINFO" + LINE_SEPARATOR);
      runInfoWriter.write("VERSION=" + Musial.CLASS_NAME + Musial.VERSION + LINE_SEPARATOR);
      runInfoWriter.write(
          "DATE=" + DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss").format(LocalDateTime.now()) + LINE_SEPARATOR);
      runInfoWriter.write("ARGUMENTS" + LINE_SEPARATOR);
      runInfoWriter.write("COMMAND=" + String.join(" ", arguments.getArguments()));
      runInfoWriter.write("MIN_QUALITY=" + arguments.getMinQuality() + LINE_SEPARATOR);
      runInfoWriter.write("MIN_COVERAGE=" + arguments.getMinCoverage() + LINE_SEPARATOR);
      runInfoWriter.write("MIN_HOM_FREQUENCY=" + arguments.getMinFrequency() + LINE_SEPARATOR);
      runInfoWriter.write("MIN_HET_FREQUENCY=" + arguments.getMinHet() + LINE_SEPARATOR);
      runInfoWriter.write("MAX_HET_FREQUENCY=" + arguments.getMaxHet() + LINE_SEPARATOR);
      runInfoWriter.write("INPUT" + LINE_SEPARATOR);
      runInfoWriter.write("REFERENCE=" + arguments.getReferenceFile().getAbsolutePath() + LINE_SEPARATOR);
      runInfoWriter.write("FEATURE=" + referenceAnalysisId + LINE_SEPARATOR);
      if (arguments.getPdInputFiles().containsKey(referenceAnalysisId)) {
        runInfoWriter.write(
            "PROTEIN_DATA=" + arguments.getPdInputFiles().get(referenceAnalysisId).getAbsolutePath() + LINE_SEPARATOR);
      }
      for (int i = 0; i < arguments.getSampleInput().size(); i++) {
        File sampleFile = arguments.getSampleInput().get(i);
        runInfoWriter.write("SAMPLE_" + (i + 1) + "=" + sampleFile.getAbsolutePath() + LINE_SEPARATOR);
      }
      runInfoWriter.close();
    } catch (FileNotFoundException e) {
      throw new MusialIOException(e.getMessage());
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outFile.getAbsolutePath());
    }
  }

}