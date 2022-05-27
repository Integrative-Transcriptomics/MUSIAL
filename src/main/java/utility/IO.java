package utility;

import com.google.common.base.Splitter;
import com.google.gson.Gson;
import datastructure.FastaContainer;
import datastructure.VariantsDictionary;
import exceptions.MusialIOException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Scanner;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.zip.GZIPInputStream;
import main.Musial;
import org.apache.commons.io.FileUtils;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;

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
  public static ArrayList<String> readLinesFromFile(String filePath) throws FileNotFoundException {
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
   * Reads a .fasta file from a {@link File} object into a {@link HashSet} of {@link FastaContainer} instances. Each
   * element contains the header and sequence of the respective entry from the .fasta file.
   *
   * @param file {@link File} specifying a .fasta file.
   * @return {@link HashSet} of {@link FastaContainer} instances.
   * @throws IOException If the respective file can not be found or accessed.
   */
  public static HashSet<FastaContainer> readFastaToSet(File file) throws IOException {
    HashSet<FastaContainer> fastaEntries = new HashSet<>();
    @SuppressWarnings("resource")
    BufferedReader br = new BufferedReader(new FileReader(file));
    String currLine;
    String currHeader = "";
    StringBuilder currSequence = new StringBuilder();
    while ((currLine = br.readLine()) != null) {
      if (currLine.startsWith(">")) {
        if (currHeader.length() > 0) {
          fastaEntries.add(new FastaContainer(currHeader, currSequence.toString()));
          currSequence = new StringBuilder();
        }
        currHeader = currLine.replace(">", "");
      } else if (!currLine.startsWith(";")) {
        currSequence.append(currLine.trim());
      }
    }
    fastaEntries.add(new FastaContainer(currHeader, currSequence.toString()));
    return fastaEntries;
  }

  public static ArrayList<FastaContainer> readFastaToList(File file) throws IOException {
    ArrayList<FastaContainer> fastaEntries = new ArrayList<>();
    fastaEntries.addAll(readFastaToSet(file));
    return fastaEntries;
  }

  public static VariantsDictionary readVariantsDictionary(File vDictFile) throws IOException {
    assert vDictFile.exists();
    // Retrieve JSON String from gzip compressed DB dump.
    try (
        BufferedReader bufferedReader = new BufferedReader(
            new InputStreamReader(Files.newInputStream(vDictFile.toPath()),
                StandardCharsets.UTF_8));
    ) {
      Gson gson = new Gson();
      VariantsDictionary vDict = gson.fromJson(bufferedReader, VariantsDictionary.class);
      vDict.novelVariants = new ConcurrentSkipListSet<>((s1, s2) -> {
        int p1 = Integer.parseInt(s1.split("@")[0]);
        int p2 = Integer.parseInt(s2.split("@")[0]);
        return Integer.compare(p1, p2);
      });
      return vDict;
    }
  }

  /**
   * Accepts a {@link File} object representing a directory that is subject to deletion.
   *
   * @param file {@link File} specifying the directory to delete.
   * @throws IOException If any file deletion procedure failed, for example the specified directory does not exist.
   */
  public static void deleteDirectory(File file) throws IOException {
    FileUtils.deleteDirectory(file);
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
        throw new MusialIOException("Failed to generate output file:\t" + file.getAbsolutePath());
      }
    } catch (IOException e) {
      throw new MusialIOException(e.getMessage());
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
   * Returns the chains sequences from a {@link Structure} instance parsed from a .pdb file.
   *
   * @param structure {@link Structure} instance to extract the chains sequences from.
   * @return {@link HashMap} storing the .pdb files amino-acid sequences per chain.
   */
  public static HashMap<String, String> getSequencesFromPdbStructure(Structure structure) {
    // Initialize results list.
    HashMap<String, String> pdbChainsSequences = new HashMap<>();
    // Parse nucleotide information from chains.
    int modelNr = 0;
    List<Chain> chains = structure.getChains(modelNr);
    for (Chain chain : chains) {
      // Skip chains representing membrane.
      if (chain.getName().equals("x")) {
        continue;
      }
      pdbChainsSequences.put(chain.getName(), chain.getAtomSequence());
    }
    return pdbChainsSequences;
  }

  /**
   * TODO
   *
   * @param pdbFile
   * @return
   * @throws IOException
   */
  public static Structure readStructure(File pdbFile) throws IOException {
    try {
      System.setErr(Musial.EMPTY_STREAM);
      Structure pdbStructure = new PDBFileReader().getStructure(pdbFile);
      System.setErr(Musial.ORIGINAL_ERR_STREAM);
      return pdbStructure;
    } finally {
      System.setErr(Musial.ORIGINAL_ERR_STREAM);
    }
  }

  public static void writeFasta(File outputFile, HashMap<String, ArrayList<String>> sequences) {
    try {
      FileWriter writer = new FileWriter(outputFile);
      for (Map.Entry<String, ArrayList<String>> sequenceEntry : sequences.entrySet()) {
        writer.write(">" + String.join("|", sequenceEntry.getValue()) + IO.LINE_SEPARATOR);
        for (String l : Splitter.fixedLength(80).split(sequenceEntry.getKey())) {
          writer.write(l + IO.LINE_SEPARATOR);
        }
        writer.flush();
      }
      writer.close();
    } catch (IOException e) {
      Logging.logWarning("Failed to write `FASTA` format file to " + outputFile.getAbsolutePath());
    }
  }

  public static void writeVcf(File outputFile, NavigableSet<String> variants, String chrom) {
    try {
      FileWriter writer = new FileWriter(outputFile);
      writer.write("##fileformat=VCFv4.2" + IO.LINE_SEPARATOR);
      writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + IO.LINE_SEPARATOR);
      for (String variant : variants) {
        String[] variantFields = variant.split("@");
        writer.write(chrom + "\t"
            + variantFields[0] + "\t"
            + ".\t"
            + variantFields[1] + "\t"
            + variantFields[2].replace("-", "") + "\t"
            + "1000\t"
            + ".\t"
            + "\t"
            + IO.LINE_SEPARATOR);
        writer.flush();
      }
      writer.close();
    } catch (IOException e) {
      Logging.logWarning("Failed to write `VCF` format file to " + outputFile.getAbsolutePath());
    }
  }

}