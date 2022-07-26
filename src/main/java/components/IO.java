package components;

import com.google.common.base.Splitter;
import com.google.gson.Gson;
import datastructure.FastaContainer;
import datastructure.NucleotideVariantAnnotationEntry;
import datastructure.VariantsDictionary;
import exceptions.MusialIOException;
import main.Musial;
import org.apache.commons.io.FileUtils;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.*;
import java.util.concurrent.ConcurrentSkipListMap;

/**
 * This class comprises static methods used for reading and writing files.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
@SuppressWarnings("unused")
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

    /**
     * Initializes a {@link VariantsDictionary} object from a variants dictionary JSON file.
     *
     * @param vDictFile {@link File} pointing to the JSON file to parse.
     * @return {@link VariantsDictionary} instance.
     * @throws IOException If any error occurs while parsing the variants dictionary JSON file.
     */
    public static VariantsDictionary readVariantsDictionary(File vDictFile) throws IOException {
        assert vDictFile.exists();
        // Retrieve JSON String from gzip compressed DB dump.
        try (
                BufferedReader bufferedReader = new BufferedReader(
                        new InputStreamReader(Files.newInputStream(vDictFile.toPath()),
                                StandardCharsets.UTF_8))
        ) {
            Gson gson = new Gson();
            return gson.fromJson(bufferedReader, VariantsDictionary.class);
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
    @SuppressWarnings("unused")
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
     * Parses a {@link Structure} instance from a {@link File} pointing to a pdb file.
     *
     * @param pdbFile {@link File} instance pointing to a pdb file.
     * @return {@link Structure} object parsed from specified pdb file.
     * @throws IOException If any error occurs while parsing the pdb file.
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

    /**
     * Writes a fasta format file to the specified output file from a {@link HashMap} instance mapping sequences to
     * lists of identifiers (used to construct the fasta entry headers).
     * <p>
     * Each key of the passed map will be used to build one fasta entry. The header of the respective entry is
     * constructed by joining all {@link String}s of the value accessible via the (sequence) key with the `|` delimiter.
     *
     * @param outputFile {@link File} object pointing to the output fasta file.
     * @param sequences  {@link HashMap} mapping sequences to {@link ArrayList} of strings.
     */
    @SuppressWarnings("unused")
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

    /**
     *
     * @param outputFile
     * @param lineContent
     * @throws MusialIOException
     */
    public static void writeFile(File outputFile, ArrayList<String> lineContent) throws MusialIOException {
        if (!Validation.isFile(outputFile)) {
            throw new MusialIOException("Failed to write to file " + outputFile.getAbsolutePath() + "; File does not exist or has no write permission. ");
        }
        try (FileWriter writer = new FileWriter(outputFile)) {
            for (String line : lineContent) {
                writer.write(line);
                writer.flush();
            }
        } catch (IOException e) {
            Logging.logWarning("Failed to write content to file " + outputFile.getAbsolutePath());
        }
    }

    /**
     * Writes a dummy .vcf format file to the specified output file from a {@link VariantsDictionary} instance.
     *
     * @param outputFile       {@link File} object pointing to the output vcf file.
     * @param vDict            {@link VariantsDictionary} containing variant information.
     * @param excludedFeatures {@link HashSet} of {@link String}s; Internal feature names to be excluded.
     * @param excludedSamples  {@link HashSet} of {@link String}s; Internal sample names to be excluded.
     */
    public static void writeVcf(File outputFile, VariantsDictionary vDict, HashSet<String> excludedFeatures, HashSet<String> excludedSamples) {
        try {
            FileWriter writer = new FileWriter(outputFile);
            writer.write("##fileformat=VCFv4.2" + IO.LINE_SEPARATOR);
            writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + IO.LINE_SEPARATOR);
            int variantPosition;
            ConcurrentSkipListMap<String, NucleotideVariantAnnotationEntry> variants;
            boolean skip;
            for (Map.Entry<Integer, ConcurrentSkipListMap<String, NucleotideVariantAnnotationEntry>> variantPositionEntry : vDict.variants.entrySet()) {
                variantPosition = variantPositionEntry.getKey();
                skip = false;
                for (String excludedFeature : excludedFeatures) {
                    if (variantPosition >= vDict.features.get(excludedFeature).start && variantPosition <= vDict.features.get(excludedFeature).end) {
                        skip = true;
                        break;
                    }
                }
                if (skip) {
                    continue;
                }
                variants = variantPositionEntry.getValue();
                for (Map.Entry<String, NucleotideVariantAnnotationEntry> variant : variants.entrySet()) {
                    if (!excludedSamples.containsAll(variant.getValue().occurrence.keySet())) {
                        writer.write(vDict.chromosome + "\t"
                                + variantPosition + "\t"
                                + ".\t"
                                + variant.getValue().annotations.get(VariantsDictionary.ATTRIBUTE_VARIANT_REFERENCE_CONTENT) + "\t"
                                + variant.getKey().replace("-", "") + "\t"
                                + "1000\t"
                                + ".\t"
                                + "\t"
                                + IO.LINE_SEPARATOR);
                        writer.flush();
                    }
                }
            }
            writer.close();
        } catch (IOException e) {
            Logging.logWarning("Failed to write `VCF` format file to " + outputFile.getAbsolutePath());
        }
    }

}