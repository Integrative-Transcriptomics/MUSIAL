package utility;

import com.aayushatharva.brotli4j.Brotli4jLoader;
import com.aayushatharva.brotli4j.encoder.Encoder;
import com.google.common.base.Splitter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.typeadapters.RuntimeTypeAdapterFactory;
import datastructure.Feature;
import datastructure.FeatureCoding;
import datastructure.MusialStorage;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import main.Musial;
import org.apache.commons.io.FileUtils;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;
import java.util.function.Function;

/**
 * This class comprises static methods used for reading and writing files.
 *
 * @author Simon Hackl
 * @version 2.2
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
     * Initializes a {@link MusialStorage} object from a variants dictionary JSON file.
     *
     * @param dumpFile {@link File} pointing to the JSON file to parse.
     * @return {@link MusialStorage} instance.
     * @throws IOException If any error occurs while parsing the variants dictionary JSON file.
     */
    public static MusialStorage loadMusialDump(File dumpFile) throws IOException {
        Function<BufferedReader, MusialStorage> loadDumpFromBufferedReader = reader -> {
            RuntimeTypeAdapterFactory<Feature> adapterFactory =
                    RuntimeTypeAdapterFactory
                            .of(Feature.class, "type", true)
                            .registerSubtype(Feature.class, "non_coding")
                            .registerSubtype(FeatureCoding.class, "coding");
            Gson gson = new GsonBuilder().setPrettyPrinting().registerTypeAdapterFactory(adapterFactory).create();
            return gson.fromJson(reader, MusialStorage.class);
        };
        if (dumpFile.getAbsolutePath().endsWith(".br")) {
            String decodedDumpFileContent = Compression.brotliDecodeBytes(FileUtils.readFileToByteArray(dumpFile));
            // Case: Dump file is compressed.
            try (
                    BufferedReader bufferedReader = new BufferedReader(new StringReader(decodedDumpFileContent))
            ) {
                return loadDumpFromBufferedReader.apply(bufferedReader);
            }
        } else {
            // Case: Dump file is not compressed.
            try (
                    BufferedReader bufferedReader = new BufferedReader(
                            new InputStreamReader(Files.newInputStream(dumpFile.toPath()),
                                    StandardCharsets.UTF_8))
            ) {
                return loadDumpFromBufferedReader.apply(bufferedReader);
            }
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
     * @throws MusialException If the directory could not be generated.
     */
    public static void generateDirectory(File file) throws MusialException {
        if (!file.isDirectory() && !file.mkdirs()) {
            throw new MusialException("(I/O) Failed to generate directory:\t" + file.getAbsolutePath());
        }
    }

    /**
     * Tries to generate the file specified by the passed {@link File} object.
     *
     * @param file {@link File} object representing a file.
     * @throws MusialException If the file could not be generated.
     */
    public static void generateFile(File file) throws MusialException {
        try {
            if (!file.createNewFile()) {
                throw new MusialException("(I/O) Failed to generate file:\t" + file.getAbsolutePath());
            }
        } catch (IOException e) {
            throw new MusialException(e.getMessage());
        }

    }

    /**
     * Copies the file pointed to with file to the file pointed to with target.
     *
     * @param file   {@link File} object, the source file.
     * @param target {@link File} object, the target file.
     * @throws MusialException If the copy procedure fails.
     */
    @SuppressWarnings("unused")
    public static void copyFile(File file, File target) throws MusialException {
        try {
            Files.copy(file.getAbsoluteFile().toPath(), target.getAbsoluteFile().toPath());
        } catch (IOException e) {
            throw new MusialException(
                    "(I/O) Failed to copy file " + file.getAbsolutePath() + " to target " + target.getAbsolutePath());
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
        List<Chain> chains = structure.getChains();
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
     * // TODO: Fix comment.
     * Writes a fasta format file to the specified output file from a {@link HashMap} instance mapping sequences to
     * lists of identifiers (used to construct the fasta entry headers).
     * <p>
     * Each key of the passed map will be used to build one fasta entry. The header of the respective entry is
     * constructed by joining all {@link String}s of the value accessible via the (sequence) key with the `|` delimiter.
     *
     * @param outputFile {@link File} object pointing to the output fasta file.
     * @param sequences  {@link HashMap} mapping sequences to {@link ArrayList} of strings.
     */
    public static void writeFasta(File outputFile, ArrayList<Tuple<String, String>> sequences) {
        try {
            FileWriter writer = new FileWriter(outputFile);
            for (Tuple<String, String> sequence : sequences) {
                writer.write(sequence.a + "\n");
                for (String l : Splitter.fixedLength(80).split(sequence.b)) {
                    writer.write(l + IO.LINE_SEPARATOR);
                }
                writer.flush();
            }
            writer.close();
        } catch (IOException e) {
            Logger.logWarning("Failed to write `FASTA` format file to " + outputFile.getAbsolutePath());
        }
    }

    /**
     * Function to write a generic file with any line content.
     *
     * @param outputFile  {@link File} object pointing to the output file.
     * @param fileContent {@link String} content to write to file.
     */
    public static void writeFile(File outputFile, String fileContent) throws MusialException {
        if (outputFile.exists()) {
            Logger.logWarning("Failed to write content to file " + outputFile.getAbsolutePath());
        } else {
            generateDirectory(outputFile.getParentFile());
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                writer.write(fileContent);
            } catch (IOException e) {
                Logger.logWarning("Failed to write content to file " + outputFile.getAbsolutePath());
            }
        }
    }

    /**
     * Compress a target file using brotli.
     *
     * @param target {@link File} target to compress.
     */
    public static void compressFile(File target) {
        try {
            String content = Files.readString(target.toPath());
            // Write compressed (brotli) output.
            Brotli4jLoader.ensureAvailability();
            byte[] compressed = Encoder.compress(content.getBytes());
            Files.write(Paths.get(target.getAbsolutePath() + ".br"), compressed);
            //noinspection ResultOfMethodCallIgnored
            target.delete();
        } catch (IOException e) {
            Logger.logWarning("Failed to compress file " + target.getAbsolutePath());
        }
    }

}