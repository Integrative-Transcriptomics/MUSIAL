package utility;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.typeadapters.RuntimeTypeAdapterFactory;
import datastructure.Feature;
import datastructure.FeatureCoding;
import datastructure.MusialStorage;
import exceptions.MusialException;
import org.apache.commons.io.FileUtils;

import java.io.*;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.function.Function;
import java.util.zip.GZIPInputStream;

/**
 * This class comprises static methods used for reading and writing files.
 *
 * @author Simon Hackl
 * @version 2.3
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
    public static MusialStorage loadMusialDump(File dumpFile) throws IOException, MusialException {
        Function<BufferedReader, MusialStorage> loadDumpFromBufferedReader = reader -> {
            RuntimeTypeAdapterFactory<Feature> adapterFactory =
                    RuntimeTypeAdapterFactory
                            .of(Feature.class, "type", true)
                            .registerSubtype(Feature.class, "non_coding")
                            .registerSubtype(FeatureCoding.class, "coding");
            Gson gson = new GsonBuilder().setPrettyPrinting().registerTypeAdapterFactory(adapterFactory).create();
            return gson.fromJson(reader, MusialStorage.class);
        };

        if (dumpFile.getAbsolutePath().endsWith(".gz")) {
            // Case: Dump file is compressed.
            try (
                    BufferedReader bufferedReader = new BufferedReader(
                            new InputStreamReader(
                                    new GZIPInputStream(
                                            Files.newInputStream(dumpFile.toPath())
                                    )
                            )
                    )
            ) {
                return loadDumpFromBufferedReader.apply(bufferedReader);
            }
        } else if (dumpFile.getAbsolutePath().endsWith(".json")) {
            // Case: Dump file is not compressed.
            try (
                    BufferedReader bufferedReader = new BufferedReader(
                            new InputStreamReader(
                                    Files.newInputStream(dumpFile.toPath())
                            )
                    )
            ) {
                return loadDumpFromBufferedReader.apply(bufferedReader);
            }
        } else {
            throw new MusialException("Unable to load MUSIAL storage from file " + dumpFile.getAbsolutePath() + " (unknown file extension).");
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

}