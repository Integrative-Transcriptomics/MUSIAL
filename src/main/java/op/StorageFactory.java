package op;

import cli.CLIBuild;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import model.Storage;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Objects;
import java.util.zip.GZIPInputStream;

/**
 * Factory class for creating and loading {@link Storage} objects.
 */
public class StorageFactory {

    /**
     * Private constructor to prevent instantiation of this utility class.
     */
    private StorageFactory() {
    }

    /**
     * A pre-configured {@link Gson} instance for JSON serialization and deserialization.
     */
    private static final Gson gson = new GsonBuilder()
            .registerTypeAdapter(Storage.class, Storage.typeAdapter()) // Register custom TypeAdapter for Storage
            .create(); // Build the Gson instance

    /**
     * Reads a {@link Storage} object from the specified file path.
     * <p>
     * This method validates the file, determines if it is GZIP-compressed based on its extension, and reads its contents. The contents are
     * then deserialized into a {@link Storage} object using Gson.
     *
     * @param path The {@link Path} to the file containing the serialized {@link Storage} object.
     * @return The deserialized {@link Storage} object.
     * @throws IOException If the file cannot be read or deserialized.
     */
    public static Storage fromPath(Path path) throws IOException {
        // Convert the path to a File object
        File file = path.toFile();

        // Validate the file to ensure it meets the required conditions.
        if (!file.canRead()) {
            throw new IOException("File %s is not readable.".formatted(file.getAbsolutePath()));
        }
        if (!file.isFile()) {
            throw new IOException("File %s is not a file.".formatted(file.getAbsolutePath()));
        }
        if (file.length() == 0) {
            throw new IOException("File %s is empty.".formatted(file.getAbsolutePath()));
        }

        // Read the file and deserialize its contents into a Storage object
        try (BufferedReader bufferedReader = new BufferedReader(
                new InputStreamReader(file.getAbsolutePath().endsWith(".gz")
                        ? new GZIPInputStream(Files.newInputStream(file.toPath())) // Handle GZIP-compressed files
                        : Files.newInputStream(file.toPath())))) { // Handle regular files
            return gson.fromJson(bufferedReader, Storage.class);
        } catch (IOException e) {
            // Throw a new IOException with a detailed error message if reading fails
            throw new IOException("Failed to read storage from file %s; %s"
                    .formatted(file.getAbsolutePath(), e.getMessage()));
        }
    }

    /**
     * Constructs a {@link Storage} object from definitions in {@link cli.CLIBuild}.
     * <p>
     * This method initializes a {@link Storage} object using parameters provided via the {@link cli.CLIBuild} instance. It sets up the
     * storage parameters, optionally populates contigs from a reference, and loads features and sample information from the CLI
     * parameters.
     *
     * @param cli The {@link cli.CLIBuild} instance containing the parameters and definitions for constructing the {@link Storage} object.
     * @return A {@link Storage} object representing the loaded data.
     * @throws IOException If an error occurs while reading files or parsing data.
     */
    public static Storage fromCli(CLIBuild cli) throws IOException {

        // Initialize storage parameters using the CLI-provided values.
        Storage.Parameters parameters = new Storage.Parameters(cli.minimalCoverage, cli.minimalFrequency, cli.maskFiltered,
                cli.skipAnnotation, cli.skipTyping, cli.maskedPositions);
        Storage storage = new Storage(parameters);

        // Populate contigs from the reference file, if provided.
        if (Objects.nonNull(cli.reference)) {
            storage.setReference(cli.reference); // Set the reference in the storage.
        }

        // Return the constructed Storage object.
        return storage;
    }

}