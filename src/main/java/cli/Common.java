package cli;

import exceptions.MusialException;
import htsjdk.samtools.util.FileExtensions;
import main.Musial;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.file.PathUtils;
import util.IO;
import util.Logging;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Stream;

/**
 * Common utility class for parsing configuration parameters and file paths.
 */
abstract class Common {

    /**
     * Parses the VCF (Variant Call Format) files from the configuration map.
     * <p>
     * This method retrieves the "vcfFiles" parameter from the configuration map and processes its entries. Each entry is checked to
     * determine if it is a directory or a file. If it is a directory, all VCF files within the directory are added to the list. If it is a
     * file, it is added directly if it has a valid VCF extension. The method ensures that at least one valid VCF file is found; otherwise,
     * an exception is thrown.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Checks if the "vcfFiles" parameter exists in the configuration map.</li>
     *     <li>Iterates over the entries in the "vcfInput" list.</li>
     *     <li>For directories, collects all valid VCF files.</li>
     *     <li>For files, validates the extension and adds them to the list.</li>
     *     <li>Throws an exception if no valid VCF files are found.</li>
     * </ul>
     * <p>
     * If the "vcfFiles" parameter is missing or no valid files are found, appropriate exceptions are thrown.
     * <p>
     * This method is set static to be reusable in the {@link CLIExpand} class.
     *
     * @param configuration A {@link Map} containing the configuration parameters. The "vcfFiles" parameter specifies the list of
     *                      directories or files to process.
     * @return A {@link List} of {@link Path} objects representing the valid VCF files.
     * @throws MusialException If the "vcfFiles" parameter is missing or no valid VCF files are found.
     * @throws IOException     If an I/O error occurs while accessing the files or directories.
     */
    static List<Path> parseInputVcfFiles(Map<String, Object> configuration) throws MusialException, IOException {
        if (!configuration.containsKey("vcfFiles")) {
            throw new MusialException("`vcfFiles` must be specified in the configuration.");
        }

        List<Path> vcfFiles = new ArrayList<>();
        //noinspection unchecked
        for (String entry : ((List<String>) configuration.get("vcfFiles"))) {
            Function<Path, Boolean> isVcf = p -> (p.toString().endsWith(FileExtensions.VCF) ||
                    p.toString().endsWith(FileExtensions.COMPRESSED_VCF));

            Path path = Path.of(entry);
            if (Files.isDirectory(path)) {
                try (Stream<Path> stream = Files.list(path)) {
                    List<Path> paths = stream.filter(isVcf::apply).toList();
                    vcfFiles.addAll(paths);
                }
            } else if (isVcf.apply(path)) {
                vcfFiles.add(path);
            }
        }
        if (vcfFiles.isEmpty()) {
            throw new MusialException("No valid VCF files found in `vcfFiles`.");
        }
        return vcfFiles;
    }

    /**
     * Parses the VCF metadata from the configuration map.
     * <p>
     * This method retrieves the "vcfMeta" parameter from the configuration map and validates the file path. If the file exists, is not a
     * directory, and is not empty, it reads the file as a nested map using the {@link IO#readTabularFileAsNestedMap(File)} utility method.
     * If the "vcfMeta" parameter is not provided or the file is invalid, an empty map is returned.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Checks if the "vcfMeta" parameter exists in the configuration map.</li>
     *     <li>Validates the file path to ensure it is a regular, non-empty file.</li>
     *     <li>Reads and parses the file into a nested map.</li>
     * </ul>
     * <p>
     * If the file is invalid, an empty map is returned.
     * <p>
     * This method is set static to be reusable in the {@link CLIExpand} class.
     *
     * @param configuration A {@link Map} containing the configuration parameters. The "vcfMeta" parameter specifies the path to the
     *                      metadata file.
     * @return A {@link Map} where the key is the sample ID, and the value is another map containing metadata attributes and their values.
     * @throws IOException If an I/O error occurs while accessing the file.
     */
    static Map<String, Map<String, String>> parseInputVcfMeta(Map<String, Object> configuration) throws IOException {
        Map<String, Map<String, String>> vcfMeta = new HashMap<>();
        // Check if an annotation file is specified in the configuration; if not, return an empty map.
        if (configuration.containsKey("vcfMeta")) {
            Path path = Path.of((String) configuration.get("vcfMeta"));
            if (PathUtils.isRegularFile(path) && !PathUtils.isDirectory(path) && !PathUtils.isEmptyFile(path)) {
                vcfMeta = IO.readTabularFileAsNestedMap(path.toFile());
            }
        }
        return vcfMeta;
    }

    /**
     * Parses the input storage file path from the command-line arguments.
     * <p>
     * This method retrieves the value of the "I" option from the provided {@link CommandLine} arguments and validates the file path. The
     * file must exist, not be a directory, and have a valid extension (either `.json` or `.json.gz`). If the file does not meet these
     * criteria, a {@link MusialException} is thrown.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Retrieves the file path from the "I" option.</li>
     *     <li>Validates that the file exists, is not a directory, and has a valid extension.</li>
     *     <li>Logs the validated file path.</li>
     * </ul>
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @return A {@link Path} object representing the validated input storage file path.
     * @throws MusialException If the file does not exist, is a directory, or has an invalid extension.
     */
    static Path parseInputStorageFile(CommandLine arguments) throws MusialException {
        Path path = Path.of(arguments.getOptionValue("I"));
        if (!Files.exists(path) || Files.isDirectory(path) || !(path.toString().endsWith(".json") || path.toString().endsWith(".json.gz"))) {
            throw new MusialException("Input storage file must be a valid .json or .json.gz file.");
        }
        Logging.logConfig("`input` set to %s.".formatted(path));
        return path;
    }

    /**
     * Parses the output path from the configuration map.
     * <p>
     * This method retrieves the "output" parameter from the configuration map and validates the file path. If the path is a directory, it
     * appends a default file name with a timestamp and extension. If the path is invalid or cannot be created, a {@link MusialException} is
     * thrown.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Checks if the "output" parameter exists in the configuration map.</li>
     *     <li>Validates the file path and creates parent directories if necessary.</li>
     *     <li>If the path is a directory, appends a default file name with a timestamp and extension.</li>
     * </ul>
     * <p>
     * If the path is invalid or cannot be created, an exception is thrown.
     *
     * @param configuration A {@link Map} containing the configuration parameters. The "output" parameter specifies the path to the output
     *                      directory or file.
     * @return A {@link Path} object representing the validated output path.
     * @throws MusialException If the "output" parameter is missing or the path is invalid.
     */
    static Path parseOutputStorageFile(Map<String, Object> configuration) throws MusialException {
        if (!configuration.containsKey("output")) {
            throw new MusialException("`output` must be specified in the configuration.");
        }

        Path path = Path.of((String) configuration.get("output"));
        try {
            // Note: The behavior of this method is slightly different than expected, as the toFile() transformation seems to remove
            //  trailing file separator symbols from the path.
            FileUtils.createParentDirectories(path.toFile());
            if (Files.isDirectory(path))
                path = path.resolve("musial-storage-%s.%s".formatted(Logging.getDate(), Musial.OUTPUT_EXTENSION));
            Logging.logConfig("`output` set to %s.".formatted(path));
        } catch (Exception e) {
            throw new MusialException("Failed to validate path %s specified for `output`.".formatted(path));
        }

        return path;
    }

    /**
     * Parses the output file path from the command-line arguments.
     * <p>
     * This method determines the output file path based on the "o" option in the provided {@link CommandLine} arguments. If the "o" option
     * is not specified, it generates a default file path using the "I" option's parent directory, appending a file name with the provided
     * suffix and extension. If the "o" option is specified as "stdout" or "print", it returns null. Otherwise, it validates the specified
     * path, creates parent directories if necessary, and appends a default file name if the path is a directory.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Checks if the "o" option is provided in the arguments.</li>
     *     <li>If not provided, generates a default path based on the "I" option's parent directory.</li>
     *     <li>If "o" is "stdout" or "print", logs the output and returns null.</li>
     *     <li>Otherwise, validates the specified path, creates parent directories, and appends a default file name if needed.</li>
     * </ul>
     * <p>
     * If the path is invalid or cannot be created, a {@link MusialException} is thrown.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @param suffix    A {@link String} representing the suffix to append to the file name.
     * @param extension A {@link String} representing the file extension to use.
     * @return A {@link Path} object representing the validated output file path, or null if the output is set to "stdout" or "print".
     * @throws MusialException If the path is invalid or cannot be created.
     */
    static Path parseOutputFile(CommandLine arguments, String suffix, String extension) throws MusialException {
        if (!arguments.hasOption("o")) {
            // Generate a default path based on the "I" option's parent directory
            Path path = Path.of(arguments.getOptionValue("I")).getParent();
            path = path.resolve("musial-%s-%s.%s".formatted(suffix, Logging.getDate(), extension));
            Logging.logConfig("`output` set to %s.".formatted(path));
            return path;
        } else {
            String specified = arguments.getOptionValue("o");
            if (specified.equals("stdout") || specified.equals("print")) {
                // Log and return null for "stdout" or "print"
                Logging.logConfig("`output` set to %s.".formatted(specified));
                return null;
            } else {
                // Validate the specified path and create parent directories if necessary
                Path path = Path.of(arguments.getOptionValue("o"));
                try {
                    FileUtils.createParentDirectories(path.toFile());
                    if (Files.isDirectory(path))
                        path = path.resolve("musial-%s-%s.%s".formatted(suffix, Logging.getDate(), extension));
                    Logging.logConfig("`output` set to %s.".formatted(path));
                    return path;
                } catch (Exception e) {
                    throw new MusialException("Failed to validate path %s specified for `output`.".formatted(path));
                }
            }
        }
    }

    /**
     * Parses the query parameters from the command-line arguments.
     * <p>
     * This method retrieves the "q" option from the provided {@link CommandLine} arguments and processes its values. If the "q" option is
     * not specified, it returns an empty set. Otherwise, it retrieves the values associated with the "q" option, converts them into a
     * {@link Set}, and returns the result.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Checks if the "q" option is provided in the arguments.</li>
     *     <li>If not provided, returns an empty set.</li>
     *     <li>If provided, retrieves the values, converts them into a set, and returns the set.</li>
     * </ul>
     * <p>
     * Note: This method does not validate the query parameters; it simply returns them as provided.
     *
     * @param arguments The {@link CommandLine} object containing the parsed command-line arguments.
     * @return A {@link Set} of {@link String} objects representing the query parameters, or an empty set if the "q" option is not
     * specified.
     */
    static Set<String> parseQuery(CommandLine arguments) {
        if (!arguments.hasOption("q")) {
            return Collections.emptySet();
        } else {
            String[] entries = arguments.getOptionValues("q");
            return new HashSet<>(Arrays.asList(entries));
        }
    }

}
