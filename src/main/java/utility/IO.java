package utility;

import com.google.common.base.Splitter;
import datastructure.Contig;
import datastructure.Feature;
import datastructure.Storage;
import datastructure.VariantInformation;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFUtils;
import main.Musial;
import org.apache.commons.codec.binary.Base64;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;
import org.apache.commons.lang3.tuple.Triple;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Utility class for input/output operations.
 * <p>
 * This final class provides static methods for various file and data handling operations,
 * such as reading, writing, compressing, and hashing files and strings. It also includes
 * methods for generating specific file formats like VCF, FASTA, and GFF.
 * <p>
 * The class is designed to be non-instantiable and serves as a collection of utility methods.
 */
public final class IO {

    /**
     * Reads the content of a file line by line and returns a list of non-empty, trimmed lines.
     * <p>
     * This method uses a {@link Scanner} to read the file with UTF-8 encoding. Each line is trimmed
     * to remove leading and trailing whitespace, and empty lines are excluded from the result.
     *
     * @param file The {@link File} object representing the file to read.
     * @return An {@link ArrayList} containing the non-empty, trimmed lines from the file.
     * @throws IOException If an I/O error occurs while reading the file.
     */
    public static ArrayList<String> readFile(File file) throws IOException {
        ArrayList<String> lines = new ArrayList<>();
        try (LineIterator lineIterator = FileUtils.lineIterator(file)) {
            while (lineIterator.hasNext()) {
                String line = lineIterator.nextLine().trim();
                if (!line.isEmpty()) {
                    lines.add(line);
                }
            }
        }
        return lines;
    }

    /**
     * Writes the specified content to a file at the given path.
     * <p>
     * This method ensures that the parent directories of the target file are created if they do not exist.
     * It then writes the provided content to the file using UTF-8 encoding. If the file already exists,
     * its content is overwritten.
     *
     * @param path    The {@link Path} where the file will be written.
     * @param content The {@link String} content to write to the file.
     * @throws IOException If an I/O error occurs during directory creation or file writing.
     */
    public static void writeFile(Path path, String content) throws IOException {
        File file = path.toFile();
        FileUtils.createParentDirectories(file); // Ensure parent directories exist.
        FileUtils.write(file, content, StandardCharsets.UTF_8, false); // Write content to the file.
    }

    /**
     * Reads a tabular file and converts its content into a nested map structure.
     * <p>
     * This method reads a tabular file where the first row contains headers and each subsequent row contains data.
     * The first column is treated as the key for the outer map, and the remaining columns are stored in an inner map
     * with their corresponding headers as keys. The file can use tab or comma as delimiters.
     *
     * @param file The {@link File} object representing the tabular file to read.
     * @return A {@link HashMap} where the outer map's key is the first column's value, and the value is another
     * {@link HashMap} containing the remaining columns as key-value pairs.
     * @throws IOException If an I/O error occurs or the file format is invalid.
     */
    public static HashMap<String, HashMap<String, String>> readTabularFileAsNestedMap(File file) throws IOException {
        HashMap<String, HashMap<String, String>> result = new HashMap<>();
        try (BufferedReader br = Files.newBufferedReader(file.toPath(), StandardCharsets.UTF_8)) {
            String[] headers = br.readLine().split("[\t,]", -1);
            if (headers.length < 2) {
                throw new IOException("Invalid file format: must have at least one identifier and one data column.");
            }
            br.lines().forEach(line -> {
                String[] values = line.split("[\t,]", -1);
                if (values.length != headers.length) {
                    throw new RuntimeException("Row length does not match header length.");
                }
                HashMap<String, String> rowMap = new HashMap<>();
                for (int i = 1; i < headers.length; i++) {
                    rowMap.put(headers[i], values[i]);
                }
                result.put(values[0], rowMap);
            });
        }
        return result;
    }

    /**
     * Detects the separator used in a list of strings.
     * <p>
     * This method analyzes the provided list of strings to determine the separator used in the content.
     * It skips lines that start with a specific sign (defined by {@link Constants#SIGN}) and checks the first
     * non-skipped line for the presence of either a tab character or a comma. If a tab is found, it returns
     * the tab separator; if a comma is found, it returns the comma separator. If neither is found, it returns
     * an empty string.
     *
     * @param content A {@link List} of {@link String} objects representing the content to analyze.
     * @return A {@link String} representing the detected separator: either a tab, a comma, or an empty string
     * if no separator is found.
     */
    public static String detectSeparator(List<String> content) {
        return content.stream()
                .filter(line -> !line.startsWith(Constants.SIGN)) // Skip lines starting with the defined sign
                .findFirst() // Find the first non-skipped line
                .map(line -> line.contains(Constants.TAB) ? Constants.TAB // Check for tab separator
                        : line.contains(Constants.COMMA) ? Constants.COMMA // Check for comma separator
                        : Constants.EMPTY) // Return empty string if no separator is found
                .orElse(Constants.EMPTY); // Return empty string if no valid line is found
    }

    /**
     * Generates the content of a plain VCF (Variant Call Format) file.
     * <p>
     * This method constructs a VCF file content as a {@link String} from a list of variants.
     * The VCF content includes the file format, source, and a header line, followed by the variant data.
     * Each variant is represented by its chromosome, position, reference base, and alternate base.
     * <p>
     * The generated VCF content follows the VCFv4.3 specification and includes the following fields:
     * <ul>
     *   <li>CHROM: Chromosome name</li>
     *   <li>POS: Position of the variant</li>
     *   <li>ID: Variant identifier (set to ".")</li>
     *   <li>REF: Reference base(s)</li>
     *   <li>ALT: Alternate base(s) (gaps are stripped)</li>
     *   <li>QUAL: Quality score (set to "100")</li>
     *   <li>FILTER: Filter status (set to ".")</li>
     *   <li>INFO: Additional information (empty)</li>
     * </ul>
     *
     * @param variants A list of {@link Tuple} objects, where each tuple contains:
     *                 <ul>
     *                   <li>A {@link Triple} with the chromosome name, position, and alternate base.</li>
     *                   <li>A {@link VariantInformation} object containing the reference base.</li>
     *                 </ul>
     * @return A {@link String} representing the VCF content.
     */
    public static String generateVcfContent(ArrayList<Tuple<Triple<String, Integer, String>, VariantInformation>> variants) {
        StringBuilder content = new StringBuilder();
        content.append("##fileformat=VCFv4.3").append(Constants.lineSeparator)
                .append("##source=MUSIAL").append(Constants.lineSeparator)
                .append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").append(Constants.lineSeparator);
        for (Tuple<Triple<String, Integer, String>, VariantInformation> variant : variants) {
            content.append(variant.a.getLeft()).append("\t") // CHROM
                    .append(variant.a.getMiddle()).append("\t") // POS
                    .append(".\t") // ID
                    .append(variant.b.getReferenceBaseString(true)).append("\t") // REF
                    .append(SequenceOperations.stripGaps(variant.a.getRight())).append("\t") // ALT
                    .append("100\t") // QUAL
                    .append(".\t") // FILTER
                    .append("\t").append(Constants.lineSeparator); // INFO
        }
        return content.toString();
    }

    /**
     * Generates the content of a reference FASTA file from the given {@link Storage} object.
     * <p>
     * This method constructs a FASTA file content as a {@link String} by iterating over the contigs
     * in the provided {@link Storage} object. Each contig's name is used as the header (prefixed with '>'),
     * and its sequence is split into lines of 80 characters for proper FASTA formatting.
     *
     * @param storage The {@link Storage} object containing the contigs and their sequences.
     * @return A {@link String} representing the content of the reference FASTA file.
     * @throws IOException              If an I/O error occurs during the generation of the FASTA content.
     * @throws IllegalArgumentException If no reference sequence information is stored in the {@link Storage} object.
     */
    public static String generateReferenceFastaContent(Storage storage) throws IOException {
        if (storage.hasMissingContigSequences())
            throw new IllegalArgumentException("No reference sequence information is stored in the specified storage.");
        StringBuilder content = new StringBuilder();
        for (Contig contig : storage.getContigs()) {
            content.append(">").append(contig.name).append(Constants.lineSeparator);
            Splitter.fixedLength(80).split(contig.getSequence()).forEach(line -> content.append(line).append(Constants.lineSeparator));
        }
        return content.toString();
    }

    /**
     * Generates the content of a GFF (General Feature Format) file from the given {@link Storage} object.
     * <p>
     * This method constructs a GFF file content as a {@link String} by iterating over the features
     * in the provided {@link Storage} object. The GFF content includes the version, processor information,
     * and the feature data. Each feature is converted to its GFF string representation using the
     * {@link Feature#toGffString()} method.
     * <p>
     * The generated GFF content follows the GFF3 specification and includes the following:
     * <ul>
     *   <li>##gff-version: Specifies the GFF version.</li>
     *   <li>##processor: Includes the software name and version used to generate the file.</li>
     *   <li>Feature data: Each feature is represented in GFF format.</li>
     * </ul>
     *
     * @param storage The {@link Storage} object containing the features to include in the GFF file.
     * @return A {@link String} representing the GFF file content.
     */
    public static String generateGffContent(Storage storage) {
        StringBuilder content = new StringBuilder();
        content.append("##gff-version 3.1.26").append(Constants.lineSeparator);
        content.append("##processor %s %s".formatted(Musial.softwareName, Musial.softwareVersion)).append(Constants.lineSeparator);
        for (Feature feature : storage.getFeatures()) {
            content.append(feature.toGffString());
        }
        return content.toString();
    }

    /**
     * Initializes a {@link VCFFileReader} instance for the passed VCF file. For this, a temporary indexed VCF file is created.
     *
     * @param file A {@link File} object pointing to a .vcf file.
     * @return A {@link VCFFileReader} instance for the passed .vcf file.
     * @throws IOException In case of an error during the initialization of the VCFFileReader.
     */
    public static VCFFileReader initializeVCFFileReader(File file) throws IOException {
        File vcfFile = VCFUtils.createTemporaryIndexedVcfFromInput(file, String.valueOf(file.hashCode()));
        vcfFile.deleteOnExit();
        return new VCFFileReader(vcfFile);
    }

    /**
     * Copies a resource from the application's classpath to a specified target {@link Path}.
     * <p>
     * This method retrieves a resource as an {@link InputStream} from the application's classpath
     * using the specified resource path. The resource is then copied to the target file path,
     * overwriting any existing file at the target location.
     *
     * @param resourceName The path to the resource within the application's classpath.
     * @param targetPath   The file path where the resource should be copied.
     * @throws MusialException If the resource cannot be found or an I/O error occurs during the copy operation.
     */
    public static void copyResourceToFile(String resourceName, Path targetPath) throws MusialException {
        try (InputStream resourceStream = Musial.class.getResourceAsStream(resourceName)) {
            Files.copy(Objects.requireNonNull(resourceStream), targetPath, StandardCopyOption.REPLACE_EXISTING);
        } catch (IOException e) {
            throw new MusialException(String.format("Failed to extract resource %s. %s", resourceName, e.getMessage()));
        }
    }

    /**
     * Compresses a string using GZIP compression and encodes the result in Base64.
     * <p>
     * This method compresses the input string using the GZIP algorithm and then encodes
     * the compressed byte array into a Base64 string. The method ensures proper resource
     * management by using a try-with-resources block for the output streams.
     *
     * @param content The {@link String} to be compressed.
     * @return A Base64-encoded {@link String} representing the GZIP-compressed content.
     * @throws IOException If an I/O error occurs during compression.
     */
    public static String gzipCompress(String content) throws IOException {
        try (ByteArrayOutputStream outputStream = new ByteArrayOutputStream(1024);
             GZIPOutputStream gzipOutputStream = new GZIPOutputStream(outputStream)) {
            gzipOutputStream.write(content.getBytes());
            gzipOutputStream.finish();
            return Base64.encodeBase64String(outputStream.toByteArray());
        }
    }

    /**
     * Decompresses a Base64-encoded GZIP-compressed string.
     * <p>
     * This method decodes the input string from Base64, decompresses the resulting GZIP-compressed data,
     * and returns the decompressed content as a string. It uses a buffer to read the decompressed data
     * in chunks and appends it to a {@link StringBuilder}.
     *
     * @param content The Base64-encoded GZIP-compressed string to decompress.
     * @return A {@link String} containing the decompressed content.
     * @throws IOException If an I/O error occurs during decompression.
     */
    public static String gzipDecompress(String content) throws IOException {
        try (GZIPInputStream gzipInputStream = new GZIPInputStream(
                new ByteArrayInputStream(Base64.decodeBase64(content)))) {
            byte[] buffer = new byte[1024];
            StringBuilder output = new StringBuilder();
            int length;
            while ((length = gzipInputStream.read(buffer)) != -1) {
                output.append(new String(buffer, 0, length));
            }
            return output.toString();
        }
    }

    /**
     * Generates the MD5 hash of the given string.
     * <p>
     * This method computes the MD5 hash of the input string and returns it as a hexadecimal string.
     * It uses the {@link org.apache.commons.codec.digest.DigestUtils#md5Hex(String)} method from the
     * Apache Commons Codec library to perform the hashing.
     *
     * @param content The {@link String} to hash.
     * @return A {@link String} representing the MD5 hash of the input content in hexadecimal format.
     */
    public static String md5Hash(String content) {
        return org.apache.commons.codec.digest.DigestUtils.md5Hex(content);
    }

}