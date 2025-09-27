package util;

import exceptions.MusialException;
import main.Musial;
import org.apache.commons.codec.binary.Base64;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.RandomStringUtils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Utility class for I/O related operations.
 * <p>
 * This class provides static methods for various file and data handling operations, such as reading, writing, compressing, and hashing
 * files and strings.
 */
public final class IO {

    /**
     * Private constructor to prevent instantiation of this utility class.
     */
    private IO() {
    }

    /**
     * A DecimalFormat instance for formatting frequencies in scientific notation.
     * <p>
     * The format uses one digit before the decimal point and two digits after, followed by an exponent (e.g., "1.23E4"). The locale is set
     * to US for consistent decimal and grouping symbols.
     */
    private static final DecimalFormat frequencyFormat = new DecimalFormat(".00E0", DecimalFormatSymbols.getInstance(Locale.US));

    /**
     * A DecimalFormat instance for formatting numbers with up to three decimal places.
     * <p>
     * The format uses up to two digits before the decimal point and three digits after (e.g., "12.345"). The locale is set to US for
     * consistent decimal and grouping symbols.
     */
    private static final DecimalFormat decimalFormat = new DecimalFormat("##.###", DecimalFormatSymbols.getInstance(Locale.US));

    /**
     * A RandomStringUtils instance for generating random strings.
     * <p>
     * This instance is used for generating random strings, typically for creating temporary file or directory names.
     */
    private static final RandomStringUtils randomStringUtils = RandomStringUtils.insecure();

    /**
     * Formats a frequency value into scientific notation.
     * <p>
     * This method formats the given frequency value using the {@link #frequencyFormat}. If the formatted value equals ".10E1", it is
     * replaced with "1.00E0" for consistency.
     *
     * @param value The frequency value to format.
     * @return A {@link String} representing the formatted frequency in scientific notation.
     */
    public static String formatFrequency(double value) {
        String formattedValue = frequencyFormat.format(value);
        if (formattedValue.equals(".10E1"))
            return "1.00";
        else if (formattedValue.equals(".00E0"))
            return "0.00";
        else
            return formattedValue;
    }

    /**
     * Formats a number to three decimal places.
     * <p>
     * This method formats the given number using the {@link #decimalFormat}. The result is a string representation of the number with up to
     * three decimal places.
     *
     * @param value The number to format.
     * @return A {@link String} representing the formatted number.
     */
    public static String formatNumber(double value) {
        return decimalFormat.format(value);
    }

    /**
     * Generates a random alphanumeric string of the specified length.
     * <p>
     * This method uses the {@link RandomStringUtils} instance to generate a random string consisting of both letters and digits. The
     * generated string is suitable for use in scenarios where a random identifier or token is needed.
     *
     * @param length The length of the random alphanumeric string to generate.
     * @return A {@link String} containing the generated random alphanumeric characters.
     */
    public static String randomAlphanumeric(int length) {
        return randomStringUtils.nextAlphanumeric(length);
    }

    /**
     * Writes the specified content to a file at the given path.
     * <p>
     * This method ensures that the parent directories of the target file are created if they do not exist. It then writes the provided
     * content to the file using UTF-8 encoding. If the file already exists, its content is overwritten.
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
     * This method reads a tabular file where the first row contains headers and each subsequent row contains data. The first column is
     * treated as the key for the outer map, and the remaining columns are stored in an inner map with their corresponding headers as keys.
     * The file can use tab or comma as delimiters.
     *
     * @param file The {@link File} object representing the tabular file to read.
     * @return A {@link HashMap} where the outer map's key is the first column's value, and the value is another {@link HashMap} containing
     * the remaining columns as key-value pairs.
     * @throws IOException If an I/O error occurs or the file format is invalid.
     */
    public static Map<String, Map<String, String>> readTabularFileAsNestedMap(File file) throws IOException {
        Map<String, Map<String, String>> result = new HashMap<>();
        try (BufferedReader br = Files.newBufferedReader(file.toPath(), StandardCharsets.UTF_8)) {
            String[] headers = br.readLine().split("[\t,]", -1);
            if (headers.length < 2) {
                throw new IOException("Invalid format of file %s, must have at least two columns, separated by `comma` or `tab`.");
            }
            br.lines().forEach(line -> {
                String[] values = line.split("[\t,]", -1);
                if (values.length != headers.length) {
                    throw new RuntimeException("Row length does not match header length.");
                }
                Map<String, String> rowMap = new HashMap<>();
                for (int i = 1; i < headers.length; i++) {
                    rowMap.put(headers[i], values[i]);
                }
                result.put(values[0], rowMap);
            });
        }
        return result;
    }

    /**
     * Copies a resource from the application's classpath to a specified target {@link Path}.
     * <p>
     * This method retrieves a resource as an {@link InputStream} from the application's classpath using the specified resource path. The
     * resource is then copied to the target file path, overwriting any existing file at the target location.
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
     * This method compresses the input string using the GZIP algorithm and then encodes the compressed byte array into a Base64 string. The
     * method ensures proper resource management by using a try-with-resources block for the output streams.
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
     * This method decodes the input string from Base64, decompresses the resulting GZIP-compressed data, and returns the decompressed
     * content as a string. It uses a buffer to read the decompressed data in chunks and appends it to a {@link StringBuilder}.
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
     * This method computes the MD5 hash of the input string and returns it as a hexadecimal string. It uses the
     * {@link org.apache.commons.codec.digest.DigestUtils#md5Hex(String)} method from the Apache Commons Codec library to perform the
     * hashing.
     *
     * @param content The {@link String} to hash.
     * @return A {@link String} representing the MD5 hash of the input content in hexadecimal format.
     */
    public static String md5Hash(String content) {
        return org.apache.commons.codec.digest.DigestUtils.md5Hex(content);
    }

}