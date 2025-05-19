package utility;

import java.io.File;
import java.io.IOException;

/**
 * Utility class for performing various validation checks.
 * <p>
 * This class provides static methods to validate strings and file objects.
 * It includes methods to check if a string represents a percentage or a positive double,
 * and methods to validate file and directory objects.
 */
public final class Validation {

    /**
     * Checks if a {@link String} represents a percentage value.
     * <p>
     * This method attempts to parse the input string as a double and checks if the value
     * lies within the range [0.0, 1.0], inclusive. If the parsing fails or the value is
     * out of range, the method returns false.
     *
     * @param s The {@link String} to validate as a percentage.
     * @return {@code true} if the string represents a valid percentage, {@code false} otherwise.
     */
    public static Boolean isPercentage(String s) {
        try {
            double d = Double.parseDouble(s);
            return (d >= 0.0 && d <= 1.0);
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Checks if a {@link String} represents a positive double value.
     * <p>
     * This method attempts to parse the input string as a double and checks if the value
     * is greater than 0.0. If the parsing fails or the value is not positive, the method
     * returns false.
     *
     * @param s The {@link String} to validate as a positive double.
     * @return {@code true} if the string represents a positive double, {@code false} otherwise.
     */
    public static Boolean isPositiveDouble(String s) {
        try {
            return Double.parseDouble(s) > 0.0;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Checks the validity of a {@link File} object.
     * <p>
     * This method validates the given {@link File} object by ensuring it:
     * <ul>
     *   <li>Represents a file (not a directory).</li>
     *   <li>Exists in the file system.</li>
     *   <li>Can be read by the application.</li>
     *   <li>Is not empty.</li>
     * </ul>
     *
     * @param file The {@link File} object to validate.
     * @throws IOException If the file is not valid; i.e., one of the conditions is not met.
     */
    public static void checkFile(File file) throws IOException {
        if (!file.canRead()) {
            throw new IOException("File %s is not readable.".formatted(file.getAbsolutePath()));
        }
        if (!file.isFile()) {
            throw new IOException("File %s is not a file.".formatted(file.getAbsolutePath()));
        }
        if (file.length() == 0) {
            throw new IOException("File %s is empty.".formatted(file.getAbsolutePath()));
        }
    }

    /**
     * Checks if a {@link File} object exists, is a directory, and can be read.
     * <p>
     * This method validates the given {@link File} object by ensuring it:
     * <ul>
     *   <li>Represents a directory (not a file).</li>
     *   <li>Exists in the file system.</li>
     *   <li>Can be read by the application.</li>
     * </ul>
     *
     * @param directory The {@link File} object to validate.
     * @throws IOException If the file is not valid; i.e., one of the conditions is not met.
     */
    public static void checkDirectory(File directory) throws IOException {
        if (!directory.canRead()) {
            throw new IOException("Directory %s is not readable.".formatted(directory.getAbsolutePath()));
        }
        if (!directory.isDirectory()) {
            throw new IOException("Directory %s is not a directory.".formatted(directory.getAbsolutePath()));
        }
    }

}
