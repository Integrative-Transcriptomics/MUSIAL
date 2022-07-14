package components;

import java.io.File;

/**
 * This class comprises static methods used for validating input data.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public class Validation {

    /**
     * Check if a {@link String} represents a positive integer.
     *
     * @param s The {@link String} to check.
     * @return Whether the passed {@link String} s passes the check.
     */
    public static Boolean isPositiveInteger(String s) {
        try {
            return Integer.parseInt(s) > 0;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Check if a {@link String} represents a percentage.
     *
     * @param s The {@link String} to check.
     * @return Whether the passed {@link String} s passes the check.
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
     * Check if a {@link String} represents a positive double.
     *
     * @param s The {@link String} to check.
     * @return Whether the passed {@link String} s passes the check.
     */
    public static Boolean isPositiveDouble(String s) {
        try {
            return Double.parseDouble(s) > 0.0;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Checks if a {@link File} object exists, is a file and can be read.
     *
     * @param file The {@link File} object to validate.
     * @return Whether the passed {@link File} file passes the check.
     */
    public static boolean isFile(File file) {
        return (file.isFile() && file.exists() && file.canRead());
    }

    /**
     * Checks if a {@link File} object exists, is a directory and can be read.
     *
     * @param file The {@link File} object to validate.
     * @return Whether the passed {@link File} file passes the check.
     */
    public static boolean isDirectory(File directory) {
        return (directory.isDirectory() && directory.exists() && directory.canRead());
    }

}
