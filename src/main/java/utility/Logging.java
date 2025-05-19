package utility;

import main.Musial;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;

/**
 * Logging utility class for printing messages to the console with different severity levels.
 * <p>
 * This class provides methods to log messages with different severity levels (INFO, DONE, ERROR, WARNING) along with
 * timestamps. It also includes a method to print software information and a set to manage warning keys to prevent
 * excessive logging.
 */
public final class Logging {

    /**
     * Set to keep track of logged warnings to avoid duplicate messages.
     */
    public static final Set<String> logDump = new HashSet<>();


    /**
     * Prints the software information to the console.
     * <p>
     * This method displays the software name, version, and license in a formatted manner.
     * The output is styled using ANSI escape codes for background and text color.
     * <p>
     * Note: The software information is retrieved from the {@link Musial} class.
     */
    public static void printSoftwareInfo() {
        System.out.printf("\033[47m\033[1;30m| %s %s | %s |\033[0m%n", Musial.softwareName, Musial.softwareVersion, Musial.softwareLicense);
    }

    /**
     * Logs an informational message to the console with a timestamp.
     * <p>
     * This method formats the message with a timestamp and a "STATUS" label styled using ANSI escape codes.
     * The formatted message is then printed to the console.
     *
     * @param msg The informational message to be logged.
     */
    public static void logInfo(String msg) {
        String timeStamp = Logging.getTimestamp();
        String message = " [" + timeStamp + "] \033[1;36mSTATUS\033[0m  " + msg;
        System.out.println(message);
    }

    /**
     * Logs a message indicating completion with a timestamp.
     * <p>
     * This method formats the message with a timestamp and a "DONE" label styled using ANSI escape codes.
     * The formatted message is then printed to the console.
     *
     * @param msg The message indicating completion to be logged.
     */
    public static void logDone(String msg) {
        String timeStamp = Logging.getTimestamp();
        String message = " [" + timeStamp + "] \033[1;92mDONE\033[0m    " + msg;
        System.out.println(message);
    }

    /**
     * Logs an error message to the console with a timestamp.
     * <p>
     * This method formats the message with a timestamp and an "ERROR" label styled using ANSI escape codes.
     * The formatted message is then printed to the console.
     *
     * @param msg The error message to be logged.
     */
    public static void logError(String msg) {
        String timeStamp = Logging.getTimestamp();
        System.out.println(" [" + timeStamp + "] " + "\033[1;31mERROR\033[0m   " + msg);
    }

    /**
     * Logs a warning message to the console with a timestamp.
     * <p>
     * This method formats the message with a timestamp and a "WARNING" label styled using ANSI escape codes.
     * The formatted message is then printed to the console.
     *
     * @param msg The warning message to be logged.
     */
    public static void logWarning(String msg) {
        String timeStamp = Logging.getTimestamp();
        System.out.println(" [" + timeStamp + "] " + "\033[1;33mWARNING\033[0m " + msg);
    }

    /**
     * Logs a warning message to the console with a timestamp, but only once for each unique key.
     * <p>
     * This method checks if the warning message with the specified key has already been logged.
     * If not, it logs the message and adds the key to the set of logged warnings.
     *
     * @param key The unique key for the warning message.
     * @param msg The warning message to be logged.
     */
    public static void logWarningOnce(String key, String msg) {
        if (logDump.add(key)) {
            Logging.logWarning(msg);
        }
    }

    /**
     * Returns the current timestamp formatted as a string.
     * <p>
     * This method retrieves the current date and time, formats it using the specified date format,
     * and returns it as a string.
     *
     * @return The current timestamp formatted as a string.
     */
    public static String getTimestamp() {
        return new SimpleDateFormat("dd-MM-yy HH:mm").format(new Date());
    }

    /**
     * Returns the current timestamp formatted as a string.
     * <p>
     * This method retrieves the current date, formats it using the specified date format,
     * and returns it as a string.
     *
     * @return The current timestamp formatted as a string.
     */
    public static String getDate() {
        return new SimpleDateFormat("dd-MM-yy").format(new Date());
    }
}
