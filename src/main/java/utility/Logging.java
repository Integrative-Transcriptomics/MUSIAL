package utility;

import main.Musial;
import org.tribuo.clustering.hdbscan.HdbscanTrainer;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.*;

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
     * Logger instance for the application.
     * <p>
     * This logger is used to log messages for the application. It is configured to use the name of the
     * {@link Musial} class as its identifier, which helps in categorizing and filtering log messages.
     */
    public static Logger logger = Logger.getLogger(Musial.class.getName());

    public static void init(Level level) {
        // Turn off the default logger handlers.
        logger.setUseParentHandlers(false);

        // Create ConsoleHandler with custom Formatter.
        ConsoleHandler ch = new ConsoleHandler();
        ch.setFormatter(new Formatter() {
            public String format(LogRecord record) {
                String prefix;
                switch (record.getLevel().getName()) {
                    case "INFO" -> prefix = "\033[0;34mINFO\033[0m    ";
                    case "CONFIG" -> prefix = "\033[0;94mCONFIG\033[0m  ";
                    case "SEVERE" -> prefix = "\033[1;91mSEVERE\033[0m  ";
                    case "WARNING" -> prefix = "\033[1;33mWARNING\033[0m ";
                    default -> prefix = "\033LOG[0m     "; // Default
                }
                prefix += "%s ".formatted(getTimestamp());
                return prefix + record.getMessage() + System.lineSeparator();
            }
        });

        // you can also adjust the level if you like
        ch.setLevel(Level.ALL);
        logger.addHandler(ch);
        logger.setLevel(level);
    }

    /**
     * Prints the software information to the console.
     * <p>
     * This method displays the software name, version, and license in a formatted manner.
     * The output is styled using ANSI escape codes for background and text color.
     * <p>
     * Note: The software information is retrieved from the {@link Musial} class.
     */
    public static void printSoftwareInfo() {
        System.out.printf("\033[47m\033[1;30m| %s %s |\033[0m%n", Musial.softwareName, Musial.softwareVersion);
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
        logger.log(Level.INFO, msg);
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
        logger.log(Level.INFO, "%s %s".formatted("\033[1;92mDONE\033[0m", msg));
    }

    /**
     * Logs a configuration message to the console with a timestamp.
     * <p>
     * This method logs messages at the CONFIG level, which is typically used for static configuration
     * information or messages that help in understanding the application's setup. The message is
     * automatically formatted with a timestamp and logged using the application's logger.
     *
     * @param msg The configuration message to be logged.
     */
    public static void logConfig(String msg) {
        logger.log(Level.CONFIG, msg);
    }

    /**
     * Logs an error message to the console with a timestamp.
     * <p>
     * This method formats the message with a timestamp and an "ERROR" label styled using ANSI escape codes.
     * The formatted message is then printed to the console.
     *
     * @param msg The error message to be logged.
     */
    public static void logSevere(String msg) {
        logger.log(Level.SEVERE, msg);
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
        logger.log(Level.WARNING, msg);
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
