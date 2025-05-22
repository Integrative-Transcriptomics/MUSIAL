package utility;

import main.Musial;

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

    /**
     * Initializes the logging system with a specified logging level.
     * <p>
     * This method configures the logger to use a custom `ConsoleHandler` with a formatter
     * that styles log messages using ANSI escape codes for different log levels (INFO, CONFIG, WARNING, SEVERE).
     * It also disables the default parent handlers to prevent duplicate logging.
     * <p>
     * The logging level can be adjusted by passing a `Level` parameter, which determines
     * the minimum severity of messages that will be logged.
     *
     * @param level The logging level to set for the logger (e.g., Level.INFO, Level.WARNING).
     */
    public static void init(Level level) {
        // Turn off the default logger handlers.
        logger.setUseParentHandlers(false);

        // Create ConsoleHandler with custom Formatter.
        ConsoleHandler ch = new ConsoleHandler();
        ch.setFormatter(new Formatter() {
            /**
             * Formats a log record with a custom prefix based on its severity level.
             * <p>
             * The prefix is styled using ANSI escape codes for visual emphasis.
             * The formatted message includes the timestamp and the log message.
             *
             * @param record The log record to format.
             * @return The formatted log message as a string.
             */
            public String format(LogRecord record) {
                String prefix;
                switch (record.getLevel().getName()) {
                    case "INFO" -> prefix = "\033[0;34mINFO\033[0m    ";
                    case "CONFIG" -> prefix = "\033[0;94mCONFIG\033[0m  ";
                    case "WARNING" -> prefix = "\033[0;33mWARNING\033[0m ";
                    case "SEVERE" -> prefix = "\033[0;93mSEVERE\033[0m  ";
                    default -> prefix = "LOG     "; // Default
                }
                prefix += "%s ".formatted(getTimestamp());
                return prefix + record.getMessage() + System.lineSeparator();
            }
        });

        // Set the handler's logging level to ALL to capture all messages.
        ch.setLevel(Level.ALL);
        logger.addHandler(ch);

        // Set the logger's logging level to the specified level.
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
     * This method logs messages at the INFO level, which is typically used for general informational
     * messages that highlight the progress of the application at a coarse-grained level.
     * The message is automatically formatted and logged using the application's logger.
     *
     * @param msg The informational message to be logged.
     */
    public static void logInfo(String msg) {
        logger.log(Level.INFO, msg);
    }

    /**
     * Logs a message indicating completion with a timestamp.
     * <p>
     * This method logs messages at the INFO level, indicating that a specific task or operation
     * has been successfully completed. The message is formatted with a "DONE" label styled using
     * ANSI escape codes for visual emphasis.
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
     * This method logs messages at the SEVERE level, which is typically used for critical error messages
     * that indicate a failure in the application. The message is automatically formatted and logged
     * using the application's logger.
     *
     * @param msg The error message to be logged.
     */
    public static void logSevere(String msg) {
        logger.log(Level.SEVERE, msg);
    }

    /**
     * Logs a critical exit message to the console with a timestamp.
     * <p>
     * This method is used to log messages indicating a critical application exit.
     * The message is formatted with a timestamp and an "EXIT" label styled using ANSI escape codes
     * for visual emphasis. The log level is set to SEVERE, which is the highest level of logging severity.
     *
     * @param msg The exit message to be logged.
     */
    public static void logExit(String msg) {
        logger.log(Level.SEVERE, "%s %s".formatted("\033[1;91mEXIT\033[0m", msg));
    }

    /**
     * Logs a warning message to the console with a timestamp.
     * <p>
     * This method logs messages at the WARNING level, which is typically used to indicate
     * potential issues or situations that require attention but are not critical errors.
     * The message is automatically formatted and logged using the application's logger.
     *
     * @param msg The warning message to be logged.
     */
    public static void logWarning(String msg) {
        logger.log(Level.WARNING, msg);
    }

    /**
     * Logs a warning message to the console with a timestamp, but only once for each unique key.
     * <p>
     * This method ensures that a warning message associated with a specific key is logged only once.
     * It uses a set (`logDump`) to track keys of already logged warnings. If the key is not present
     * in the set, the warning message is logged, and the key is added to the set.
     * <p>
     * This is useful for avoiding repetitive logging of the same warning message.
     *
     * @param key The unique key identifying the warning message.
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
     * Returns the current date formatted as a string.
     * <p>
     * This method retrieves the current date, formats it using the pattern "dd-MM-yy",
     * and returns it as a string. The format includes the day, month, and year in a two-digit format.
     *
     * @return The current date formatted as "dd-MM-yy".
     */
    public static String getDate() {
        return new SimpleDateFormat("dd-MM-yy").format(new Date());
    }
}
