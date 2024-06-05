package utility;

import main.Musial;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

/**
 * Comprises static methods used for CL logging.
 *
 * @author Simon Hackl
 */
public final class Logger {

    /**
     * Date formatter.
     */
    private static final DateFormat DATE_FORMAT = new SimpleDateFormat("dd-MM-yy HH:mm:ss");

    /**
     * Print CL text yielding software information.
     */
    public static void printSoftwareInfo() {
        ArrayList<String> textContents = new ArrayList<>() {{
            add("[" + Musial.NAME + " " + Musial.VERSION);
            add("License: " + Musial.LICENSE);
            add("Contact: " + Musial.CONTACT + "]");
        }};
        System.out.println("\033[47m\033[1;30m" + String.join(" | ", textContents) + "\033[0m");
    }

    /**
     * Print status message with time-stamp.
     *
     * @param msg {@link String} to be printed.
     */
    public static void logStatus(String msg) {
        String timeStamp = Logger.getTimestampDayTime();
        String message = " [" + timeStamp + "] STATUS  " + msg;
        System.out.println(message);
    }

    /**
     * Print error message with time-stamp.
     *
     * @param msg {@link String} to be printed.
     */
    public static void logError(String msg) {
        String timeStamp = Logger.getTimestampDayTime();
        System.out.println(" [" + timeStamp + "] " + "\033[1;31mERROR\033[0m" + "   " + msg);
    }

    /**
     * Print warning message with time-stamp.
     *
     * @param msg {@link String} to be printed.
     */
    public static void logWarning(String msg) {
        String timeStamp = Logger.getTimestampDayTime();
        System.out.println(" [" + timeStamp + "] " + "\033[1;33mWARNING\033[0m" + " " + msg);
    }

    /**
     * Generates, formats and returns a {@link String} representing the current date and time.
     *
     * @return A {@link String} representing the current date and time in format `dd-MM-yyyy HH:mm:ss`.
     */
    public static String getTimestampDayTime() {
        return DATE_FORMAT.format(new Date());
    }
}
