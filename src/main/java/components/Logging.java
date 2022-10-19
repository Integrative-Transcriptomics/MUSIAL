package components;

import cli.CLIColors;
import main.Musial;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

/**
 * This class comprises static methods used for logging.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class Logging {

    public static final ArrayList<String> CACHE = new ArrayList<>();

    /**
     * Prints a simple {@link String} message without any additional information.
     *
     * @param msg {@link String} to be printed.
     */
    public static void logMsg(String msg) {
        if (!Musial.SILENT) {
            System.out.println(msg);
        }
        Logging.cache(msg);
    }

    /**
     * Adds any {@link String} content to the {@link Logging#CACHE}.
     *
     * @param content {@link String} content to be added.
     */
    public static void cache(String content) {
        Logging.CACHE.add(content);
    }

    /**
     * Prints a colorized banner yielding software information.
     */
    public static void logSoftwareInfo() {
        ArrayList<String> textContents = new ArrayList<>() {{
            add(" " + Musial.NAME + " " + Musial.VERSION);
            add("License: " + Musial.LICENSE);
            add("Contact: " + Musial.CONTACT + " ");
        }};
        if (!Musial.SILENT) {
            String textFormat = CLIColors.WHITE_BACKGROUND + CLIColors.BLACK_BOLD;
            System.out.println(textFormat + String.join(" | ", textContents) + CLIColors.RESET);
        }
        Logging.cache(String.join(" | ", textContents));
    }

    /**
     * Prints a colorized status message with time-stamp.
     *
     * @param msg {@link String} to be printed.
     */
    public static void logStatus(String msg) {
        String timeStamp = Logging.getTimestampDayTime();
        if (!Musial.SILENT) {
            String textFormat = CLIColors.BLUE_BACKGROUND_BRIGHT + CLIColors.BLACK_BOLD;
            System.out.println(textFormat + " (" + timeStamp + ") STATUS  " + CLIColors.RESET + " " + msg);
        }
        Logging.cache("(" + timeStamp + ") STATUS >" + msg);
    }

    /**
     * Prints a colorized error message with time-stamp.
     *
     * @param msg {@link String} to be printed.
     */
    public static void logError(String msg) {
        String timeStamp = Logging.getTimestampDayTime();
        if (!Musial.SILENT) {
            String textFormat = CLIColors.RED_BACKGROUND_BRIGHT + CLIColors.BLACK_BOLD;
            System.out.println(textFormat + " (" + timeStamp + ") ERROR   " + CLIColors.RESET + " " + msg);
        }
        Logging.cache("(" + timeStamp + ") ERROR >" + msg);
    }

    /**
     * Prints a colorized warning message with time-stamp.
     *
     * @param msg {@link String} to be printed.
     */
    public static void logWarning(String msg) {
        String timeStamp = Logging.getTimestampDayTime();
        if (!Musial.SILENT) {
            String textFormat = CLIColors.YELLOW_BACKGROUND_BRIGHT + CLIColors.BLACK_BOLD;
            System.out.println(textFormat + " (" + timeStamp + ") WARNING " + CLIColors.RESET + " " + msg);
        }
        Logging.cache("(" + timeStamp + ") WARNING >" + msg);
    }

    /**
     * Returns a colorized {@link String} for pretty console output.
     *
     * @param parameter {@link String}; The content to colorize.
     * @return Colorized input {@link String}.
     */
    public static String colorParameter(String parameter) {
        return CLIColors.YELLOW_UNDERLINED + parameter + CLIColors.RESET;
    }

    /**
     * Generates a colorized {@link String} DONE tag.
     *
     * @return Colorized DONE {@link String}.
     */
    public static String getDoneTag() {
        return CLIColors.GREEN_BACKGROUND_BRIGHT + CLIColors.BLACK_BOLD + "DONE" + CLIColors.RESET;
    }

    /**
     * Generates a colorized {@link String} START tag.
     *
     * @return Colorized START {@link String}.
     */
    public static String getStartTag() {
        return CLIColors.CYAN_BACKGROUND_BRIGHT + CLIColors.BLACK_BOLD + "START" + CLIColors.RESET;
    }

    /**
     * Generates a colorized {@link String} tag with custom text content.
     *
     * @param content {@link String} content to use in tag.
     * @return Colorized custom {@link String}.
     */
    public static String getCustomTag(String content) {
        return CLIColors.PURPLE_BACKGROUND_BRIGHT + CLIColors.BLACK_BOLD + content + CLIColors.RESET;
    }

    /**
     * Generates, formats and returns a {@link String} representing the current date and time.
     *
     * @return A {@link String} representing the current date and time in format `dd-MM-yyyy HH:mm:ss`.
     */
    public static String getTimestampDayTime() {
        return new SimpleDateFormat("dd-MM-yyyy HH:mm:ss").format(new Date());
    }

    /**
     * Generates, formats and returns a {@link String} representing the current time.
     *
     * @return A {@link String} representing the current date and time in format `HH:mm:ss`.
     */
    public static String getTimestampTime() {
        return new SimpleDateFormat("HH:mm:ss").format(new Date());
    }
}
