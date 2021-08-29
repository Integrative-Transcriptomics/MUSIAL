package utility;

import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * This class comprises static methods used for status logging.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class Logging {

  /**
   * @param message A {@link String} printed as info message.
   */
  public static void logStatus(String message) {
    System.out.println( "\n[" + Logging.getTimestampDayTime( ) + "] STATUS: " + message  );
  }

  public static void logError(String message) {
    System.out.println( "\n[" + Logging.getTimestampDayTime( ) + "] ERROR: " + message  );
    System.exit(-1);
  }

  public static String getDoneMessage() {
    return "Done [" + Logging.getTimestampTime() + "]";
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
