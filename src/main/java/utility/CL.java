package utility;

import java.io.File;
import java.io.IOException;

/**
 * This class comprises static methods used for calling command line commands.
 *
 * @author Alexander Seitz
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public class CL {

  /**
   * Runs a command line command from within the application. Currently any errors are ignored.
   *
   * @param command {@link String[]} specifying the command, i.e. white spaces of the command are not considered.
   * @param errorLog {@link String} specifying a path to which errors should be written.
   * @param output {@link String} specifying a path to which the results of the command should be written.
   * @param runInDir {@link String} specifying a path to a directory from which the command should be run.
   */
  public static void runCommand(String[] command, String errorLog, String output, String runInDir){
    try {
      ProcessBuilder pb = new ProcessBuilder(command);
      if(errorLog.length()>0){
        pb = pb.redirectError(new File(errorLog));
      }
      if(output.length()>0){
        pb = pb.redirectOutput(new File(output));
      }
      if(runInDir.length()>0) {
        pb.directory(new File(runInDir));
      }
      Process process = pb.start();
      process.waitFor();
    } catch (Exception e) {
      Logging.logError( e.getMessage() );
    }
  }


}
