package utility;

import exceptions.MusialException;

import java.io.File;

/**
 * Utility class for operating system-related functionalities.
 * <p>
 * This class provides methods for executing command-line commands and handling
 * file and directory operations. It is designed to facilitate interaction with
 * the underlying operating system from within the application.
 * </p>
 */
public class OS {

    /**
     * Executes a command-line command using a {@link ProcessBuilder}.
     * <p>
     * This method runs the specified command in a separate process. It allows redirecting
     * error logs, output, and setting the working directory for the process.
     * </p>
     *
     * @param command  An array of strings representing the command and its arguments.
     * @param errorLog The file path to redirect error logs. If empty, errors are not redirected.
     * @param output   The file path to redirect standard output. If empty, output is not redirected.
     * @param runInDir The directory in which the command should be executed. If empty, the default directory is used.
     */
    public static void runCommand(String[] command, String errorLog, String output, String runInDir) throws MusialException {
        try {
            ProcessBuilder pb = new ProcessBuilder(command);
            if (errorLog.length() > 0) {
                pb = pb.redirectError(new File(errorLog));
            }
            if (output.length() > 0) {
                pb = pb.redirectOutput(new File(output));
            }
            if (runInDir.length() > 0) {
                pb.directory(new File(runInDir));
            }
            Process process = pb.start();
            process.waitFor();
        } catch (Exception e) {
            throw new MusialException("Error executing command %s; %s.".formatted(String.join(" ", command), e.getMessage()));
        }
    }

}
