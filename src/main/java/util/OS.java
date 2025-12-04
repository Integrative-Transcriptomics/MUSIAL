package util;

import java.io.File;

/**
 * Utility class for operating system-related functionalities.
 * <p>
 * This class provides methods for executing command-line commands and handling file and directory operations. It is designed to facilitate
 * interaction with the underlying operating system from within the application.
 */
public final class OS {

    /**
     * Private constructor to prevent instantiation of this utility class.
     */
    private OS() {
    }

    /**
     * Executes a command-line command using a {@link ProcessBuilder}.
     * <p>
     * This method runs the specified command in a separate process. It allows redirecting error logs, output, and setting the working
     * directory for the process.
     *
     * @param command  An array of strings representing the command and its arguments.
     * @param errorLog The file path to redirect error logs. If empty, errors are not redirected.
     * @param output   The file path to redirect standard output. If empty, output is not redirected.
     * @param runInDir The directory in which the command should be executed. If empty, the default directory is used.
     * @throws RuntimeException If an error occurs while executing the command.
     */
    public static void runCommand(String[] command, String errorLog, String output, String runInDir) throws RuntimeException {
        try {
            ProcessBuilder pb = new ProcessBuilder(command);
            if (!errorLog.isEmpty()) {
                pb = pb.redirectError(new File(errorLog));
            }
            if (!output.isEmpty()) {
                pb = pb.redirectOutput(new File(output));
            }
            if (!runInDir.isEmpty()) {
                pb.directory(new File(runInDir));
            }
            Process process = pb.start();
            process.waitFor();
        } catch (Exception e) {
            throw new RuntimeException("Error executing command %s; %s.".formatted(String.join(" ", command), e.getMessage()));
        }
    }

}
