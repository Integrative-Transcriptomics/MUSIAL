/**
 * Encapsulates management of command line argument parsing and validation for each {@link main.MusialTask}.
 * <p>
 * Each implemented class should be paired with one {@link main.MusialTask} and {@link task} class, and implement the {@link cli.CLI}
 * interface. Implemented classes are intended to be initialized in the {@link main.Musial#main} method.
 * <hr>
 * For example the {@link cli.CLIBuild} class is responsible for parsing and validating the {@code build} task command line arguments and is
 * paired with the {@link main.MusialTask#BUILD} task and {@link task.ExecutorBuild} class.
 */
package cli;