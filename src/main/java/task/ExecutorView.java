package task;

import cli.CLIView;
import model.Storage;
import op.StorageFactory;
import op.StorageTable;
import util.Logging;

import java.io.IOException;
import java.util.Objects;

/**
 * The {@code ExecutorView} class is responsible for executing the {@code view} task for inspecting data.
 * <p>
 * This class initializes the CLI view, loads the storage, and manages the storage table. It provides functionality to apply queries,
 * populate the storage table based on the content type, and either print or write the output to a file.
 */
public class ExecutorView {

    /**
     * Command-line interface for the view task.
     */
    private final CLIView cli;

    /**
     * Storage instance to get genomic data.
     *
     * @noinspection FieldCanBeLocal
     */
    private final Storage storage;

    /**
     * Storage table for managing and displaying data.
     */
    private final StorageTable storageTable;

    /**
     * Constructs an instance of the {@link ExecutorView} class.
     * <p>
     * This constructor initializes the CLI view, loads the storage from the input path, and creates a storage table for managing data.
     *
     * @param cli The {@link CLIView} object containing user input and options.
     * @throws IOException If an error occurs while loading the storage.
     */
    public ExecutorView(CLIView cli) throws IOException {
        this.cli = cli;
        Logging.logInfo("Load storage.");
        this.storage = StorageFactory.fromPath(cli.input);
        this.storageTable = new StorageTable(storage);
        Logging.logDone("");
    }

    /**
     * Executes the main logic of the {@code view} task.
     * <p>
     * This method applies the query filters if provided, populates the storage table based on the content type, and either prints the table
     * to the console or writes it to the specified output file.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Applies query filters to the storage table if the query is not empty.</li>
     *     <li>Populates the storage table based on the content type (features, samples, or variants).</li>
     *     <li>Prints the table to the console if no output file is specified, or writes it to the output file.</li>
     * </ul>
     *
     * @throws IllegalStateException If the content type is unexpected.
     */
    public void run() {
        Logging.logInfo("View %s.".formatted(cli.content));
        if (!cli.query.isEmpty()) {
            this.storageTable.setFilters(cli.query);
        }
        switch (cli.content) {
            case FEATURES -> storageTable.populateFromFeatures();
            case SAMPLES -> storageTable.populateFromSamples();
            case VARIANTS -> storageTable.populateFromVariants();
            default -> throw new IllegalStateException("Unexpected value: " + cli.content);
        }

        if (Objects.isNull(cli.output)) {
            this.storageTable.print(true);
        } else {
            this.storageTable.write(cli.output.toAbsolutePath().toString(), '\t');
        }
        Logging.logDone("");
    }
}
