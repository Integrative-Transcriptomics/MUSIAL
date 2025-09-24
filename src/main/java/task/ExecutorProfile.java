package task;

import cli.CLIProfile;
import model.Storage;
import op.StorageFactory;
import op.StorageTable;
import util.Logging;

import java.io.IOException;
import java.util.Objects;

/**
 * The {@code ExecutorProfile} class is responsible for executing the {@code profile} task to profile samples.
 * <p>
 * This class initializes the CLI profile, loads the storage, and manages the storage table. It provides functionality to apply queries,
 * populate the storage table based on the content type, and either print or write the output to a file.
 */
public class ExecutorProfile {

    /**
     * Command-line interface for the profile task.
     */
    private final CLIProfile cli;

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
     * Constructs an instance of the {@link ExecutorProfile} class.
     * <p>
     * This constructor initializes the CLI profile, loads the storage from the input path, and creates a storage table for managing data.
     *
     * @param cli The {@link CLIProfile} object containing user input and options.
     * @throws IOException If an error occurs while loading the storage.
     */
    public ExecutorProfile(CLIProfile cli) throws IOException {
        this.cli = cli;
        Logging.logInfo("Load storage.");
        this.storage = StorageFactory.fromPath(cli.input);
        this.storageTable = new StorageTable(storage);
        Logging.logDone("");
    }

    /**
     * Executes the main logic of the {@code ExecutorProfile} class.
     * <p>
     * This method applies query filters if provided, populates the storage table based on the content type, and either prints the table to
     * the console or writes it to the specified output file.
     * <p>
     * The method performs the following steps:
     * <ul>
     *     <li>Checks if the query is not empty and applies the filters to the storage table.</li>
     *     <li>Profiles the storage table based on the content type (VARIANTS, ALLELES, or PROTEOFORMS).</li>
     *     <li>Prints the table to the console if no output file is specified, or writes it to the output file.</li>
     * </ul>
     * <p>
     * If the content type is unexpected, an {@link IllegalStateException} is thrown.
     */
    public void run() {
        Logging.logInfo("Profile %s.".formatted(cli.content));
        if (!cli.query.isEmpty()) {
            this.storageTable.setFilters(cli.query);
        }
        switch (cli.content) {
            case VARIANTS -> storageTable.populateWithVariantsProfile(cli.reduced);
            case ALLELES -> storageTable.populateWithAlleleProfile(cli.reduced);
            case PROTEOFORMS -> storageTable.populateWithProteoformProfile(cli.reduced);
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
