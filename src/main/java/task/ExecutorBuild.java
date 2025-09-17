package task;

import cli.CLIBuild;
import exceptions.MusialException;
import model.Storage;
import op.*;
import util.Logging;

import java.io.IOException;

/**
 * The {@code ExecutorBuild} class is responsible for executing the build task for genomic data analysis.
 * <p>
 * This class orchestrates the process of loading genomic features, processing variant call files (VCF), annotating variants, inferring
 * sequence types, and computing statistics.
 */
public class ExecutorBuild {

    /**
     * Command-line interface for the build task.
     */
    private final CLIBuild cli;

    /**
     * Storage instance to hold genomic data.
     */
    private final Storage storage;

    /**
     * Updater for the storage instance.
     */
    private final StorageUpdater storageUpdater;

    /**
     * Loader for genomic features.
     */
    private final FeatureLoader featureLoader;

    /**
     * Annotator for genetic variants.
     */
    private final VariantAnnotator variantAnnotator;

    /**
     * Constructs an instance of the {@link ExecutorBuild} class.
     * <p>
     * This constructor initializes the necessary components for the build task, including:
     * <ul>
     *   <li>Storage for genomic data, created based on the command-line interface (CLI) input.</li>
     *   <li>A {@link StorageUpdater} to manage updates to the storage.</li>
     *   <li>A {@link FeatureLoader} to load and validate genomic features.</li>
     *   <li>A {@link VCFProcessor} to process variant call files (VCF).</li>
     *   <li>A {@link VariantAnnotator} to annotate genetic variants.</li>
     * </ul>
     *
     * @param cli The {@link CLIBuild} instance containing the command-line arguments and options.
     * @throws IOException If an I/O error occurs during the initialization of storage or other components.
     */
    public ExecutorBuild(CLIBuild cli) throws IOException {
        this.cli = cli;
        Logging.logInfo("Initialize storage.");
        storage = StorageFactory.fromCli(this.cli);
        storageUpdater = new StorageUpdater(storage);
        featureLoader = new FeatureLoader(storage, cli.featureList, cli.features);
        variantAnnotator = new VariantAnnotator(storage);
        Logging.logDone("");
    }

    /**
     * Runs the {@code build} task for genomic data analysis.
     * <p>
     * This method performs the following steps:
     * <ul>
     *   <li>Loads and validates genomic features from the input files.</li>
     *   <li>Processes VCF files to load variant calls into the storage.</li>
     *   <li>Runs variant annotation using SnpEff if applicable.</li>
     *   <li>Infers sequence types if reference sequences are available.</li>
     *   <li>Computes statistics for the genomic data in the storage.</li>
     *   <li>Writes the processed storage data to the specified output file.</li>
     * </ul>
     *
     * @throws MusialException If an error occurs during the processing of genomic data.
     * @throws IOException     If an I/O error occurs during file operations.
     * @noinspection DuplicatedCode
     */
    public void run() throws MusialException, IOException {
        // Load and validate genomic features.
        Logging.logInfo("Load and validate genomic features.");
        featureLoader.loadFeatures();
        featureLoader.validateFeatures();
        Logging.logDone("Processed %d of %d annotated features.".formatted(featureLoader.getLoadedFeatureCount(),
                cli.featureList.size()));

        // Process VCF files and load variants into storage.
        try (VCFProcessor vcfProcessor = new VCFProcessor(cli.vcfFiles, storage, !storage.hasReference())) {
            Logging.logInfo("Analyze VCF files.");
            vcfProcessor.processFiles();
            Logging.logDone("Processed %d variant calls from %d VCF file(s). %d calls were ignored, %d calls were filtered.".formatted(
                    vcfProcessor.getProcessedCallsCount(), cli.vcfFiles.size(), vcfProcessor.getIgnoredCallsCount(),
                    vcfProcessor.getFilteredCallsCount()));

            // Update variants from the processed VCF data.
            Logging.logInfo("Update variants.");
            vcfProcessor.updateVariants();
            Logging.logDone("");
        }

        // Update sample attributes from metadata.
        storageUpdater.updateSampleAttributes(cli.vcfMeta);

        // Check and run SnpEff annotation if applicable.
        if (storage.parameters.skipAnnotation()) {
            Logging.logInfo("Skip variant annotation per user request.");
        } else if (!storage.hasReference()) {
            Logging.logWarning("Skip variant annotation; no reference sequence is available.");
        } else if (storage.getFeatures().isEmpty()) {
            Logging.logWarning("Skip variant annotation; no features are available.");
        } else if (storage.getFeatures().stream().allMatch(f -> f.type.equals("region"))) {
            Logging.logWarning("Skip variant annotation; all features are of type region.");
        } else if (!storage.hasNovelVariants()) {
            Logging.logWarning("Skip variant annotation; no novel variants to annotate.");
        } else {
            Logging.logInfo("Run variant annotation with SnpEff.");
            variantAnnotator.runSnpEff(cli.output.getParent());
            Logging.logDone("");
        }

        // Infer sequence types if reference sequences are available.
        if (storage.parameters.skipTyping()) {
            Logging.logInfo("Skip sequence typing per user request.");
        } else if (storage.getFeatures().isEmpty()) {
            Logging.logWarning("Skip sequence typing; no features are available.");
        } else {
            Logging.logInfo("Run sequence typing.");
            storageUpdater.updateSequenceTypes();
            Logging.logDone("");
        }

        // Compute statistics for the storage.
        Logging.logInfo("Compute statistics.");
        storageUpdater.updateStatistics();
        Logging.logDone("");

        // Write the storage data to the specified output file.
        Logging.logInfo("Write storage to file: " + cli.output.toAbsolutePath());
        StorageIO.toJSON(storage, cli.output);
        Logging.logDone("Storage contains %d features, %d samples, and %d variants.".formatted(storage.getFeatures().size(),
                storage.getSamples().size(), storage.getVariantsCount()));
    }

}
