package task;

import cli.CLIExpand;
import exceptions.MusialException;
import model.Storage;
import op.*;
import util.Logging;

import java.io.IOException;
import java.nio.file.Path;

public class ExecutorExpand {

    /**
     * Command-line interface for the expand task.
     */
    private final CLIExpand cli;

    /**
     * Storage instance to hold genomic data.
     */
    private final Storage storage;

    /**
     * Updater for the storage instance.
     */
    private final StorageUpdater storageUpdater;

    /**
     * Processor for VCF files.
     */
    private final VCFProcessor vcfProcessor;

    /**
     * Annotator for genetic variants.
     */
    private final VariantAnnotator variantAnnotator;

    public ExecutorExpand(CLIExpand cli) throws IOException {
        this.cli = cli;
        Logging.logInfo("Load storage.");
        this.storage = StorageFactory.fromPath(cli.input);
        this.storageUpdater = new StorageUpdater(storage);
        vcfProcessor = new VCFProcessor(cli.vcfFiles, storage, storage.hasReference());
        variantAnnotator = new VariantAnnotator(storage);
    }

    /**
     * @noinspection DuplicatedCode
     */
    public void run() throws IOException, MusialException {
        // Process VCF files and load variants into storage.
        Logging.logInfo("Load variant calls.");
        vcfProcessor.analyzeFiles();
        Logging.logInfo("Processed %d variant calls from %d VCF file(s). %d calls were ignored, %d calls were filtered.".formatted(
                vcfProcessor.getProcessedCalls(), cli.vcfFiles.size(), vcfProcessor.getIgnoredCalls(), vcfProcessor.getFilteredCalls()));
        Logging.logInfo("Total of %d samples were loaded.".formatted(storage.getSamples().size()));
        storageUpdater.updateSampleAttributes(cli.vcfMeta);
        storageUpdater.updateVariants();
        Logging.logInfo("Total of %d variants were loaded.".formatted(storage.getVariantsCount()));

        // Determine working path for output files.
        Path path;
        if (cli.overwrite) {
            path = cli.input;
        } else {
            path = cli.output;
        }

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
            variantAnnotator.runSnpEff(path.getParent());
        }

        // Infer sequence types if reference sequences are available.
        if (storage.parameters.skipTyping()) {
            Logging.logInfo("Skip sequence typing per user request.");
        } else if (storage.getFeatures().isEmpty()) {
            Logging.logWarning("Skip sequence typing; no features are available.");
        } else {
            Logging.logInfo("Run sequence typing.");
            storageUpdater.updateSequenceTypes();
        }

        // Compute statistics for the storage.
        Logging.logInfo("Recompute statistics.");
        storageUpdater.updateStatistics();

        // Write the storage data to the specified output file.
        Logging.logInfo("Write storage to file: " + path.toAbsolutePath());
        StorageIO.toJSON(storage, path);
    }

}
