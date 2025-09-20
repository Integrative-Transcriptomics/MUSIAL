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
     * Annotator for genetic variants.
     */
    private final VariantAnnotator variantAnnotator;

    /**
     * Initial count of variants in the storage before expansion.
     */
    private final long initialVariantCount;

    /**
     * Initial count of samples in the storage before expansion.
     */
    private final int initialSampleCount;

    public ExecutorExpand(CLIExpand cli) throws IOException {
        this.cli = cli;
        Logging.logInfo("Load storage.");
        this.storage = StorageFactory.fromPath(cli.input);
        this.initialVariantCount = storage.getVariantsCount();
        this.initialSampleCount = storage.getSamples().size();
        this.storageUpdater = new StorageUpdater(storage);
        variantAnnotator = new VariantAnnotator(storage);
        Logging.logDone("");
    }

    /**
     * Runs the {@code expand} task to expand the genomic data storage with new variants and samples.
     * <p>
     * This method performs the following steps:
     * <ul>
     *   <li>Processes VCF files to load variants into the storage.</li>
     *   <li>Reloads existing variant calls from the storage.</li>
     *   <li>Detaches samples already present in the storage to avoid duplication.</li>
     *   <li>Updates variants and sample attributes based on the processed VCF data.</li>
     *   <li>Determines the output path for writing the updated storage data.</li>
     *   <li>Runs variant annotation using SnpEff if applicable.</li>
     *   <li>Performs sequence typing if reference sequences are available.</li>
     *   <li>Recomputes statistics for the storage.</li>
     *   <li>Writes the updated storage data to the specified output file or logs the changes in dry-run mode.</li>
     * </ul>
     *
     * @throws IOException     If an I/O error occurs during file processing or storage operations.
     * @throws MusialException If an error specific to the application logic occurs.
     */
    public void run() throws IOException, MusialException {
        // Process VCF files and load variants into storage.
        if (!cli.vcfFiles.isEmpty()) {
            try (VCFProcessor vcfProcessor = new VCFProcessor(cli.vcfFiles, storage, !storage.hasReference())) {
                Logging.logInfo("Analyze VCF files.");
                vcfProcessor.processFiles();
                Logging.logDone("Processed %d variant calls from %d VCF file(s). %d calls were ignored, %d calls were filtered.".formatted(
                        vcfProcessor.getProcessedCallsCount(), cli.vcfFiles.size(), vcfProcessor.getIgnoredCallsCount(),
                        vcfProcessor.getFilteredCallsCount()));

                // Reload existing variant calls from the storage.
                Logging.logInfo("Load existing variant calls.");
                int loadedCount = vcfProcessor.loadVariantCallsFromStorage();
                Logging.logDone("Loaded %d existing variant calls.".formatted(loadedCount));

                // Detach samples that are already present in the storage to avoid duplication.
                vcfProcessor.getSamples().forEach(sampleIdentifier -> {
                    if (storage.hasSample(sampleIdentifier)) {
                        // Retain sample attributes from existing storage.
                        cli.vcfMeta.put(sampleIdentifier, storage.getSample(sampleIdentifier).getAttributes());

                        // Detach existing sample to avoid duplication.
                        storage.detachSample(sampleIdentifier);
                    }
                });

                // Update variants from the processed VCF data.
                Logging.logInfo("Update variants.");
                vcfProcessor.updateVariants();
                Logging.logDone("");
            }
        } else {
            Logging.logInfo("No VCF files provided, skip VCF file analysis.");
        }

        // Update sample attributes from metadata.
        storageUpdater.updateSampleAttributes(cli.vcfMeta);

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
        } else if (storage.getActiveVariantsCount() == 0) {
            Logging.logWarning("Skip variant annotation; no active variants to annotate.");
        } else {
            Logging.logInfo("Run variant annotation with SnpEff.");
            variantAnnotator.runSnpEff(path.getParent());
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
        Logging.logInfo("Recompute statistics.");
        storageUpdater.updateStatistics();
        Logging.logDone("");

        // Write the storage data to the specified output file.
        if (cli.dry) {
            Logging.logDone("(Dry) Storage expandable by %d samples and %d variants.".formatted(storage.getSamples().size() - initialSampleCount,
                    storage.getVariantsCount() - initialVariantCount));
        } else {
            Logging.logInfo("Write storage to file: " + path.toAbsolutePath());
            StorageIO.toJSON(storage, path);
            Logging.logDone("Storage expanded by %d samples and %d variants.".formatted(storage.getSamples().size() - initialSampleCount,
                    storage.getVariantsCount() - initialVariantCount));
        }
    }

}
