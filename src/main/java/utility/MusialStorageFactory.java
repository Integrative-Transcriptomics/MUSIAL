package utility;

import cli.BuildConfiguration;
import datastructure.Feature;
import datastructure.MusialStorage;
import exceptions.MusialException;
import runnables.SampleAnalyzer;

import java.io.File;
import java.io.IOException;

/**
 * Comprises static methods to analyze samples in order to build a MUSIAL storage.
 * <p>
 * For a specified {@link Feature} a pool of {@link SampleAnalyzer} is created
 * building a {@link MusialStorage} instance concurrently by parsing in the specified `vcf` files.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.0
 */
public final class MusialStorageFactory {

    /**
     * Constructs and returns a new {@link MusialStorage} instance from the information specified by the {@link BuildConfiguration} instance.
     *
     * @param configuration {@link BuildConfiguration} instance yielding parameter specification for the MUSIAL BUILD module.
     * @return {@link MusialStorage} instance.
     * @throws MusialException .
     */
    public static MusialStorage build(BuildConfiguration configuration)
            throws MusialException {
        if (configuration.output.exists()) {
            throw new MusialException("(I/O) Specified file '" + configuration.output.getAbsolutePath() + "' to output already exists.");
        } else {
            return new MusialStorage(configuration.referenceSequence, configuration.minimalCoverage, configuration.minimalHomozygousFrequency, configuration.minimalHeterozygousFrequency,
                    configuration.maximalHeterozygousFrequency, configuration.minimalQuality);
        }
    }

    /**
     * Parses an existing {@link MusialStorage} instance from a .json(.br) file.
     *
     * @param inputFile {@link File} instance pointing to the file that contains the {@link MusialStorage} content.
     * @return {@link MusialStorage} instance.
     * @throws MusialException .
     */
    public static MusialStorage load(File inputFile)
            throws MusialException, IOException {
        if (!inputFile.exists()) {
            throw new MusialException("(I/O) Specified file '" + inputFile.getAbsolutePath() + "' to load MUSIAL dump from does not exist.");
        } else {
            Logger.logStatus("Load MUSIAL dump from '" + inputFile.getAbsolutePath() + "'");
            return IO.loadMusialDump(inputFile);
        }
    }

}
