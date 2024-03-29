package datastructure;

import cli.BuildConfiguration;
import exceptions.MusialException;
import utility.IO;
import utility.Logger;

import java.io.File;
import java.io.IOException;

/**
 * Factory class for {@link MusialStorage} instances. Implements methods to build new or load existing instances from local files.
 *
 * @author Simon Hackl
 * @version 2.3
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
            throw new MusialException("(I/O) Specified output file '" + configuration.output.getAbsolutePath() + "' already exists.");
        } else {
            return new MusialStorage(configuration.referenceSequence, configuration.minimalCoverage, configuration.minimalFrequency, configuration.features.values(), configuration.samples.values(), configuration.excludedPositions);
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
