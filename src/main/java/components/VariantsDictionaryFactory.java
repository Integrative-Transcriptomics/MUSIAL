package components;

import cli.ModuleBuildParameters;
import datastructure.FeatureEntry;
import datastructure.VariantsDictionary;
import exceptions.MusialException;
import runnables.SampleAnalyzerRunnable;

import java.io.File;
import java.io.IOException;

/**
 * Comprises static methods to analyze samples in order to build a new variants dictionary.
 * <p>
 * For a specified {@link FeatureEntry} a pool of {@link SampleAnalyzerRunnable} is created
 * building a {@link VariantsDictionary} instance concurrently by parsing in the specified `vcf` files.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class VariantsDictionaryFactory {

    /**
     * Constructs and returns a new {@link VariantsDictionary} instance from the information specified by the {@link ModuleBuildParameters} instance.
     *
     * @param configuration {@link ModuleBuildParameters} instance yielding parameter specification for the MUSIAL BUILD module.
     * @return {@link VariantsDictionary} instance.
     * @throws IOException     Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialException Thrown if any method fails wrt. internal logic, i.e. assignment of proteins to genomes.
     */
    public static VariantsDictionary build(ModuleBuildParameters configuration)
            throws IOException, MusialException {
        if (configuration.outputFile.exists()) {
            throw new MusialException("(I/O) Specified file (" + configuration.outputFile.getAbsolutePath() + ") to build variants dictionary already exists.");
        } else {
            Logging.logStatus("Generate new variants dictionary at " + Logging.colorParameter(configuration.outputFile.getAbsolutePath()));
            return new VariantsDictionary(configuration.minCoverage, configuration.minHomFrequency, configuration.minHetFrequency,
                    configuration.maxHetFrequency, configuration.minQuality);
        }
    }

    /**
     * Parses an existing {@link VariantsDictionary} file and returns it.
     *
     * @param inputFile {@link File} instance pointing to the file that contains the variants dictionary content.
     * @return {@link VariantsDictionary} instance.
     * @throws MusialException Thrown if any method fails wrt. internal logic, i.e. assignment of proteins to genomes.
     */
    public static VariantsDictionary load(File inputFile)
            throws MusialException, IOException {
        if (!inputFile.exists()) {
            throw new MusialException("(I/O) Specified file (" + inputFile.getAbsolutePath() + ") to load variants dictionary does not exist.");
        } else {
            Logging.logStatus("Load variants dictionary from " + Logging.colorParameter(inputFile.getAbsolutePath()));
            return IO.readVariantsDictionary(inputFile);
        }
    }

}
