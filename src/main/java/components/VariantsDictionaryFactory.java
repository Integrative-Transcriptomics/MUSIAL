package components;

import cli.CLIParametersUpdateVDict;
import datastructure.FeatureEntry;
import datastructure.VariantsDictionary;
import exceptions.MusialIntegrityException;

import java.io.IOException;

import runnables.SampleAnalyzerRunnable;

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
     * Constructs and returns a new {@link VariantsDictionary} instance from the information specified in the passed {@link CLIParametersUpdateVDict} instance.
     *
     * @param cliarguments {@link CLIParametersUpdateVDict} instance yielding parameter specification for the MUSIAL update variants dictionary module.
     * @return {@link VariantsDictionary} instance.
     * @throws InterruptedException     Thrown if any {@link Runnable} implementing class instance is interrupted by the user.
     * @throws IOException              Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialIntegrityException Thrown if any method fails wrt. internal logic, i.e. assignment of proteins to genomes.
     */
    public static VariantsDictionary build(CLIParametersUpdateVDict cliarguments)
            throws InterruptedException, IOException, MusialIntegrityException {
        if (cliarguments.outputFile.exists()) {
            Logging.logStatus("Update existing variants dictionary at " + cliarguments.outputFile + ")");
            return IO.readVariantsDictionary(cliarguments.outputFile);
        } else {
            Logging.logStatus("Generate new variants dictionary at " + cliarguments.outputFile + ")");
            String chromosome = null;
            for (FeatureEntry featureEntry : cliarguments.features.values()) {
                if (chromosome == null) {
                    chromosome = featureEntry.chromosome;
                } else if (!featureEntry.chromosome.equals(chromosome)) {
                    throw new MusialIntegrityException(
                            "All features specified for one variant dictionary have to be located on the same chromosome; features " +
                                    "were specified for " + chromosome + " and " + featureEntry.chromosome + ".");
                }
            }
            return new VariantsDictionary(cliarguments.minCoverage, cliarguments.minHomFrequency, cliarguments.minHetFrequency,
                    cliarguments.maxHetFrequency, cliarguments.minQuality, chromosome);
        }
    }

}
