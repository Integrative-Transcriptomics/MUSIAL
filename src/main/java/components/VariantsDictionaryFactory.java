package components;

import cli.CLIColors;
import cli.ModuleParametersBuild;
import datastructure.FeatureEntry;
import datastructure.VariantsDictionary;
import exceptions.MusialException;
import runnables.SampleAnalyzerRunnable;

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
     * Constructs and returns a new {@link VariantsDictionary} instance from the information specified by the {@link ModuleParametersBuild} instance.
     *
     * @param cliarguments {@link ModuleParametersBuild} instance yielding parameter specification for the MUSIAL update variants dictionary module.
     * @return {@link VariantsDictionary} instance.
     * @throws IOException     Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialException Thrown if any method fails wrt. internal logic, i.e. assignment of proteins to genomes.
     */
    public static VariantsDictionary build(ModuleParametersBuild cliarguments)
            throws IOException, MusialException {
        if (cliarguments.outputFile.exists()) {
            Logging.logStatus("Update existing variants dictionary at " + CLIColors.YELLOW_UNDERLINED + cliarguments.outputFile + CLIColors.RESET);
            return IO.readVariantsDictionary(cliarguments.outputFile);
        } else {
            Logging.logStatus("Generate new variants dictionary at " + CLIColors.YELLOW_UNDERLINED + cliarguments.outputFile + CLIColors.RESET);
            return new VariantsDictionary(cliarguments.minCoverage, cliarguments.minHomFrequency, cliarguments.minHetFrequency,
                    cliarguments.maxHetFrequency, cliarguments.minQuality);
        }
    }

}
