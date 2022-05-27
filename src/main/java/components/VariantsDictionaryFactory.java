package components;

import cli.CLIParametersUpdateVDict;
import datastructure.FeatureEntry;
import datastructure.VariantsDictionary;
import exceptions.MusialIntegrityException;
import java.io.IOException;
import runnables.SampleAnalyzerRunnable;
import utility.IO;
import utility.Logging;

/**
 * Comprises static methods to analyze samples in order to build a new variant database.
 * <p>
 * For a specified {@link FeatureEntry} a pool of {@link SampleAnalyzerRunnable} is created
 * building a {@link VariantsDictionary} instance concurrently by parsing in the specified `vcf` files.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class VariantsDictionaryFactory {

  public static VariantsDictionary build(CLIParametersUpdateVDict arguments)
      throws InterruptedException, IOException, MusialIntegrityException {
    if (arguments.outputFile.exists()) {
      Logging.logStatus("Read-in existing variants dictionary from " + arguments.outputFile);
      return IO.readVariantsDictionary(arguments.outputFile);
    } else {
      Logging.logStatus("Build new variants dictionary at " + arguments.outputFile);
      String chromosome = null;
      for (FeatureEntry featureEntry : arguments.features.values()) {
        if (chromosome == null) {
          chromosome = featureEntry.chromosome;
        } else if (!featureEntry.chromosome.equals(chromosome)) {
          throw new MusialIntegrityException(
              "All features specified for one variant dictionary have to be located on the same chromosome; features " +
                  "were specified for " + chromosome + " and " + featureEntry.chromosome + ".");
        }
      }
      return new VariantsDictionary(arguments.minCoverage, arguments.minHomFrequency, arguments.minHetFrequency,
          arguments.maxHetFrequency, arguments.minQuality, chromosome);
    }
    /*
    pb.maxHint(sampleEntries.size());
    pb.stepTo(0);
    ExecutorService executor = Executors.newFixedThreadPool(arguments.threads);
    for (SampleEntry sampleEntry : sampleEntries) {
      executor.execute(
          new SampleAnalyzerRunnable(
              sampleEntry,
              sampleEntry.vcfFileReader.query(locus, 1, referenceEntry.getSequence().length()),
              VariantsDictionary,
              pb
          )
      );
    }
    executor.shutdown();
    //noinspection ResultOfMethodCallIgnored
    executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    pb.setExtraMessage(Logging.getDoneMessage());
    return VariantsDictionary;
     */
  }

}
