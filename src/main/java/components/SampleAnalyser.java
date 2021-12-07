package components;

import com.google.common.collect.Sets;
import datastructure.FeatureAnalysisEntry;
import datastructure.SampleAnalysisEntry;
import datastructure.VariantPositionsTable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CopyOnWriteArraySet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import main.Musial;
import me.tongfei.progressbar.ProgressBar;
import runnables.SampleAnalyserRunnable;
import utility.Logging;

/**
 * Comprises static methods to analyze samples.
 * <p>
 * In MUSIAL this is done by multi-threading. In detail the cartesian product of a set of
 * {@link FeatureAnalysisEntry} and a set of {@link SampleAnalysisEntry} instances is generated. For more detail
 * view the respective class documentations.
 * <p>
 * This results in a set of pairs of {@link FeatureAnalysisEntry} and {@link SampleAnalysisEntry} objects, each
 * specifying a sample that is analyzed against a reference genome region (which may be a single gene or the whole
 * genome). For each such pair a {@link SampleAnalyserRunnable} is initialized and run.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class SampleAnalyser {

  /**
   * Runs multi-threaded sample analysis.
   * <p>
   * Constructs each possible pair of the passed referenceEntries and sampleEntries elements. For each pair a new
   * single-threaded {@link SampleAnalyserRunnable} is initialized.
   *
   * @param referenceEntries Set of {@link FeatureAnalysisEntry} instances. Each specifies a locus on the reference
   *                         genome.
   * @param sampleEntries    Set of {@link SampleAnalysisEntry} instances. Each specifies the `.vcf` file content from
   *                         one sample.
   * @param arguments        {@link ArgumentsParser} containing arguments parsed from command line.
   * @return {@link VariantPositionsTable} containing the information returned from each single
   * {@link SampleAnalyserRunnable}.
   * @throws InterruptedException If the thread was interrupted.
   */
  public static VariantPositionsTable run(HashSet<FeatureAnalysisEntry> referenceEntries,
                                          CopyOnWriteArraySet<SampleAnalysisEntry> sampleEntries,
                                          ArgumentsParser arguments)
      throws InterruptedException {
    ProgressBar progress = Musial.buildProgress();
    try (progress) {
      ExecutorService executor = Executors.newFixedThreadPool(arguments.getNumThreads());
      Set<List<Object>> runEntries = Sets.cartesianProduct(referenceEntries, sampleEntries);
      VariantPositionsTable variantPositionsTable = new VariantPositionsTable(sampleEntries.size());
      progress.maxHint(runEntries.size() + 1);
      for (List<Object> runEntry : runEntries) {
        FeatureAnalysisEntry featureAnalysisEntry = (FeatureAnalysisEntry) runEntry.get(0);
        SampleAnalysisEntry sampleAnalysisEntry = (SampleAnalysisEntry) runEntry.get(1);
        variantPositionsTable.addSampleToReference(featureAnalysisEntry.name, sampleAnalysisEntry.sampleName);
        executor.execute(
            new SampleAnalyserRunnable(
                sampleAnalysisEntry.sampleName,
                featureAnalysisEntry,
                sampleAnalysisEntry.vcfFileReader
                    .query(featureAnalysisEntry.referenceSequenceLocation,
                        featureAnalysisEntry.locationStart,
                        featureAnalysisEntry.locationEnd),
                variantPositionsTable,
                arguments.getMinCoverage(),
                arguments.getMinFrequency(),
                arguments.getMinHet(),
                arguments.getMaxHet(),
                arguments.getMinQuality(),
                progress
            )
        );
      }
      executor.shutdown();
      //noinspection ResultOfMethodCallIgnored
      executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
      addReferenceToVariantPositionsTable(variantPositionsTable, referenceEntries, progress);
      progress.setExtraMessage(Logging.getDoneMessage());
      return variantPositionsTable;
    }
  }

  /**
   * Adds information from the reference genome specified by the single elements of referenceEntries to the passed
   * variantPositionsTable.
   *
   * @param variantPositionsTable {@link VariantPositionsTable} containing information about the variant sites of a
   *                              set of samples against possibly multiple reference loci.
   * @param referenceEntries      Set of {@link FeatureAnalysisEntry} instances, each specifying a reference locus.
   * @param progress              {@link ProgressBar} instance to visualize progress information to the user.
   */
  private static void addReferenceToVariantPositionsTable(
      VariantPositionsTable variantPositionsTable, HashSet<FeatureAnalysisEntry> referenceEntries,
      ProgressBar progress) {
    for (FeatureAnalysisEntry referenceEntry : referenceEntries) {
      char[] referenceSequenceBases = referenceEntry.getReferenceSequence().toCharArray();
      int referenceSequenceStart = referenceEntry.locationStart;
      for (int i = 0; i < referenceSequenceBases.length; i++) {
        int position = referenceSequenceStart + i;
        variantPositionsTable
            .putVariablePosition(referenceEntry.name, "Reference", String.valueOf(position),
                referenceSequenceBases[i], Double.NaN, Double.NaN, Double.NaN, true,
                new HashMap<>());
      }
      progress.step();
    }
  }

}
