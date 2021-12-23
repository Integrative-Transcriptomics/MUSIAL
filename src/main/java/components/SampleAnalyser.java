package components;

import com.google.common.collect.Sets;
import datastructure.ReferenceFeatureEntry;
import datastructure.SampleAnalysisEntry;
import datastructure.VariantContentTable;
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
 * {@link ReferenceFeatureEntry} and a set of {@link SampleAnalysisEntry} instances is generated. For more detail
 * view the respective class documentations.
 * <p>
 * This results in a set of pairs of {@link ReferenceFeatureEntry} and {@link SampleAnalysisEntry} objects, each
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
   * Constructs each possible pair of the passed referenceFeatureEntries and sampleEntries elements. For each pair a new
   * single-threaded {@link SampleAnalyserRunnable} is initialized.
   *
   * @param referenceFeatureEntries Set of {@link ReferenceFeatureEntry} instances. Each specifies a locus on the reference
   *                                genome.
   * @param sampleEntries           Set of {@link SampleAnalysisEntry} instances. Each specifies the `.vcf` file content from
   *                                one sample.
   * @param arguments               {@link ArgumentsParser} containing arguments parsed from command line.
   * @return {@link VariantContentTable} containing the information returned from each single
   * {@link SampleAnalyserRunnable}.
   * @throws InterruptedException If the thread was interrupted.
   */
  public static VariantContentTable run(HashSet<ReferenceFeatureEntry> referenceFeatureEntries,
                                        CopyOnWriteArraySet<SampleAnalysisEntry> sampleEntries,
                                        ArgumentsParser arguments)
      throws InterruptedException {
    VariantContentTable variantContentTable = new VariantContentTable(sampleEntries.size());
    Set<List<Object>> runEntries = Sets.cartesianProduct(referenceFeatureEntries, sampleEntries);
    ProgressBar progress1 = Musial.buildProgress();
    progress1.maxHint(referenceFeatureEntries.size());
    progress1.stepTo(0);
    Logging.logStatus("Storing reference sequence information (" + referenceFeatureEntries.size() + ").");
    try (progress1) {
      addReferenceToVariantContentTable(variantContentTable, referenceFeatureEntries, progress1);
      progress1.setExtraMessage(Logging.getDoneMessage());
    }
    ProgressBar progress2 = Musial.buildProgress();
    progress2.maxHint(runEntries.size());
    progress2.stepTo(0);
    Logging.logStatus("Analysing sample variant call files per reference feature (" + runEntries.size() + ").");
    try (progress2) {
      ExecutorService executor = Executors.newFixedThreadPool(arguments.getNumThreads());
      for (List<Object> runEntry : runEntries) {
        ReferenceFeatureEntry referenceFeatureEntry = (ReferenceFeatureEntry) runEntry.get(0);
        SampleAnalysisEntry sampleAnalysisEntry = (SampleAnalysisEntry) runEntry.get(1);
        variantContentTable.addSampleToReference(referenceFeatureEntry.name, sampleAnalysisEntry.sampleName);
        executor.execute(
            new SampleAnalyserRunnable(
                sampleAnalysisEntry.sampleName,
                referenceFeatureEntry,
                sampleAnalysisEntry.vcfFileReader
                    .query(referenceFeatureEntry.referenceSequenceLocation,
                        referenceFeatureEntry.locationStart,
                        referenceFeatureEntry.locationEnd),
                variantContentTable,
                arguments.getMinCoverage(),
                arguments.getMinFrequency(),
                arguments.getMinHet(),
                arguments.getMaxHet(),
                arguments.getMinQuality(),
                progress2
            )
        );
      }
      executor.shutdown();
      //noinspection ResultOfMethodCallIgnored
      executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
      progress2.setExtraMessage(Logging.getDoneMessage());
      return variantContentTable;
    }
  }

  /**
   * Adds information from the reference genome specified by the single elements of referenceFeatureEntries to the passed
   * variantContentTable.
   *
   * @param variantContentTable     {@link VariantContentTable} containing information about the variant sites of a
   *                                set of samples against possibly multiple reference loci.
   * @param referenceFeatureEntries Set of {@link ReferenceFeatureEntry} instances, each specifying a reference locus.
   * @param progress                {@link ProgressBar} instance to visualize progress information to the user.
   */
  private static void addReferenceToVariantContentTable(
      VariantContentTable variantContentTable, HashSet<ReferenceFeatureEntry> referenceFeatureEntries,
      ProgressBar progress) {
    for (ReferenceFeatureEntry referenceFeatureEntry : referenceFeatureEntries) {
      char[] referenceSequenceBases = referenceFeatureEntry.getReferenceSequence().toCharArray();
      int referenceSequenceStart = referenceFeatureEntry.locationStart;
      for (int i = 0; i < referenceSequenceBases.length; i++) {
        int position = referenceSequenceStart + i;
        variantContentTable
            .putVariantPosition(referenceFeatureEntry.name, "Reference", String.valueOf(position),
                referenceSequenceBases[i], Double.NaN, Double.NaN, Double.NaN, true,
                new HashMap<>());
      }
      progress.step();
    }
  }

}
