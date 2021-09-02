package tools;

import com.google.common.collect.Sets;
import datastructure.ReferenceAnalysisEntry;
import datastructure.SampleAnalysisEntry;
import datastructure.VariablePosition;
import datastructure.VariablePositionsTable;
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
 * {@link ReferenceAnalysisEntry} and a set of {@link SampleAnalysisEntry} instances is generated. For more detail
 * view the respective class documentations.
 * <p>
 * This results in a set of pairs of {@link ReferenceAnalysisEntry} and {@link SampleAnalysisEntry} objects, each
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
   * @param referenceEntries Set of {@link ReferenceAnalysisEntry} instances. Each specifies a locus on the reference
   *                         genome.
   * @param sampleEntries    Set of {@link SampleAnalysisEntry} instances. Each specifies the `.vcf` file content from
   *                         one sample.
   * @param arguments        {@link ArgumentsParser} containing arguments parsed from command line.
   * @return {@link VariablePositionsTable} containing the information returned from each single
   * {@link SampleAnalyserRunnable}.
   */
  public static VariablePositionsTable run(HashSet<ReferenceAnalysisEntry> referenceEntries,
                                           CopyOnWriteArraySet<SampleAnalysisEntry> sampleEntries,
                                           ArgumentsParser arguments)
      throws InterruptedException {
    ProgressBar progress = Musial.buildProgress();
    try (progress) {
      ExecutorService executor = Executors.newFixedThreadPool(arguments.getNumThreads());
      Set<List<Object>> runEntries = Sets.cartesianProduct(referenceEntries, sampleEntries);
      VariablePositionsTable variablePositionsTable = new VariablePositionsTable();
      progress.maxHint(runEntries.size() + referenceEntries.size());
      for (List<Object> runEntry : runEntries) {
        ReferenceAnalysisEntry referenceAnalysisEntry = (ReferenceAnalysisEntry) runEntry.get(0);
        SampleAnalysisEntry sampleAnalysisEntry = (SampleAnalysisEntry) runEntry.get(1);
        executor.execute(
            new SampleAnalyserRunnable(
                sampleAnalysisEntry.sampleName,
                referenceAnalysisEntry,
                sampleAnalysisEntry.vcfFileReader
                    .query(referenceAnalysisEntry.analysisSequenceLocation,
                        referenceAnalysisEntry.analysisSequenceStart,
                        referenceAnalysisEntry.analysisSequenceEnd),
                variablePositionsTable,
                arguments.getMinCoverage(),
                arguments.getMinFrequency(),
                arguments.getMinQuality(),
                progress
            )
        );
      }
      executor.shutdown();
      //noinspection ResultOfMethodCallIgnored
      executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
      addReferenceToVariablePositionsTable(variablePositionsTable, referenceEntries, progress);
      progress.setExtraMessage(Logging.getDoneMessage());
      return variablePositionsTable;
    }
  }

  /**
   * Adds information from the reference genome specified by the single elements of referenceEntries to the passed
   * variablePositionsTable.
   *
   * @param variablePositionsTable {@link VariablePositionsTable} containing information about the variant sites of a
   *                               set of samples against possibly multiple reference loci.
   * @param referenceEntries       Set of {@link ReferenceAnalysisEntry} instances, each specifying a reference locus.
   */
  private static void addReferenceToVariablePositionsTable(
      VariablePositionsTable variablePositionsTable, HashSet<ReferenceAnalysisEntry> referenceEntries,
      ProgressBar progress) {
    for (ReferenceAnalysisEntry referenceEntry : referenceEntries) {
      char[] referenceSequenceBases = referenceEntry.analysisSequence.toCharArray();
      int referenceSequenceStart = referenceEntry.analysisSequenceStart;
      for (int i = 0; i < referenceSequenceBases.length; i++) {
        variablePositionsTable.putVariablePosition(
            referenceEntry.analysisIdentifier,
            "Reference",
            String.valueOf(referenceSequenceStart + i),
            new VariablePosition(
                referenceSequenceBases[i],
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY,
                "",
                1.0
            )
        );
      }
      progress.step();
    }
  }

}
