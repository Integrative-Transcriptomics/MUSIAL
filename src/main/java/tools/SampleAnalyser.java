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
 *
 */
public final class SampleAnalyser {

  /**
   * @param referenceEntries
   * @param sampleEntries
   * @param arguments
   * @return
   * @throws InterruptedException
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
      executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
      addReferenceToVariablePositionsTable(variablePositionsTable, referenceEntries, progress);
      progress.setExtraMessage(Logging.getDoneMessage());
      return variablePositionsTable;
    }
  }

  /**
   * @param variablePositionsTable
   * @param referenceEntries
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
