package components;

import datastructure.FastaEntry;
import datastructure.ReferenceFeatureEntry;
import datastructure.SampleAnalysisEntry;
import exceptions.MusialBioException;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.CopyOnWriteArraySet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import main.Musial;
import me.tongfei.progressbar.ProgressBar;
import runnables.SamplePreprocessorRunnable;
import utility.IO;
import utility.Logging;

/**
 * Comprises static methods to generate sets of {@link ReferenceFeatureEntry} and {@link SampleAnalysisEntry}
 * instances.
 * <p>
 * Each {@link ReferenceFeatureEntry} defines one genomic location for downstream analysis: These may comprise single genes
 * on different chromosomes or plasmids, contigs or full genomes.
 * <p>
 * Each {@link SampleAnalysisEntry} defines a input sample by specifying the sample name and input `.vcf` file
 * {@link htsjdk.variant.vcf.VCFFileReader}.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class InputPreprocessor {

  /**
   * Preprocesses the specified reference `.fasta` file.
   *
   * @param arguments {@link ArgumentsParser} storing arguments parsed from the command line.
   * @return A {@link HashSet} of {@link ReferenceFeatureEntry} instances.
   * @throws IOException        If any access to the reference `.fasta` file failed.
   * @throws MusialBioException If any faulty information is processed during reference feature entry generation.
   */
  public static HashSet<ReferenceFeatureEntry> preprocessReferenceInput(ArgumentsParser arguments)
      throws IOException, MusialBioException {
    ProgressBar progress = Musial.buildProgress();
    try (progress) {
      File referenceSequenceFile = arguments.getReferenceFile();
      ArrayList<ReferenceFeatureEntry> referenceFeatures = arguments.getReferenceFeatures();
      HashSet<ReferenceFeatureEntry> analysisEntries = new HashSet<>();
      HashSet<FastaEntry> fastaEntries;
      if (referenceFeatures.size() == 0) {
        // If no gene features were specified by the user, the reference fasta is parsed and each entry is processed
        // into one ReferenceAnalysisEntry.
        Logging.logStatus("Reading reference genome.");
        fastaEntries = IO.readFasta(referenceSequenceFile);
        progress.maxHint(fastaEntries.size());
        progress.stepTo(0);
        for (FastaEntry fastaEntry : fastaEntries) {
          String entryName = fastaEntry.getHeader().split(" ")[0].trim();
          String entrySequence = fastaEntry.getSequence();
          ReferenceFeatureEntry referenceFeatureEntry = new ReferenceFeatureEntry(entryName, entryName, false,
              entryName, 1,
              entrySequence.length());
          referenceFeatureEntry.setReferenceSequence(entrySequence);
          analysisEntries.add(referenceFeatureEntry);
          progress.step();
        }
      } else {
        // Else in each fasta entry it is searched for each specified gene feature. Finally for each gene feature one
        // ReferenceAnalysisEntry is adjusted.
        Logging.logStatus("Reading reference gene features (" + referenceFeatures.size() + ").");
        fastaEntries = IO.readFasta(referenceSequenceFile);
        progress.maxHint(referenceFeatures.size());
        progress.stepTo(0);
        for (ReferenceFeatureEntry referenceFeature : referenceFeatures) {
          for (FastaEntry fastaEntry : fastaEntries) {
            String entryLocation = fastaEntry.getHeader().split(" ")[0].trim();
            if (entryLocation.equals(referenceFeature.referenceSequenceLocation)) {
              String referenceFeatureSequence =
                  fastaEntry.getSequence(referenceFeature.locationStart, referenceFeature.locationEnd);
              referenceFeature.setReferenceSequence(referenceFeatureSequence);
              analysisEntries.add(referenceFeature);
            }
          }
          progress.step();
        }
      }
      progress.setExtraMessage(Logging.getDoneMessage());
      return analysisEntries;
    }
  }

  /**
   * Preprocesses the specified samples `.vcf` file.
   *
   * @param arguments {@link ArgumentsParser} storing arguments parsed from the command line.
   * @return A thread safe {@link CopyOnWriteArraySet} of {@link SampleAnalysisEntry} instances.
   * @throws InterruptedException If the process was interrupted.
   */
  public static CopyOnWriteArraySet<SampleAnalysisEntry> preprocessSampleInput(ArgumentsParser arguments)
      throws InterruptedException {
    ProgressBar progress = Musial.buildProgress();
    try (progress) {
      Logging.logStatus("Validating/generating sample index files (" + arguments.getSampleNames().size() + ").");
      progress.maxHint(arguments.getSampleNames().size());
      progress.stepTo(0);
      ExecutorService executor = Executors.newFixedThreadPool(arguments.getNumThreads());
      CopyOnWriteArraySet<SampleAnalysisEntry> vcfFileReaderPool = new CopyOnWriteArraySet<>();
      for (int i = 0; i < arguments.getSampleInput().size(); i++) {
        executor.execute(new SamplePreprocessorRunnable(
            arguments.getSampleInput().get(i),
            arguments.getSampleNames().get(i),
            vcfFileReaderPool,
            progress
        ));
      }
      executor.shutdown();
      //noinspection ResultOfMethodCallIgnored
      executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
      progress.setExtraMessage(Logging.getDoneMessage());
      return vcfFileReaderPool;
    }
  }
}
