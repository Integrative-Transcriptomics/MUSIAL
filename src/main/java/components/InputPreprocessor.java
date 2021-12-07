package tools;

import datastructure.FastaEntry;
import datastructure.FeatureAnalysisEntry;
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
 * Comprises static methods to generate sets of {@link FeatureAnalysisEntry} and {@link SampleAnalysisEntry}
 * instances.
 * <p>
 * Each {@link FeatureAnalysisEntry} defines one genomic location for downstream analysis: These may comprise single genes
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
   * @return A {@link HashSet} of {@link FeatureAnalysisEntry} instances.
   * @throws IOException If any access to the reference `.fasta` file failed.
   */
  public static HashSet<FeatureAnalysisEntry> preprocessReferenceInput(ArgumentsParser arguments)
      throws IOException, MusialBioException {
    ProgressBar progress = Musial.buildProgress();
    try (progress) {
      File referenceSequenceFile = arguments.getReferenceFile();
      ArrayList<FeatureAnalysisEntry> geneFeatures = arguments.getIncludedGeneFeatures();
      HashSet<FeatureAnalysisEntry> analysisEntries = new HashSet<>();
      HashSet<FastaEntry> fastaEntries = IO.readFasta(referenceSequenceFile);
      if (geneFeatures.size() == 0) {
        // If no gene features were specified by the user, the reference fasta is parsed and each entry is processed
        // into one ReferenceAnalysisEntry.
        progress.maxHint(fastaEntries.size());
        for (FastaEntry fastaEntry : fastaEntries) {
          String entryName = adjustFastaHeader(fastaEntry.getHeader());
          String entryLocation = fastaEntry.getHeader().split(" ")[0].trim();
          String entrySequence = fastaEntry.getSequence();
          FeatureAnalysisEntry featureAnalysisEntry = new FeatureAnalysisEntry(entryName, entryLocation, 1,
              entrySequence.length());
          featureAnalysisEntry.setReferenceSequence(entrySequence);
          analysisEntries.add(featureAnalysisEntry);
          progress.step();
        }
      } else {
        // Else in each fasta entry it is searched for each specified gene feature. Finally for each gene feature one
        // ReferenceAnalysisEntry is adjusted.
        progress.maxHint((long) fastaEntries.size() * geneFeatures.size());
        for (FastaEntry fastaEntry : fastaEntries) {
          String entryLocation = fastaEntry.getHeader().split(" ")[0].trim();
          for (FeatureAnalysisEntry geneFeature : geneFeatures) {
            if (entryLocation.equals(geneFeature.referenceSequenceLocation)) {
              String referenceFeatureSequence =
                  fastaEntry.getSequence(geneFeature.locationStart, geneFeature.locationEnd);
              geneFeature.setReferenceSequence(referenceFeatureSequence);
              analysisEntries.add(geneFeature);
            }
            progress.step();
          }
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
      progress.maxHint(arguments.getSampleNames().size());
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

  /**
   * Adjusts a fasta entry header in such a way that no symbols that are not allowed for file names are contained
   * anymore.
   *
   * @param header {@link String} representing a header line parsed from a `.fasta` file.
   * @return A adjusted fasta header than can be used safely as filename.
   */
  private static String adjustFastaHeader(String header) {
    return header
        .trim()
        .substring(1)
        .replace('/', '_')
        .replace('\\', '_')
        .replace('\0', '_')
        .replace(',', '_')
        .replace('.', '_')
        .replace(':', '_')
        .replace('?', '_')
        .replace('*', '_')
        .replace('<', '_')
        .replace('>', '_')
        .replace('"', '_')
        .replace('|', '_');
  }

}
