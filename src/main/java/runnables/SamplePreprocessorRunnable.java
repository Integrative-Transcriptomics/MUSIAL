package runnables;

import datastructure.SampleAnalysisEntry;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.nio.file.Path;
import java.util.concurrent.CopyOnWriteArraySet;
import main.Musial;
import me.tongfei.progressbar.ProgressBar;
import utility.Logging;

/**
 * Implementation of the {@link Runnable} interface to preprocess a sample input.
 * <p>
 * Runs a sample preprocessing in a single thread.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class SamplePreprocessorRunnable implements Runnable {

  /**
   * {@link File} instance pointing to the samples input `.vcf` file.
   */
  private final File sampleInput;
  /**
   * The internal name of the sample.
   */
  private final String sampleName;
  /**
   * A thread safe {@link CopyOnWriteArraySet<SampleAnalysisEntry>} to store the preprocessed sample input.
   */
  private final CopyOnWriteArraySet<SampleAnalysisEntry> sampleEntries;
  /**
   * An {@link ProgressBar} instance to output user information during runtime.
   */
  private final ProgressBar progress;

  /**
   * Constructor of {@link SamplePreprocessorRunnable}.
   *
   * @param sampleInput   {@link File} instance pointing to the `.vcf` file of the respective sample.
   * @param sampleName    {@link String} the internal name of the sample.
   * @param sampleEntries {@link CopyOnWriteArraySet<SampleAnalysisEntry>} thread safe set to store results of the
   *                      analysis.
   * @param progress      {@link ProgressBar} to indicate runtime information.
   */
  public SamplePreprocessorRunnable(File sampleInput, String sampleName,
                                    CopyOnWriteArraySet<SampleAnalysisEntry> sampleEntries, ProgressBar progress) {
    this.sampleInput = sampleInput;
    this.sampleName = sampleName;
    this.sampleEntries = sampleEntries;
    this.progress = progress;
  }

  /**
   * Runs the threaded analysis of the data specified with the instances properties.
   */
  @Override
  public void run() {
    try {
      // Check for the existence of a `.tbi.gz` index file of the input `.vcf` file.
      if (!new File(sampleInput.getAbsolutePath() + ".tbi.gz").exists()) {
        // If none is present, an index is created...
        TabixIndex tabixIndex = IndexFactory.createTabixIndex(
            sampleInput,
            new VCFCodec(),
            null
        );
        // ...and written to the same directory as the input `.vcf` file.
        tabixIndex.write(Path.of(sampleInput.getAbsolutePath() + ".tbi.gz"));
      }
      // A VCFFileReader can now be initialized and used to create a SampleAnalysisEntry.
      sampleEntries.add(new SampleAnalysisEntry(
          new VCFFileReader(sampleInput, new File(sampleInput.getAbsolutePath() + ".tbi.gz")), sampleName
      ));
    } catch (Exception e) {
      if (Musial.debug) {
        e.printStackTrace();
      } else {
        Logging.logWarning("An error occurred during sample preprocessing: " + e.getMessage());
      }
    } finally {
      progress.step();
    }
  }

}
