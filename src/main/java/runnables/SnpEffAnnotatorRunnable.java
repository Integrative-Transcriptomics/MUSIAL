package runnables;

import java.io.File;
import me.tongfei.progressbar.ProgressBar;
import tools.ArgumentsParser;
import utility.CL;

/**
 * Implementation of the {@link Runnable} interface to run SnpEff.
 * <p>
 * Runs SnpEff on one sample in a single thread.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class SnpEffAnnotatorRunnable implements Runnable {

  /**
   * Path to which SnpEff was extracted.
   */
  private final String snpEffPath;
  /**
   * Path to which SnpEff results are written.
   */
  private final String snpEffResultsPath;
  /**
   * The internal name of the sample to analyze.
   */
  private final String sampleName;
  /**
   * The index of the sample, i.e. at which index the samples input File and Name are stored in the parsed arguments.
   */
  private final int sampleIndex;
  /**
   * The name of the reference.
   */
  private final String referenceName;
  /**
   * {@link File} instance pointing to the samples input `.vcf` file.
   */
  private final File sampleInput;
  /**
   * Arguments parsed from the command line.
   */
  private final ArgumentsParser arguments;
  /**
   * An {@link ProgressBar} instance to output user information during runtime.
   */
  private final ProgressBar progress;

  /**
   * Constructor of {@link SnpEffAnnotatorRunnable}.
   *
   * @param snpEffPath        {@link String} representing the path at which SnpEff runnable is located.
   * @param snpEffResultsPath {@link String} representing the path to which results should be written.
   * @param sampleName        {@link String} name of the sample.
   * @param sampleIndex       {@link Integer} index at which the sample name and input are stored in the lists of the passed
   *                          {@link ArgumentsParser} instance.
   * @param referenceName     {@link String} name of the reference to which respect the sample is analyzed.
   * @param sampleInput       {@link File} pointing to the samples input `.vcf` file.
   * @param arguments         {@link ArgumentsParser} instance, contains arguments parsed from command line.
   * @param progress          {@link ProgressBar} instance to output user information during runtime.
   */
  public SnpEffAnnotatorRunnable(String snpEffPath, String snpEffResultsPath, String sampleName, int sampleIndex,
                                 String referenceName, File sampleInput, ArgumentsParser arguments,
                                 ProgressBar progress) {
    this.snpEffPath = snpEffPath;
    this.snpEffResultsPath = snpEffResultsPath;
    this.sampleName = sampleName;
    this.sampleIndex = sampleIndex;
    this.referenceName = referenceName;
    this.sampleInput = sampleInput;
    this.arguments = arguments;
    this.progress = progress;
  }

  /**
   * Runs the threaded analysis of the data specified with the instances properties.
   */
  @Override
  public void run() {
    // String array representation of a command line command.
    String[] runSNPEff = {"java", "-jar", "snpEff.jar", "-v", "-s", snpEffResultsPath + "/" + sampleName +
        "_snpEff_summary" +
        ".html",
        referenceName, sampleInput.getAbsolutePath()};
    // Definition of the output file name.
    String sampleOutputPath = snpEffResultsPath + "/" + sampleName + "_snpEff_annotated.vcf";
    /*
    Runs the specified command.
     */
    CL.runCommand(runSNPEff, snpEffResultsPath + "/" + sampleName + "_snpEff.log",
        sampleOutputPath,
        snpEffPath);
    // Replace the original input `.vcf` file with the annotated version.
    arguments.setSampleInput(new File(sampleOutputPath), sampleIndex);
    progress.step();
  }
}
