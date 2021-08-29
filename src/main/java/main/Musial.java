package main;

import datastructure.ReferenceAnalysisEntry;
import datastructure.SampleAnalysisEntry;
import datastructure.VariablePositionsTable;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Objects;
import java.util.Properties;
import java.util.concurrent.CopyOnWriteArraySet;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
import tools.ArgumentsParser;
import tools.InputPreprocessor;
import tools.OutputWriter;
import tools.SampleAnalyser;
import tools.SnpEffAnnotator;
import utility.Logging;

/**
 * Main class of MUSIAL (MUlti Sample varIant AnaLysis), a tool to calculate SNV, gene, and whole genome alignments,
 * together with other relevant statistics based on vcf files.
 *
 * @author Alexander Seitz
 * @author Simon Hackl
 * @version 2.0
 * @TODO: 25.08.2021 : Genes from negative sense strand.
 * @TODO: 25.08.2021 : Missing components; Heterozygous calling, phylogenetics, SnpEff output.
 * @TODO: 27.08.2021 : Change output if multiple quality criteria are not met?
 * @TODO: 27.08.2021 : Multi threading for single genomic location?
 * @TODO: 28.08.2021 : Which fields to use for snv statistics.
 * @TODO: 28.08.2021 : Remove per position per sample frequency. Does not make any sense?
 */
public final class Musial {

  public static String CLASS_NAME = "";
  public static String VERSION = "";
  public static final PrintStream ORIGINAL_OUT_STREAM = System.out;
  public static final PrintStream ORIGINAL_ERR_STREAM = System.err;
  public static final PrintStream EMPTY_STREAM = new PrintStream(new OutputStream() {
    public void write(int b) {
    }
  });
  private static final ProgressBarBuilder progressBarBuilder =
      new ProgressBarBuilder().setStyle(ProgressBarStyle.ASCII).setMaxRenderedLength(75);

  /**
   * Conducts the main pipeline of the tool (see comments and single methods for more detail).
   *
   * @param args Arguments parsed from the command line.
   */
  public static void main(String[] args) {
    // 1. Load metadata from resource files.
    loadMetadata();
    // 2. Parse command line arguments.
    ArgumentsParser arguments = parseCLArguments(args);
    // 3. (Optional) Run SnpEff.
    if (Objects.requireNonNull(arguments).isRunSnpEff()) {
      runSnpEff(arguments);
    }
    // 4. Process input reference sequence and (optional) gene features.
    HashSet<ReferenceAnalysisEntry> referenceEntries = preprocessReferenceInput(arguments);
    // 5. Generate a pool of VCFFileReader instances, one for each sample.
    CopyOnWriteArraySet<SampleAnalysisEntry> sampleEntries = preprocessSampleInput(arguments);
    // 6. Run the analysis
    VariablePositionsTable variablePositionsTable = invokeSampleAnalyserRunner(referenceEntries, sampleEntries,
        arguments);
    // 7. Generate output.
    writeOutput(arguments, variablePositionsTable);
    // (DEV ONLY) Check memory usage.
    Runtime runtime = Runtime.getRuntime();
    runtime.gc();
    long memory = runtime.totalMemory() - runtime.freeMemory();
    System.out.println("Used memory is megabytes: " + (memory / (1024L * 1024L)));
    System.exit(0);
  }

  /**
   * Loads metadata, such as the projects title and version from resources and prints the information to stdout.
   */
  private static void loadMetadata() {
    Properties properties = new Properties();
    try {
      // Load title from resources.
      InputStream in = Musial.class.getResourceAsStream("/title.properties");
      properties.load(in);
      Musial.CLASS_NAME = properties.getProperty("title");
      // Load version from resources.
      in = Musial.class.getResourceAsStream("/version.properties");
      properties.load(in);
      Musial.VERSION = properties.getProperty("version");
      assert in != null;
      // Pretty print information into stdout.
      System.out.println("=".repeat(18));
      System.out.println(" ".repeat(3) + " " + Musial.CLASS_NAME + " " + Musial.VERSION + " " + " ".repeat(3));
      System.out.println("=".repeat(18));
    } catch (Exception e) {
      Logging.logError(e.getCause() + ": " + e.getMessage());
    }
  }

  /**
   * Parses command line arguments by instantiating a {@link ArgumentsParser} object.
   * <p>
   * The command line arguments specified by the user and accepted by the `main` method are passed and used to
   * instantiate an {@link ArgumentsParser} object which is returned if all argument validation steps were successful.
   * If any error arises during the parsing step the program exits.
   *
   * @param args Arguments parsed from the command line.
   * @return An instance of the {@link ArgumentsParser} class.
   */
  private static ArgumentsParser parseCLArguments(String[] args) {
    try {
      return new ArgumentsParser(args);
    } catch (Exception e) {
      Logging.logError(e.getCause() + ": " + e.getMessage());
    }
    return null;
  }

  /**
   * Runs SnpEff (additional annotation of .vcf files)
   * <p>
   * The command line arguments specified by the user and accepted by the `main` method are passed to the static method
   * of the {@link SnpEffAnnotator} class. The method does not return any value, but all annotated .vcf files generated
   * with SnpEff will be stored in the directory <output>/snpEff/results. If any errors arise during the process the
   * program exits.
   *
   * @param arguments Arguments parsed from the command line.
   */
  private static void runSnpEff(ArgumentsParser arguments) {
    Logging.logStatus("Running SnpEff");
    try {
      SnpEffAnnotator.runSnpEff(arguments);
    } catch (Exception e) {
      Logging.logError(e.getCause() + ": " + e.getMessage());
    }
  }

  /**
   * Parses the specified reference data and generates a set of analysis entries.
   * <p>
   * Each {@link ReferenceAnalysisEntry} represents a genomic region of the reference data, for example a single gene, contig,
   * plasmid or full genome. Each such entry contains naming as well as reference sequence information.
   *
   * @param arguments Arguments parsed from the command line.
   * @return A {@link HashSet} containing {@link ReferenceAnalysisEntry} instances specifying genomic regions of the specified
   * reference data that are subject to further analysis.
   */
  private static HashSet<ReferenceAnalysisEntry> preprocessReferenceInput(ArgumentsParser arguments) {
    Logging.logStatus("Preparing reference data for analysis");
    try {
      return InputPreprocessor.preprocessReferenceInput(arguments);
    } catch (Exception e) {
      Logging.logError(e.getCause() + ": " + e.getMessage());
    }
    return null;
  }

  /**
   * Analyzes the specified sample data and generates a set of analysis entries.
   * <p>
   * Each {@link SampleAnalysisEntry} represents a sample name and a {@link htsjdk.variant.vcf.VCFFileReader}
   * instance pointing to the specified input `.vcf` file of the sample.
   *
   * @param arguments Arguments parsed from the command line.
   * @return A {@link CopyOnWriteArraySet} containing {@link SampleAnalysisEntry} instances, each specifying the
   * input for one sample to analyze.
   */
  private static CopyOnWriteArraySet<SampleAnalysisEntry> preprocessSampleInput(ArgumentsParser arguments) {
    Logging.logStatus("Preparing sample data for analysis");
    try {
      return InputPreprocessor.preprocessSampleInput(arguments);
    } catch (Exception e) {
      Logging.logError(e.getCause() + ": " + e.getMessage());
    }
    return null;
  }

  /**
   * Invokes multi-threaded analysis of the input samples.
   * <p>
   * The single steps are implemented in the {@link SampleAnalyser} and {@link runnables.SampleAnalyserRunnable}
   * classes.
   *
   * @param analysisEntries       A set of {@link ReferenceAnalysisEntry}.
   * @param sampleAnalysisEntries A set of {@link SampleAnalysisEntry}.
   * @param arguments             Arguments parsed from the command line.
   * @return A {@link VariablePositionsTable} containing the information processed from the specified input.
   */
  private static VariablePositionsTable invokeSampleAnalyserRunner(HashSet<ReferenceAnalysisEntry> analysisEntries,
                                                                   CopyOnWriteArraySet<SampleAnalysisEntry> sampleAnalysisEntries,
                                                                   ArgumentsParser arguments) {
    Logging.logStatus("Analysing samples");
    try {
      return SampleAnalyser.run(analysisEntries, sampleAnalysisEntries, arguments);
    } catch (Exception e) {
      Logging.logError(e.getCause() + ": " + e.getMessage());
    }
    return null;
  }

  /**
   * Generate output files.
   * <p>
   * The single steps are implemented in the {@link OutputWriter} and {@link utility.IO} classes.
   *
   * @param arguments              Arguments parsed from the command line.
   * @param variablePositionsTable The {@link VariablePositionsTable} containing the information to write output with.
   */
  private static void writeOutput(ArgumentsParser arguments, VariablePositionsTable variablePositionsTable) {
    Logging.logStatus("Generating output");
    try {
      OutputWriter.writeOutput(arguments, variablePositionsTable);
    } catch (Exception e) {
      e.printStackTrace();
      Logging.logError(e.getMessage());
    }
  }

  /**
   * Generates and returns {@link ProgressBar} instance.
   *
   * @return Instance of {@link ProgressBar} with the options specified by the progressBarBuilder property.
   */
  public static ProgressBar buildProgress() {
    return progressBarBuilder.build();
  }

}
