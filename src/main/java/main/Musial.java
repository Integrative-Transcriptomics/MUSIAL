package main;

import datastructure.FeatureAnalysisEntry;
import datastructure.SampleAnalysisEntry;
import datastructure.VariablePositionsTable;
import exceptions.MusialCLAException;
import exceptions.MusialDuplicationException;
import exceptions.MusialFaultyDataException;
import exceptions.MusialIOException;
import java.io.IOException;
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
import org.apache.commons.cli.ParseException;
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
 */
public final class Musial {

  /**
   * The class name of the software, i.e. MUSIAL.
   */
  public static String CLASS_NAME = "";
  /**
   * The version of the software.
   */
  public static String VERSION = "";
  /**
   * Original out stream.
   */
  public static final PrintStream ORIGINAL_OUT_STREAM = System.out;
  /**
   * Original error stream.
   */
  public static final PrintStream ORIGINAL_ERR_STREAM = System.err;
  /**
   * Alternative output stream to dump log massages.
   */
  public static final PrintStream EMPTY_STREAM = new PrintStream(new OutputStream() {
    public void write(int b) {
    }
  });
  /**
   * Factory to generate {@link ProgressBar} instances with predefined properties.
   */
  private static final ProgressBarBuilder progressBarBuilder =
      new ProgressBarBuilder().setStyle(ProgressBarStyle.ASCII).setMaxRenderedLength(75);

  /**
   * Conducts the main pipeline of the tool (see comments and single methods for more detail).
   *
   * @param args Arguments parsed from the command line.
   */
  public static void main(String[] args) {
    try {
      // 1. Load metadata from resource files.
      loadMetadata();
      // 2. Parse command line arguments.
      ArgumentsParser arguments = parseCLArguments(args);
      // 3. (Optional) Run SnpEff.
      if (Objects.requireNonNull(arguments).isRunSnpEff()) {
        runSnpEff(arguments);
      }
      // 4. Process input reference sequence and (optional) gene features.
      HashSet<FeatureAnalysisEntry> referenceEntries = preprocessReferenceInput(arguments);
      // 5. Generate a pool of VCFFileReader instances, one for each sample.
      CopyOnWriteArraySet<SampleAnalysisEntry> sampleEntries = preprocessSampleInput(arguments);
      // 6. Run the analysis
      VariablePositionsTable variablePositionsTable = invokeSampleAnalyserRunner(referenceEntries, sampleEntries,
          arguments);
      // 7. Generate output.
      writeOutput(arguments, variablePositionsTable, referenceEntries);
    } catch (Exception e) {
      e.printStackTrace();
      Logging.logError(e.getCause() + ": " + e.getMessage());
    }
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
  private static ArgumentsParser parseCLArguments(String[] args)
      throws MusialCLAException, MusialIOException, ParseException, MusialFaultyDataException,
      MusialDuplicationException {
    return new ArgumentsParser(args);
  }

  /**
   * Runs SnpEff (additional annotation of .vcf files)
   * <p>
   * The command line arguments specified by the user and accepted by the `main` method are passed to the static method
   * of the {@link SnpEffAnnotator} class. The method does not return any value, but all annotated .vcf files generated
   * with SnpEff will be stored in the directory [output]/snpEff/results. If any errors arise during the process the
   * program exits.
   *
   * @param arguments Arguments parsed from the command line.
   */
  private static void runSnpEff(ArgumentsParser arguments) throws MusialIOException, IOException, InterruptedException {
    Logging.logStatus("Running SnpEff");
    SnpEffAnnotator.runSnpEff(arguments);
  }

  /**
   * Parses the specified reference data and generates a set of analysis entries.
   * <p>
   * Each {@link FeatureAnalysisEntry} represents a genomic region of the reference data, for example a single gene, contig,
   * plasmid or full genome. Each such entry contains naming as well as reference sequence information.
   *
   * @param arguments Arguments parsed from the command line.
   * @return A {@link HashSet} containing {@link FeatureAnalysisEntry} instances specifying genomic regions of the specified
   * reference data that are subject to further analysis.
   */
  private static HashSet<FeatureAnalysisEntry> preprocessReferenceInput(ArgumentsParser arguments)
      throws IOException, MusialFaultyDataException {
    Logging.logStatus("Preparing reference data for analysis");
    return InputPreprocessor.preprocessReferenceInput(arguments);
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
  private static CopyOnWriteArraySet<SampleAnalysisEntry> preprocessSampleInput(ArgumentsParser arguments)
      throws InterruptedException {
    Logging.logStatus("Preparing sample data for analysis");
    return InputPreprocessor.preprocessSampleInput(arguments);
  }

  /**
   * Invokes multi-threaded analysis of the input samples.
   * <p>
   * The single steps are implemented in the {@link SampleAnalyser} and {@link runnables.SampleAnalyserRunnable}
   * classes.
   *
   * @param analysisEntries       A set of {@link FeatureAnalysisEntry}.
   * @param sampleAnalysisEntries A set of {@link SampleAnalysisEntry}.
   * @param arguments             Arguments parsed from the command line.
   * @return A {@link VariablePositionsTable} containing the information processed from the specified input.
   */
  private static VariablePositionsTable invokeSampleAnalyserRunner(HashSet<FeatureAnalysisEntry> analysisEntries,
                                                                   CopyOnWriteArraySet<SampleAnalysisEntry> sampleAnalysisEntries,
                                                                   ArgumentsParser arguments)
      throws InterruptedException {
    Logging.logStatus("Analysing samples");
    return SampleAnalyser.run(analysisEntries, sampleAnalysisEntries, arguments);
  }

  /**
   * Generate output files.
   * <p>
   * The single steps are implemented in the {@link OutputWriter} and {@link utility.IO} classes.
   *
   * @param arguments              Arguments parsed from the command line.
   * @param variablePositionsTable The {@link VariablePositionsTable} containing the information to write output with.
   * @param featureAnalysisEntries represents a genomic region of the reference data, for example a single gene, contig,
   *                               plasmid or full genome. Each such entry contains naming as well as reference sequence information.
   */
  private static void writeOutput(ArgumentsParser arguments, VariablePositionsTable variablePositionsTable,
                                  HashSet<FeatureAnalysisEntry> featureAnalysisEntries)
      throws MusialIOException {
    Logging.logStatus("Generating output");
    OutputWriter.writeOutput(arguments, variablePositionsTable, featureAnalysisEntries);
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
