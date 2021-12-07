package tools;

import datastructure.FeatureAnalysisEntry;
import exceptions.MusialBioException;
import exceptions.MusialCLAException;
import exceptions.MusialIOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import main.Musial;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.io.FilenameUtils;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;
import utility.IO;
import utility.Validation;

/**
 * Parser for command line arguments.
 * <p>
 * Used to parse, validate and access all passed command line options or print help information to the user.
 *
 * @author Simon Hackl
 * @version 2.0
 * @since 2.0
 */
public final class ArgumentsParser {

  /**
   * Class name property of {@link Musial} used to format help message.
   */
  private static final String CLASS_NAME = Musial.CLASS_NAME + " [OPTIONS]";
  /**
   * The specified command line arguments.
   */
  private final String[] arguments;
  /**
   * Number of threads to use for computations.
   */
  private int numThreads = 1;
  /**
   * Minimum (read) coverage in order to call a SNV.
   */
  private Double minCoverage = 5d;
  /**
   * Minimum allele frequency (with respect to reads) in order to call a SNV.
   */
  private Double minFrequency = 0.9;
  /**
   * Minimum (mapping?) quality in order to call a SNV.
   */
  private Double minQuality = 30d;
  /**
   * {@link File} object representing the reference genome sequence .fasta file.
   */
  private final File referenceFile;
  /**
   * {@link ArrayList} of {@link File} objects pointing to .vcf files, one for each input sample.
   */
  private final List<File> sampleInput = Collections.synchronizedList(new ArrayList<>());
  /**
   * {@link ArrayList} of {@link String} objects representing the names for each sample. The i-th entry corresponds
   * to the name of the i-th .vcf file represented in the `sampleInput` property.
   */
  private final List<String> sampleNames = Collections.synchronizedList(new ArrayList<>());
  /**
   * {@link File} object representing the reference genome annotation .gff file.
   */
  private File annotationInput;
  /**
   * {@link File} object representing the output directory.
   */
  private final File outputDirectory;
  /**
   * Whether SNPEff is called on each input samples .vcf file before further analysis.
   */
  private boolean runSnpEff = false;
  /**
   * {@link ArrayList} of {@link FeatureAnalysisEntry} objects representing gene names, their parent sequence and locations,
   * that are exclusively analysed instead of the whole genome.
   */
  private final ArrayList<FeatureAnalysisEntry> includedGeneFeatures = new ArrayList<>();
  /**
   * Minimum allele frequency (with respect to reads) in order to call heterozygous SNVs.
   */
  private Double minHet = 0.45;
  /**
   * Maximum allele frequency (with respect to reads) in order to call heterozygous SNVs.
   */
  private Double maxHet = 0.55;
  /**
   * {@link HashMap} mapping {@link String} objects representing gene names to {@link File} objects pointing to `
   * .pdb` files, one for each gene name.
   */
  private final HashMap<String, File> pdInputFiles = new HashMap<>();

  /**
   * Constructor of the ArgumentParser class.
   * <p>
   * The ArgumentParser is used to parse and validate the user command line input for Musial. Parsed arguments are
   * stored - and accessible - via class properties.
   *
   * @param args {@link String} {@link Array} containing the command line arguments.
   * @throws MusialCLAException         If any error occurs during parsing or validating an
   *                                    {@link MusialCLAException} is thrown.
   * @throws ParseException             If any parsing error occurs.
   * @throws MusialBioException  If any faulty data was specified (i.e. features on the reference sequence with
   *                                    length 0).
   */
  public ArgumentsParser(String[] args)
      throws MusialCLAException, ParseException, MusialBioException,
      MusialIOException {
    // Store original command line arguments.
    this.arguments = args;
    // Add option to print help message.
    Options helpOptions = new Options();
    helpOptions.addOption("h", "help", false, "Display help information.");
    // Add non-help options.
    Options options = new Options();
    options.addOption("h", "help", false, "Display help information.");
    // Add options used to specify the filtering of SNVs.
    options.addOption("nt", "numThreads", true, "Number of threads to use [" + this.numThreads + "]");
    options.addOption("mc", "minCoverage", true, "Minimum coverage to call a SNV [" + this.minCoverage + "]");
    options.addOption("mf", "minFrequency", true, "Minimum (read) frequency to call a SNV [" + this.minFrequency + "]");
    options.addOption("mq", "minQual", true, "Minimum quality to call a SNV [" + this.minQuality + "]");
    // Add options used to specify input data.
    options.addOption(Option.builder("r")
        .longOpt("reference")
        .desc("Path to a .fasta file; The reference genome sequence.")
        .required()
        .hasArg()
        .build());
    options.addOption(Option.builder("s")
        .longOpt("sampleSpecification")
        .desc("Text file specifying the sample input files and sample meta-information. Each line has to contain at " +
            "least the path to the samples .vcf file and may contain the following optional comma-separated fields " +
            "in the given order: <NAME> (a string representing the name to use for the sample).")
        .hasArg()
        .build());
    options.addOption(Option.builder("a")
        .longOpt("referenceAnnotation")
        .desc("Path to a .gff file; The reference genome feature annotations.")
        .hasArg()
        .hasArgs()
        .build());
    options.addOption(Option.builder("o")
        .longOpt("output")
        .desc("Path to the directory at which results files are generated.")
        .required()
        .hasArg()
        .build());
    // Add options to specify optional behaviour/additional computations of MUSIAL.
    options.addOption("sf", "snpEff", false, "Run SNPEff. [" + this.runSnpEff + "]");
    options.addOption(Option.builder("gf")
        .longOpt("geneFeatures")
        .desc("Text file specifying gene features to analyze exclusively instead of the whole genome. Each line has " +
            "to contain at least the identifier to match from the specified .gff file and may contain the following " +
            "optional comma-separated fields in the given order: <QUERY> (a string used to query the reference " +
            "genome feature annotation to search for the features locus, should match the .gff NAME attribute), " +
            "<NAME> (a string representing the name to " +
            "use for the feature), <PDBPATH> (a string representing the path to a .pdb file specifying the protein " +
            "product of the feature).")
        .hasArg()
        .build());
    options.addOption(Option.builder("hf")
        .longOpt("heterozygousFrequencies")
        .desc("If set, heterozygous calls will be respected by accounting for read frequencies in the set interval. " +
            "Expected " +
            " min,max. [" + this.minHet + "," + this.maxHet +
            "]")
        .numberOfArgs(2)
        .valueSeparator(',')
        .build());
    // Instantiate a formatter for the help message and a default command line parser.
    HelpFormatter helpformatter = new HelpFormatter();
    CommandLineParser parser = new DefaultParser();
    CommandLine cmd;
    // First it is checked if the user specified only the -h option.
    try {
      cmd = parser.parse(helpOptions, args);
      if (cmd.hasOption('h')) {
        printHelp(options, helpformatter);
        System.exit(1);
      }
    } catch (ParseException ignored) {
    }
    // If the user has specified the -h option together with other options the same behaviour is applied.
    try {
      cmd = parser.parse(options, args);
      if (cmd.hasOption('h')) {
        printHelp(options, helpformatter);
        System.exit(1);
      }
    } catch (ParseException ignored) {
    }
    // Else it is tried to parse and validate command line options.
    cmd = parser.parse(options, args);
    // In the following each possible command line argument is validated, if any faulty input is detected the
    // application will exit with an exception.
    // Validation of SNV filtering parameters:
    if (cmd.hasOption("nt")) {
      if (Validation.isPositiveInteger(cmd.getOptionValue("nt"))) {
        this.numThreads = Integer.parseInt(cmd.getOptionValue("nt"));
      } else {
        throw new MusialCLAException("`-nt` Number of threads is no positive integer:\t" + cmd.getOptionValue("nt"));
      }
    }
    if (cmd.hasOption("mc")) {
      if (Validation.isPositiveInteger(cmd.getOptionValue("mc"))) {
        this.minCoverage = (double) Integer.parseInt(cmd.getOptionValue("mc"));
      } else {
        throw new MusialCLAException("`-mc` Value for min. coverage is no positive integer:\t" + cmd.getOptionValue(
            "mc"));
      }
    }
    if (cmd.hasOption("mf")) {
      if (Validation.isPercentage(cmd.getOptionValue("mf"))) {
        this.minFrequency = Double.parseDouble(cmd.getOptionValue("mf"));
      } else {
        throw new MusialCLAException(
            "`-mf` Value for min. frequency is no double in the interval [0.0,1.0]:\t" + cmd.getOptionValue(
                "mf"));
      }
    }
    if (cmd.hasOption("mq")) {
      if (Validation.isPositiveDouble(cmd.getOptionValue("mq"))) {
        this.minQuality = Double.parseDouble(cmd.getOptionValue("mq"));
      } else {
        throw new MusialCLAException(
            "`-mq` Value for min. quality is no positive double:\t" + cmd.getOptionValue("mq"));
      }
    }
    // Validation of data input parameters:
    this.referenceFile = new File(cmd.getOptionValue('r'));
    if (Validation.validateInputFile(this.referenceFile)) {
      throw new MusialCLAException("`-r` The specified reference input file does not exist or has no read " +
          "permission:\t" + this.referenceFile.getAbsolutePath());
    }
    // Assert that the sample input option is specified.
    File sampleFile = new File(cmd.getOptionValue('s'));
    try {
      parseSampleInput(sampleFile);
    } catch (FileNotFoundException e) {
      throw new MusialCLAException(
          "`-s` The specified sample input file does not exist or has no read permission:\t" +
              sampleFile.getAbsolutePath());
    }
    if (cmd.hasOption("a")) {
      this.annotationInput = new File(cmd.getOptionValue("a"));
      if (Validation.validateInputFile(this.annotationInput)) {
        throw new MusialCLAException("`-a` The specified input file does not exist or has no read permission:\t" +
            this.annotationInput.getAbsolutePath());
      }
    }
    this.outputDirectory = new File(cmd.getOptionValue('o'));
    boolean outputDirectoryExists = this.outputDirectory.exists();
    if (!outputDirectoryExists) {
      outputDirectoryExists = this.outputDirectory.mkdirs();
    }
    if (!this.outputDirectory.isDirectory()) {
      throw new MusialCLAException("`-o` The specified output directory is no directory.");
    }
    if (Objects.requireNonNull(this.outputDirectory.list()).length > 0) {
      throw new MusialCLAException("`-o` The specified output directory is not empty.");
    }
    if (!outputDirectoryExists) {
      throw new MusialCLAException("`-o` The specified output directory does not exist and could not be created.");
    }
    // Validate optional behaviour options:
    if (cmd.hasOption("sf")) {
      if (!cmd.hasOption("a")) {
        throw new MusialCLAException(
            "`-sf` For SNPEff analysis a reference annotation has to be specified (see 'a' option).");
      } else {
        this.runSnpEff = true;
      }
    }
    if (cmd.hasOption("gf")) {
      if (!cmd.hasOption("a")) {
        throw new MusialCLAException(
            "`-gf` To use gene feature analysis a reference annotation has to be specified (see 'g' option).");
      }
      File featureFile = new File(cmd.getOptionValue("gf"));
      try {
        parseGeneFeatures(featureFile);
      } catch (FileNotFoundException e) {
        throw new MusialCLAException(
            "`-gf` The specified gene feature input file does not exist or has no read permission:\t" +
                featureFile.getAbsolutePath());
      }
    }
    if (cmd.hasOption("hf")) {
      String[] percentages = cmd.getOptionValues("hf");
      if (percentages != null) {
        if (percentages.length == 2) {
          if (Validation.isPercentage(percentages[0]) && Validation.isPercentage(percentages[1])) {
            double min = Double.parseDouble(percentages[0]);
            double max = Double.parseDouble(percentages[1]);
            if (min > max) {
              throw new MusialCLAException(
                  "`-hf` The specified minimum percentage is larger than the specified maximum percentage:\t" +
                      cmd.getOptionValue("hf"));
            }
            this.minHet = min;
            this.maxHet = max;
          } else {
            throw new MusialCLAException(
                "`-hf` One of the specified values is no percentage:\t" + cmd.getOptionValue("hf"));
          }
        } else {
          throw new MusialCLAException(
              "`-hf` Wrong number of arguments, expected two:\t" + cmd.getOptionValue("hf"));
        }
      }
    }
  }

  /**
   * @return The parsed number of threads.
   */
  public int getNumThreads() {
    return numThreads;
  }

  /**
   * @return The parsed min. coverage for variant calls.
   */
  public Double getMinCoverage() {
    return minCoverage;
  }

  /**
   * @return The parsed min. frequency for variant calls.
   */
  public Double getMinFrequency() {
    return minFrequency;
  }

  /**
   * @return The parsed min. quality for variant calls.
   */
  public Double getMinQuality() {
    return minQuality;
  }

  /**
   * @return The parsed `.fasta` reference input {@link File}.
   */
  public File getReferenceFile() {
    return referenceFile;
  }

  /**
   * @return The parsed `.vcf` sample input files as {@link List<File>}.
   */
  public List<File> getSampleInput() {
    return sampleInput;
  }

  /**
   * @return The parsed sample names as {@link List<String>}.
   */
  public List<String> getSampleNames() {
    return sampleNames;
  }

  /**
   * @return The parsed `.gff` input {@link File}.
   */
  public File getAnnotationInput() {
    return annotationInput;
  }

  /**
   * @return The specified output directory as {@link File}.
   */
  public File getOutputDirectory() {
    return outputDirectory;
  }

  /**
   * @return Whether SnpEff should be run.
   */
  public boolean isRunSnpEff() {
    return runSnpEff;
  }

  /**
   * @return The included gene features as {@link ArrayList<FeatureAnalysisEntry>}.
   */
  public ArrayList<FeatureAnalysisEntry> getIncludedGeneFeatures() {
    return includedGeneFeatures;
  }

  /**
   * @return The minimum (read) frequency to call a heterozygous SNV.
   */
  public Double getMinHet() {
    return minHet;
  }

  /**
   * @return The maximum (read) frequency to call a heterozygous SNV.
   */
  public Double getMaxHet() {
    return maxHet;
  }

  /**
   * @return {@link HashMap} mapping gene names to `.pdb` files, one for each parsed gene name.
   */
  public HashMap<String, File> getPdInputFiles() {
    return pdInputFiles;
  }

  /**
   * Prints help information for command line arguments.
   *
   * @param options       The specified command line argument {@link Options}.
   * @param helpFormatter A {@link HelpFormatter} for displaying help information.
   */
  private void printHelp(Options options, HelpFormatter helpFormatter) {
    helpFormatter.printHelp(ArgumentsParser.CLASS_NAME, options);
  }

  /**
   * Populates the {@link ArgumentsParser} sampleInput and sampleNames property with the content from a {@link File}
   * object. The {@link File} object is expected to contain per line: A path to a .vcf file and optional
   * comma-separated a sample name.
   *
   * @param file A {@link File} pointing to a .txt file specifying the sample input.
   */
  private void parseSampleInput(File file) throws FileNotFoundException, MusialIOException {
    List<String> fileContent = IO.getLinesFromFile(file.getAbsolutePath());
    for (String fileLine : fileContent) {
      String[] fileLineSplit = fileLine.split(",");
      int fileLineTags = fileLineSplit.length;
      File sampleVcfFile;
      String sampleName;
      switch (fileLineTags) {
        case 1 -> {
          sampleVcfFile = new File(fileLineSplit[0].trim());
          sampleName = FilenameUtils.removeExtension(sampleVcfFile.getName());
        }
        case 2 -> {
          sampleVcfFile = new File(fileLineSplit[0].trim());
          sampleName = fileLineSplit[1].trim();
        }
        default -> {
          sampleVcfFile = new File("");
          sampleName = "";
        }
      }
      if (Validation.validateInputFile(sampleVcfFile)) {
        throw new MusialIOException(
            "The following input file does not exist or has no read permission:\t" + sampleVcfFile.getAbsolutePath());
      }
      sampleInput.add(sampleVcfFile);
      sampleNames.add(sampleName);
    }
  }

  /**
   * Populates the {@link ArgumentsParser} includedGeneFeatures and pdbInputFiles property with the content from a
   * {@link File} object. The {@link File} object is expected to contain per line: A String to match the feature from
   * the Name attribute in the .gff files entries, (optional) an internal name for the feature, (optional) path to a
   * .pdb file containing protein structure information of the feature.
   *
   * @param file A {@link File} pointing to a .txt file specifying the gene feature input.
   */
  private void parseGeneFeatures(File file)
      throws FileNotFoundException, MusialCLAException, MusialIOException,
      MusialBioException {
    List<String> fileContent = IO.getLinesFromFile(file.getAbsolutePath());
    ArrayList<String> featureQueries = new ArrayList<>();
    ArrayList<String> featureNames = new ArrayList<>();
    HashMap<String, String> featurePDBInformation = new HashMap<>();
    for (String fileLine : fileContent) {
      String[] fileLineSplit = fileLine.split(",");
      int fileLineTags = fileLineSplit.length;
      String featureQuery;
      String featureName;
      String featurePDBPath;
      switch (fileLineTags) {
        case 1 -> {
          featureQuery = fileLineSplit[0].trim();
          featureQueries.add(featureQuery);
          featureNames.add(featureQuery);
        }
        case 2 -> {
          featureQuery = fileLineSplit[0].trim();
          featureName = fileLineSplit[1].trim();
          featureQueries.add(featureQuery);
          featureNames.add(featureName);
        }
        case 3 -> {
          featureQuery = fileLineSplit[0].trim();
          featureName = fileLineSplit[1].trim();
          featurePDBPath = fileLineSplit[2].trim();
          featureQueries.add(featureQuery);
          featureNames.add(featureName);
          featurePDBInformation.put(featureName, featurePDBPath);
          if (Validation.validateInputFile(new File(featurePDBPath))) {
            throw new MusialIOException(
                "The following input file does not exist or has no read permission:\t" +
                    featurePDBPath);
          }
        }
      }
    }
    FeatureList referenceAnnotationFeatures;
    try {
      referenceAnnotationFeatures = IO.readGFF(this.annotationInput);
    } catch (IOException e) {
      throw new MusialCLAException(
          "`-a` The specified reference genome annotation could not be read:\t" + e.getMessage());
    }
    HashSet<String> parsedGeneNames = new HashSet<>();
    for (int i = 0; i < featureQueries.size(); i++) {
      String featureQuery = featureQueries.get(i);
      String featureName = featureNames.get(i);
      String featurePDBPath = "";
      if (featurePDBInformation.containsKey(featureQuery)) {
        featurePDBPath = featurePDBInformation.get(featureQuery);
      } else if (featurePDBInformation.containsKey(featureName)) {
        featurePDBPath = featurePDBInformation.get(featureName);
      }
      FeatureList matchedFeatures = referenceAnnotationFeatures.selectByAttribute("Name", featureQuery);
      if (matchedFeatures.size() == 0) {
        throw new MusialCLAException(
            "`-gf` The following specified gene was not found in the reference genome " +
                "annotation:\t" + featureQuery);
      } else if (matchedFeatures.size() > 1) {
        throw new MusialCLAException(
            "`-gf` More than one feature was identified for the following gene in the reference " +
                "genome annotation:\t" + featureQuery);
      } else if (parsedGeneNames.contains(featureQuery)) {
        throw new MusialCLAException(
            "`-gf` The following gene name was specified more than once:\t" + featureQuery);
      } else {
        FeatureI matchedFeature = matchedFeatures.get(0);
        Location featureCoordinates = matchedFeature.location();
        String featureParentSequence = matchedFeature.seqname();
        /*
        TODO:
        Investigate shift of start position by plus one.The returned value does not match .gff.
        This is currently fixed in the FeatureAnalysisEntry.java class.
         */
        this.includedGeneFeatures.add(new FeatureAnalysisEntry(featureName, featureParentSequence,
            featureCoordinates.getBegin(),
            featureCoordinates.getEnd()));
        this.pdInputFiles.put(featureName, new File(featurePDBPath));
        parsedGeneNames.add(featureQuery);
      }
    }
  }

  /**
   * Sets the input `.vcf` file of a sample to a new {@link File}.
   *
   * @param newFile     The {@link File} the samples input should be set to.
   * @param sampleIndex The index of the sample in the samplesInput property.
   */
  public void setSampleInput(File newFile, int sampleIndex) {
    this.sampleInput.set(sampleIndex, newFile);
  }

  /**
   * @return The original command line arguments as {@link Array} of {@link String}.
   */
  public String[] getArguments() {
    return arguments;
  }
}