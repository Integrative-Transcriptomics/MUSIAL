package tools;

import datastructure.GeneFeature;
import exceptions.MusialCLAException;
import exceptions.MusialDuplicationException;
import exceptions.MusialFaultyDataException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
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
 * <p>
 * TODO: The following options are deprecated since version 2.0:
 * options.addOption("u", "uncovered", false, "write uncovered positions for each sample");
 * options.addOption("lc", "lowCov", false, "calculate low coverage regions");
 * options.addOption("ca", "covAdd", true, "minCoverage to make additional call ["+this.minCovAdd+"]");
 * options.addOption("d", "compare", true, "calculate differences compared to other SNVTable");
 * options.addOption("s", "sharedAlleleFreq", false, "write shared allele frequencies for each position");
 * options.addOption("ogf", "outgroupFile", true, "vcfFile for outgroup");
 * options.addOption("ogn", "outgroupName", true, "name for outgroup (contained in input)");
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
   * Number of threads to use for computations.
   */
  private int numThreads = 1;
  /**
   * Whether temporary files should be deleted or not.
   */
  private boolean clean = false;
  /**
   * Minimum (read) coverage in order to call a SNV.
   */
  private Double minCoverage = 5d;
  /**
   * Minimum allele frequency in order to call a SNV.
   */
  private Double minFrequency = 0.9;
  /**
   * Minimum (mapping?) quality in order to call a SNV.
   */
  private Double minQuality = 30d;
  /**
   * {@link File} object representing the reference genome sequence .fasta file.
   */
  private final File referenceInput;
  /**
   * Whether the reference genome sequence should be added to the alignment.
   */
  private boolean addReference = false;
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
   * {@link ArrayList} of {@link Integer} objects representing reference genome positions that should not be
   * considered for SNV calling.
   */
  private final ArrayList<Integer> excludedPositions = new ArrayList<>();
  /**
   * Whether a phylogenetic tree should be computed considering the samples as taxa.
   */
  private boolean computePhylogenetics = false;
  /**
   * Whether SNPEff is called on each input samples .vcf file before further analysis.
   */
  private boolean runSnpEff = false;
  /**
   * {@link ArrayList} of {@link GeneFeature} objects representing gene names, their parent sequence and locations,
   * that are exclusively analysed instead of the whole genome.
   */
  private final ArrayList<GeneFeature> includedGeneFeatures = new ArrayList<>();
  /**
   * TODO: Write doc. once finally implemented.
   */
  private Double minHet = 0.45;
  /**
   * TODO: Write doc. once finally implemented.
   */
  private Double maxHet = 0.55;
  /**
   * Whether heterozygous SNV calls should be considered.
   */
  private Boolean callHeterozygous = false;
  /**
   * Whether config files for the interactive visualization extension should be generated.
   */
  private Boolean generateIveConfigs = false;
  /**
   * {@link ArrayList} of {@link File} objects pointing to .pdb files, one for each gene name.
   */
  private final ArrayList<File> pdInputFiles = new ArrayList<>();
  /**
   * Whether information parsed from `.pdb` files should be included in the analysis.
   */
  private boolean includeStructureInformation = false;

  /**
   * Constructor of the ArgumentParser class.
   * <p>
   * The ArgumentParser is used to parse and validate the user command line input for Musial. Parsed arguments are
   * stored - and accessible - via class properties.
   *
   * @param args {@link String} {@link Array} containing the command line arguments.
   * @throws MusialCLAException If any error occurs during parsing or validating an
   *                            {@link MusialCLAException} is thrown.
   */
  public ArgumentsParser(String[] args)
      throws MusialCLAException, ParseException, MusialFaultyDataException, MusialDuplicationException {
    // Add option to print help message.
    Options helpOptions = new Options();
    helpOptions.addOption("h", "help", false, "Display help information.");
    // Add non-help options.
    Options options = new Options();
    options.addOption("h", "help", false, "Display help information.");
    // Add options used to specify the filtering of SNVs.
    options.addOption("t", "threads", true, "Number of threads to use [" + this.numThreads + "]");
    options.addOption("c", "coverage", true, "Minimum coverage to call a SNV [" + this.minCoverage + "]");
    options.addOption("f", "frequency", true, "Minimum frequency to call a SNV [" + this.minFrequency + "]");
    options.addOption("q", "minQual", true, "Minimum quality to call a SNV [" + this.minQuality + "]");
    options.addOption("x", "clean", false, "Whether temporary files should be deleted [" + this.clean +
        "]");
    // Add options used to specify input data.
    options.addOption(Option.builder("r")
        .longOpt("reference")
        .desc("Reference genome sequence in .fasta format.")
        .required()
        .hasArg()
        .build());
    options.addOption("ar", "addReference", false,
        "Add reference sequence information to the output. [" + this.addReference + "]");
    options.addOption(Option.builder("il")
        .longOpt("inputList")
        .desc("Input files (i.e. filepaths) in .vcf format as list separated by ','.")
        .hasArg()
        .hasArgs()
        .valueSeparator(',')
        .build());
    options.addOption(Option.builder("if")
        .longOpt("inputFile")
        .desc("Input files in .vcf format from a file, each line has to contain one filepath.")
        .hasArg()
        .hasArgs()
        .build());
    options.addOption(Option.builder("ie")
        .longOpt("inputEager")
        .desc("Output directory of an EAGER run used to detect input .vcf files.")
        .hasArg()
        .build());
    options.addOption(Option.builder("sn")
        .longOpt("sampleNames")
        .desc("List of strings used as sample names separated by ','. If none are given the input file names are used" +
            " as sample names.")
        .hasArg()
        .hasArgs()
        .valueSeparator(',')
        .build());
    options.addOption("g", "gff", true, "Genome annotation file in .gff format.");
    options.addOption(Option.builder("o")
        .longOpt("output")
        .desc("the output folder")
        .required()
        .hasArg()
        .build());
    options.addOption(Option.builder("e")
        .longOpt("exclude")
        .desc("Exclude given positions in alignment. Expects a list of integers separated by ','.")
        .hasArg()
        .hasArgs()
        .valueSeparator(',')
        .build());
    // Add options to specify optional behaviour/additional computations of MUSIAL.
    options.addOption("p", "phylogenetics", false, "Calculate a phylogenetic tree. [" + this.computePhylogenetics +
        "]");
    options.addOption("sf", "snpEff", false, "Run SNPEff. [" + this.runSnpEff + "]");
    options.addOption(Option.builder("gn")
        .longOpt("geneNames")
        .desc("List of strings, separated by ',', of gene names to analyze exclusively instead of the whole genome. " +
            "Requires the specification of a .gff reference genome sequence annotation and each passed gene should be" +
            " identifiable via this .gff file.")
        .hasArg()
        .hasArgs()
        .valueSeparator(',')
        .build());
    options.addOption(Option.builder("gf")
        .longOpt("geneFile")
        .desc("Gene names to analyze exclusively instead of the whole genome parsed from a file, each line has to " +
            "contain one gene name. Requires the specification of a .gff reference genome sequence annotation and " +
            "each passed gene should be identifiable via this .gff file.")
        .hasArg()
        .hasArgs()
        .build());
    options.addOption(Option.builder("pf")
        .longOpt("proteinFile")
        .desc("Paths to protein data `.pdb` files parsed from a file. Each line has to contain one path and the " +
            "number of `.pdb` files is expected to match any passed number of genes to analyze. In advance the first " +
            "file should represent the protein structure of the first gene passed.")
        .hasArg()
        .hasArgs()
        .build());
    options.addOption("ive", "generateIveConfig", false,
        "Generate config files for the interactive visualization extension. [" + this.generateIveConfigs +
            "]");
    options.addOption(Option.builder("y")
        .longOpt("heterozygous")
        .desc("Call heterozygous positions, if set. Expected format is min,max. [" + this.minHet + "," + this.maxHet +
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
    // Next it is tried to parse and validate command line options.
    cmd = parser.parse(options, args);
    // In the following each possible command line argument is validated, if any faulty input is detected the
    // application will exit with an exception.
    // Validation of SNV filtering parameters:
    if (cmd.hasOption("t")) {
      if (Validation.isPositiveInteger(cmd.getOptionValue("t"))) {
        this.numThreads = Integer.parseInt(cmd.getOptionValue("t"));
      } else {
        throw new MusialCLAException("`-t` Number of threads is no positive integer:\t" + cmd.getOptionValue("t"));
      }
    }
    /*
    if (cmd.hasOption("m")) {
      if (Validation.isPositiveInteger(cmd.getOptionValue("m"))) {
        this.memory = Integer.parseInt(cmd.getOptionValue("m"));
      } else {
        throw new MusialCLAException("`-m` Maximal memory is no positive integer:\t" + cmd.getOptionValue(
            "m"));
      }
    }
    */
    if (cmd.hasOption("c")) {
      if (Validation.isPositiveInteger(cmd.getOptionValue("c"))) {
        this.minCoverage = (double) Integer.parseInt(cmd.getOptionValue("c"));
      } else {
        throw new MusialCLAException("`-c` Value for min. coverage is no positive integer:\t" + cmd.getOptionValue(
            "c"));
      }
    }
    if (cmd.hasOption("f")) {
      if (Validation.isPercentage(cmd.getOptionValue("f"))) {
        this.minFrequency = Double.parseDouble(cmd.getOptionValue("f"));
      } else {
        throw new MusialCLAException(
            "`-f` Value for min. frequency is no double in the interval [0.0,1.0]:\t" + cmd.getOptionValue(
                "f"));
      }
    }
    if (cmd.hasOption("q")) {
      if (Validation.isPositiveDouble(cmd.getOptionValue("q"))) {
        this.minQuality = Double.parseDouble(cmd.getOptionValue("q"));
      } else {
        throw new MusialCLAException(
            "`-q` Value for min. quality is no positive double:\t" + cmd.getOptionValue("q"));
      }
    }
    if (cmd.hasOption("x")) {
      this.clean = true;
    }
    // Validation of data input parameters:
    this.referenceInput = new File(cmd.getOptionValue('r'));
    if (!Validation.validateInputFile(this.referenceInput)) {
      throw new MusialCLAException("`-r` The specified input file does not exist or has no read permission:\t" +
          this.referenceInput.getAbsolutePath());
    }
    if (cmd.hasOption("ar")) {
      this.addReference = true;
    }
    // Assert that exactly one input-option is specified.
    String inputType = "";
    if (!(cmd.hasOption("il") || cmd.hasOption("if") || cmd.hasOption("ie")) ||
        (cmd.hasOption("il") && cmd.hasOption("if")) ||
        (cmd.hasOption("il") && cmd.hasOption("ie")) ||
        (cmd.hasOption("if") && cmd.hasOption("ie"))) {
      throw new MusialCLAException("Exactly one input option of il, if or ie has to be specified.");
    } else {
      if (cmd.hasOption("il")) {
        inputType = "il";
        populateSampleInput(cmd.getOptionValues("il"));
      } else if (cmd.hasOption("if")) {
        inputType = "if";
        try {
          populateSampleInput(IO.getLinesFromFile(cmd.getOptionValue("if")));
        } catch (FileNotFoundException e) {
          throw new MusialCLAException(
              "`-" + inputType + "` The specified input file does not exist:\t" + cmd.getOptionValue(
                  "if"));
        }
      } else if (cmd.hasOption("ie")) {
        inputType = "ie";
        try {
          ArrayList<ArrayList<String>> eagerOutputInformation = IO.getInputFromEagerOutput(cmd.getOptionValue("ie"));
          populateSampleInput(eagerOutputInformation.get(0));
          populateSampleNames(eagerOutputInformation.get(1));
        } catch (Exception e) {
          throw new MusialCLAException(
              "`-" + inputType + "` During collection of input information from the specified EAGER output (" +
                  cmd.getOptionValue("ie") + ") an error occurred:\t" + e.getMessage());
        }
      }
    }
    if (this.sampleInput.size() == 0) {
      throw new MusialCLAException("`-" + inputType + "` No input .vcf files were found.");
    }
    for (File file : this.sampleInput) {
      if (!Validation.validateInputFile(file)) {
        throw new MusialCLAException("`-" + inputType + "` The specified input file does not exist or has no read " +
            "permission:\t" +
            file.getAbsolutePath());
      }
    }
    if (this.sampleNames.size() == 0) {
      if (cmd.hasOption("sn")) {
        populateSampleNames(cmd.getOptionValues("sn"));
      } else {
        populateSampleNames(this.sampleInput.stream().map(file -> FilenameUtils.removeExtension(file.getName()))
            .collect(Collectors.toList()));
      }
    }
    if (this.sampleInput.size() != this.sampleNames.size()) {
      throw new MusialCLAException("The specified sample input files and names do not match in their size.");
    }
    if (cmd.hasOption("g")) {
      this.annotationInput = new File(cmd.getOptionValue("g"));
      if (!Validation.validateInputFile(this.annotationInput)) {
        throw new MusialCLAException("`-g` The specified input file does not exist or has no read permission:\t" +
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
    if (cmd.hasOption("e")) {
      for (String optionValue : cmd.getOptionValues("e")) {
        if (!Validation.isPositiveInteger(optionValue)) {
          throw new MusialCLAException(
              "`-e` The specified excluded positions contain non-integer values:\t" + optionValue);
        } else {
          this.excludedPositions.add(Integer.parseInt(optionValue));
        }
      }
    }
    // Validate optional behaviour options:
    if (cmd.hasOption("p")) {
      this.computePhylogenetics = true;
    }
    if (cmd.hasOption("sf")) {
      if (!cmd.hasOption("g")) {
        throw new MusialCLAException(
            "`-sf` For the SNPEff analysis you need to input a reference annotation (see 'g' option)" +
                ".");
      } else {
        this.runSnpEff = true;
      }
    }
    boolean hasGeneNames = false;
    String geneInputType = "";
    List<String> geneNames = new ArrayList<>();
    if (cmd.hasOption("gn")) {
      if (!cmd.hasOption("g")) {
        throw new MusialCLAException(
            "`-gn` For name based gene alignment analysis you need to input a reference annotation " +
                "(see 'g' option).");
      }
      geneNames = Arrays.asList(cmd.getOptionValues("gn"));
      hasGeneNames = true;
      geneInputType = "gn";
    }
    if (cmd.hasOption("gf")) {
      if (cmd.hasOption("gn")) {
        throw new MusialCLAException(
            "`-gf` Included gene names and loci were already processed from the 'gn' option, please " +
                "specify only one of 'gn' and 'gf'.");
      }
      if (!cmd.hasOption("g")) {
        throw new MusialCLAException(
            "`-gf` For name based gene alignment analysis you need to input a reference annotation " +
                "(see 'g' option).");
      }
      try {
        geneNames = IO.getLinesFromFile(cmd.getOptionValue("gf"));
        hasGeneNames = true;
        geneInputType = "gf";
      } catch (FileNotFoundException e) {
        throw new MusialCLAException(
            "`-gf` The specified input file does not exist:\t" + cmd.getOptionValue("gf"));
      }
    }
    if (hasGeneNames) {
      FeatureList referenceAnnotationFeatures;
      try {
        referenceAnnotationFeatures = IO.readGFF(this.annotationInput);
      } catch (IOException e) {
        throw new MusialCLAException(
            "`-" + geneInputType + "` The specified reference genome annotation could not be read:\t" + e.getMessage());
      }
      HashSet<String> parsedGeneNames = new HashSet<>();
      for (String geneName : geneNames) {
        FeatureList geneFeatures = referenceAnnotationFeatures.selectByAttribute("Name", geneName);
        if (geneFeatures.size() == 0) {
          throw new MusialCLAException(
              "`-" + geneInputType + "` The following specified gene was not found in the reference genome " +
                  "annotation:\t" + geneName);
        } else if (geneFeatures.size() > 1) {
          throw new MusialCLAException(
              "`-" + geneInputType + "` More than one feature was identified for the following gene in the reference " +
                  "genome annotation:\t" + geneName);
        } else if (parsedGeneNames.contains(geneName)) {
          throw new MusialDuplicationException(
              "`-" + geneInputType + "` The following gene name was specified more than once:\t" + geneName);
        } else {
          FeatureI geneFeature = geneFeatures.get(0);
          Location geneLocation = geneFeature.location();
          String geneSeqName = geneFeature.seqname();
          this.includedGeneFeatures.add(new GeneFeature(geneSeqName, geneName, geneLocation.start(),
              geneLocation.end()));
          parsedGeneNames.add(geneName);
        }
      }
    }
    if ( cmd.hasOption("ive") ) {
      this.generateIveConfigs = true;
    }
    if ( cmd.hasOption("pf") ) {
      List<String> proteinFilePaths;
      try {
        proteinFilePaths = IO.getLinesFromFile(cmd.getOptionValue("pf"));
      } catch (FileNotFoundException e) {
        throw new MusialCLAException("`-pf` The specified input file does not exist:\t" + cmd.getOptionValue("pf"));
      }
      if ( proteinFilePaths.size() != geneNames.size() ) {
        throw new MusialCLAException("`-pf` The number of specified protein files (" + proteinFilePaths.size() + ") " +
            "and gene names (" + geneNames.size() + ") do not match.");
      }
      for (String proteinFilePath : proteinFilePaths) {
        this.pdInputFiles.add(new File(proteinFilePath));
      }
      this.includeStructureInformation = true;
    }
    if (cmd.hasOption("y")) {
      this.callHeterozygous = true;
      String[] percentages = cmd.getOptionValues("y");
      if (percentages != null) {
        if (percentages.length == 2) {
          if (Validation.isPercentage(percentages[0]) && Validation.isPercentage(percentages[1])) {
            double min = Double.parseDouble(percentages[0]);
            double max = Double.parseDouble(percentages[1]);
            if (min > max) {
              throw new MusialCLAException(
                  "`-y` The specified minimum percentage is larger than the specified maximum percentage:\t" +
                      cmd.getOptionValue("y"));
            }
            this.minHet = min;
            this.maxHet = max;
          } else {
            throw new MusialCLAException(
                "`-y` One of the specified values is no percentage:\t" + cmd.getOptionValue("y"));
          }
        } else {
          throw new MusialCLAException(
              "`-y` Wrong number of arguments, expected two:\t" + cmd.getOptionValue("y"));
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
   * @return Whether reference information should be added.
   */
  public boolean isAddReference() {
    return addReference;
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
  public File getReferenceInput() {
    return referenceInput;
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
   * @return The specified excluded positions as {@link ArrayList<Integer>}.
   */
  public ArrayList<Integer> getExcludedPositions() {
    return excludedPositions;
  }

  /**
   * @return Whether a phylogeny should be computed.
   */
  public boolean isComputePhylogenetics() {
    return computePhylogenetics;
  }

  /**
   * @return Whether SnpEff should be run.
   */
  public boolean isRunSnpEff() {
    return runSnpEff;
  }

  /**
   * @return The included gene features as {@link ArrayList<GeneFeature>}.
   */
  public ArrayList<GeneFeature> getIncludedGeneFeatures() {
    return includedGeneFeatures;
  }

  /**
   * @return ?
   */
  public Double getMinHet() {
    return minHet;
  }

  /**
   * @return ?
   */
  public Double getMaxHet() {
    return maxHet;
  }

  /**
   * @return Whether heterozygous calls should be made.
   */
  public Boolean getCallHeterozygous() {
    return callHeterozygous;
  }

  /**
   * @return Whether temporary files should be deleted after running the tool.
   */
  public boolean isClean() {
    return clean;
  }

  /**
   * @return Whether config files for the interactive visualization extension should be generated.
   */
  public Boolean isGenerateIveConfigs() {
    return generateIveConfigs;
  }

  /**
   * @return Whether protein structure information is passed and can be included into the analysis.
   */
  public boolean isIncludeStructureInformation() {
    return includeStructureInformation;
  }

  /**
   * @return {@link ArrayList<File>} pointing to `.pdb` files, one for each parsed gene name.
   */
  public ArrayList<File> getPdInputFiles() {
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
   * Populates the {@link ArgumentsParser::sampleInput} property with {@link File} objects pointing to input .vcf files.
   *
   * @param li A {@link String[]} containing paths to .vcf files.
   */
  private void populateSampleInput(String[] li) {
    for (String inputFile : li) {
      this.sampleInput.add(new File(inputFile));
    }
  }

  /**
   * Populates the {@link ArgumentsParser::sampleInput} property with {@link File} objects pointing to input .vcf files.
   *
   * @param li A {@link List<String>} containing paths to .vcf files.
   */
  private void populateSampleInput(List<String> li) {
    for (String inputFile : li) {
      this.sampleInput.add(new File(inputFile));
    }
  }

  /**
   * Populates the {@link ArgumentsParser::sampleNames} property with {@link String} objects specifying sample names.
   *
   * @param li A {@link String[]} containing sample names.
   */
  private void populateSampleNames(List<String> li) {
    this.sampleNames.addAll(li);
  }

  /**
   * Populates the {@link ArgumentsParser::sampleNames} property with {@link String} objects specifying sample names.
   *
   * @param li A {@link List<String>} containing sample names.
   */
  private void populateSampleNames(String[] li) {
    this.sampleNames.addAll(Arrays.asList(li));
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
}