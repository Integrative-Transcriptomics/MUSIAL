package cli;

import com.google.gson.Gson;
import com.google.gson.internal.LinkedTreeMap;
import datastructure.FeatureEntry;
import datastructure.SampleEntry;
import datastructure.VariantsDictionary;
import exceptions.MusialBioException;
import exceptions.MusialCLException;
import exceptions.MusialIOException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.Array;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.function.Consumer;

import main.Musial;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.json.simple.JSONObject;
import utility.IO;
import utility.Logging;
import utility.Validation;

/**
 * Parses command line interface arguments for the `MUSIAL generateVDict` module.
 * <p>
 * Used to parse and validate passed command line options or print help information to the user. Parsed arguments are
 * stored to be accessible for other components of the tool.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class CLIParametersUpdateVDict implements CLIParameters {

    /**
     * The unprocessed command line content passed by the user.
     */
    public String[] ARGUMENTS = null;
    /**
     * Minimum (read) coverage in order to call a SNV.
     */
    public Double minCoverage;
    /**
     * Minimum allele frequency (with respect to reads) in order to call a SNV.
     */
    public Double minHomFrequency;
    /**
     * Minimum allele frequency (with respect to reads) in order to call heterozygous SNVs.
     */
    public Double minHetFrequency;
    /**
     * Maximum allele frequency (with respect to reads) in order to call heterozygous SNVs.
     */
    public Double maxHetFrequency;
    /**
     * Minimum quality, i.e. the value given in the QUAL field of .vcf files, in order to call a SNV.
     */
    public Double minQuality;
    /**
     * {@link File} object specifying the reference genome sequence .fasta file.
     */
    public File referenceFASTA;
    /**
     * TODO
     */
    public File referenceGFF;
    /**
     * {@link ConcurrentSkipListMap} of {@link String}, {@link File} pairs specifying one `VCF` file per sample to
     * analyze.
     */
    public final ConcurrentSkipListMap<String, SampleEntry> samples = new ConcurrentSkipListMap<>();
    /**
     * {@link ConcurrentSkipListMap} of {@link String}, {@link File} pairs. Each key specifies one feature on the
     * passed reference `FASTA` to be analyzed. Each assigned value points to an optional `PDB` file per feature and
     * induces the analysis variants allocated to the protein structure; values are expected to be `null` or `""` if
     * variant to protein structure allocation shall not be run.
     */
    public final ConcurrentSkipListMap<String, FeatureEntry> features = new ConcurrentSkipListMap<>();
    /**
     * {@link File} object specifying the output file.
     */
    public File outputFile;
    /**
     * TODO
     */
    private FeatureList referenceFeatures = null;

    /**
     * Constructor of the {@link CLIParametersUpdateVDict} class.
     * <p>
     * Used to parse and validate command line interface input. Parsed arguments are stored - and accessible by other
     * components - via class properties.
     *
     * @param args {@link String} {@link Array} containing the command line arguments.
     * @throws MusialCLException If any error occurs during parsing or validating an
     *                           {@link MusialCLException} is thrown.
     * @throws MusialIOException If any error occurs during input file validation.
     */
    public CLIParametersUpdateVDict(String[] args)
            throws MusialCLException, MusialIOException, MusialBioException, FileNotFoundException {
        // Store original command line arguments.
        this.ARGUMENTS = args;
        // Initialize `Option` object with all parameters.
        Options options = new Options();
        options.addOption(Option.builder("vds")
                .longOpt("vDictSpec")
                .desc("TODO")
                .hasArg()
                .required()
                .build());
        options.addOption("d", "debug", false, "Enables debug mode, i.e. stacktrace is printed when errors occur. [false]");
        // Instantiate a formatter for the help message and a default command line parser.
        HelpFormatter helpformatter = new HelpFormatter();
        CommandLineParser parser = new DefaultParser();
        Consumer<Options> checkHelp = (Options checkOptions) -> {
            try {
                parser.parse(checkOptions, args);
            } catch (ParseException e) {
                printHelp(options, helpformatter);
                System.out.println();
                Logging.logError(e.getMessage());
            }
        };
        // If `-h`/`--help` was specified with other options apply the same behaviour.
        checkHelp.accept(options);
        // Else, command line interface arguments are parsed.
        CommandLine cmd;
        JSONObject parameters;
        try {
            cmd = parser.parse(options, args);
            Musial.DEBUG = cmd.hasOption("d");
            // Parse specified variants dictionary to `JSON` object.
            try (
                    BufferedReader bufferedReader = new BufferedReader(
                            new InputStreamReader(
                                    Files.newInputStream(Path.of(cmd.getOptionValue("vds"))),
                                    StandardCharsets.UTF_8)
                    );
            ) {
                Gson gson = new Gson();
                parameters = gson.fromJson(bufferedReader, JSONObject.class);
            } catch (Exception e) {
                throw new MusialCLException(
                        "Failed to read variants dictionary specification " + cmd.getOptionValue("vds") + ": " + e.getMessage());
            }
            // Parse minCoverage
            if (Validation.isPositiveDouble(String.valueOf(parameters.get("minCoverage")))) {
                this.minCoverage = (Double) parameters.get("minCoverage");
            } else {
                throw new MusialCLException(
                        "Invalid `minCoverage` " + parameters.get("minCoverage") + "; expected positive float.");
            }
            // Parse minFrequency
            if (Validation.isPercentage(String.valueOf(parameters.get("minHomFrequency")))) {
                this.minHomFrequency = (Double) parameters.get("minHomFrequency");
            } else {
                throw new MusialCLException(
                        "Invalid `minHomFrequency` " + parameters.get("minHomFrequency") +
                                "; expected float between 0.0 and 1.0.");
            }
            // Parse minHetFrequency
            if (Validation.isPercentage(String.valueOf(parameters.get("minHetFrequency")))) {
                this.minHetFrequency = (Double) parameters.get("minHetFrequency");
            } else {
                throw new MusialCLException(
                        "Invalid `minHetFrequency` " + parameters.get("minHetFrequency") +
                                "; expected float between 0.0 and 1.0.");
            }
            // Parse maxHetFrequency
            if (Validation.isPercentage(String.valueOf(parameters.get("maxHetFrequency")))) {
                this.maxHetFrequency = (Double) parameters.get("maxHetFrequency");
            } else {
                throw new MusialCLException(
                        "Invalid `maxHetFrequency` " + parameters.get("maxHetFrequency") +
                                "; expected float between 0.0 and 1.0.");
            }
            // Parse minQuality
            if (Validation.isPositiveDouble(String.valueOf(parameters.get("minQuality")))) {
                this.minQuality = (Double) parameters.get("minQuality");
            } else {
                throw new MusialCLException(
                        "Invalid `minQuality` " + parameters.get("minQuality") + "; expected positive float.");
            }
            // Parse No. threads to use
            if (Validation.isPositiveInteger(String.valueOf(Math.round((Double) parameters.get("threads"))))) {
                Musial.THREADS = (int) Math.round((Double) parameters.get("threads"));
            } else {
                throw new MusialCLException(
                        "Invalid `threads` " + parameters.get("threads") + "; expected positive integer.");
            }
            // Parse reference genome (.fasta)
            File r = new File((String) parameters.get("referenceFASTA"));
            if (Validation.isFile(r)) {
                this.referenceFASTA = r;
            } else {
                throw new MusialCLException(
                        "Invalid `referenceFASTA` " + parameters.get("referenceFASTA") + "; failed to read file.");
            }
            // Parse reference annotation (.gff/.gff3)
            File a = new File((String) parameters.get("referenceGFF"));
            if (Validation.isFile(a)) {
                this.referenceGFF = a;
            } else {
                throw new MusialCLException(
                        "Invalid `referenceGFF` " + parameters.get("referenceGFF") + "; failed to read file.");
            }
            // Parse output file
            File o = new File((String) parameters.get("outputFile"));
            if (Validation.isDirectory(new File(o.getParent()))) {
                this.outputFile = o;
            } else {
                throw new MusialCLException(
                        "Invalid `outputFile` " + parameters.get("outputFile") + "; unable to access parent directory.");
            }
            // Parse samples
            //noinspection rawtypes
            LinkedTreeMap samples = (LinkedTreeMap) parameters.get("samples");
            //noinspection rawtypes
            LinkedTreeMap sampleEntry;
            for (Object sampleKey : samples.keySet()) {
                //noinspection rawtypes
                sampleEntry = (LinkedTreeMap) samples.get(sampleKey);
                File sampleVcfFile = new File((String) sampleEntry.get("vcfFile"));
                if (sampleEntry.get("vcfFile") != null && !Validation.isFile(sampleVcfFile)) {
                    throw new MusialCLException("Failed to access specified `vcf` file for sample " + sampleKey);
                }
                //noinspection unchecked
                addSample(
                        (String) sampleKey,
                        sampleVcfFile,
                        (Map<String, String>) sampleEntry.get("annotations")
                );
            }
            // Parse features
            //noinspection rawtypes
            LinkedTreeMap features = (LinkedTreeMap) parameters.get("features");
            //noinspection rawtypes
            LinkedTreeMap featureEntry;
            for (Object featureKey : features.keySet()) {
                //noinspection rawtypes
                featureEntry = (LinkedTreeMap) features.get(featureKey);
                File featurePdbFile;
                if ( featureEntry.containsKey("pdbFile") ) {
                    featurePdbFile = new File((String) featureEntry.get("pdbFile"));
                    if (featureEntry.get("pdbFile") != null && !Validation.isFile(featurePdbFile)) {
                        throw new MusialCLException("Failed to access specified `pdb` file for feature " + featureKey + " " + featurePdbFile);
                    }
                } else {
                    featurePdbFile = null;
                }
                //noinspection unchecked
                addFeature(
                        (String) featureKey,
                        (String) featureEntry.get("geneName"),
                        featurePdbFile,
                        (Map<String, String>) featureEntry.get("annotations")
                );
            }
        } catch (ParseException e) {
            Logging.logError("An error occurred during command line interface parameter parsing: " + e.getMessage());
        }
    }

    /**
     * Populates the {@link CLIParametersUpdateVDict} sampleInput and sampleNames property.
     * <p>
     * If the passed {@link File} object points to a file, the properties are populated using the content from this file.
     * The {@link File} object is expected to contain per line: A path to a .vcf file and optional
     * comma-separated a sample name.
     * <p>
     * If the passed {@link File} object points to a directory, the directory and all its sub-directories are searched
     * for .vcf files which are used as sample input.
     *
     * @param file A {@link File} pointing to a .txt file specifying the sample input or a directory.
     * @throws MusialIOException If the passed file does not exist, has no content or no read permission.
     */
    private void addSample(String name, File file, Map<String, String> annotations) throws MusialIOException {
        if (file.isFile()) {
            SampleEntry sampleEntry = new SampleEntry(file, name);
            for (Map.Entry<String, String> annotation : annotations.entrySet()) {
                sampleEntry.annotations.put(annotation.getKey(), annotation.getValue());
            }
            this.samples.put(name, sampleEntry);
        } else {
            throw new MusialIOException(
                    "Specified sample input file does not exist or has no read permission:\t" + file.getAbsolutePath());
        }
    }

    /**
     * TODO
     *
     * @param name
     * @param geneName
     * @param file
     * @param annotations
     * @throws MusialIOException
     * @throws MusialCLException
     * @throws MusialBioException
     */
    private void addFeature(String name, String geneName, File file, Map<String, String> annotations)
            throws MusialIOException, MusialCLException, MusialBioException {
        if (this.referenceFeatures == null) {
            try {
                this.referenceFeatures = IO.readGFF(this.referenceGFF);
            } catch (IOException e) {
                throw new MusialIOException(
                        "Failed to read specified reference genome annotation " + this.referenceGFF.getAbsolutePath() + "\t" +
                                e.getMessage());
            }
        }
        FeatureList matchedFeatures = this.referenceFeatures.selectByAttribute("Name", geneName);
        if (matchedFeatures.size() == 0) {
            throw new MusialCLException("Feature " + geneName + " not found in the specified reference genome annotation.");
        } else if (matchedFeatures.size() > 1) {
            throw new MusialCLException("Feature " + geneName + " found more than once in the specified reference genome " +
                    "annotation.");
        } else if (this.features.containsKey(geneName)) {
            throw new MusialCLException("Feature " + geneName + " was specified multiple times.");
        } else {
            FeatureI matchedFeature = matchedFeatures.get(0);
            Location featureCoordinates = matchedFeature.location();
            String featureParentSequence = matchedFeature.seqname();
      /* FIXME: Starting positions are shifted by minus one, i.e. the returned values do not match with the ones of
          the `gff` files. This is currently fixed in the FeatureEntry.java class.
       */
            FeatureEntry featureEntry = new FeatureEntry(name, featureParentSequence, featureCoordinates.getBegin(),
                    featureCoordinates.getEnd());
            featureEntry.annotations.put("geneName", geneName);
            for (Map.Entry<String, String> annotationEntry : annotations.entrySet()) {
                featureEntry.annotations.put(annotationEntry.getKey(), annotationEntry.getValue());
            }
            this.features.put(name, featureEntry);
        }
        if (file != null && file.exists()) {
            if (Validation.isFile(file)) {
                this.features.get(name).pdbFile = file;
            } else {
                throw new MusialIOException(
                        "Specified feature `PDB` file " + file.getAbsolutePath() + " does not exist or has no read permission.");
            }
        }
    }

}