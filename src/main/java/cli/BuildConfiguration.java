package cli;

import com.google.gson.internal.LinkedTreeMap;
import datastructure.Feature;
import datastructure.FeatureCoding;
import datastructure.Sample;
import exceptions.MusialException;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.reference.BlockCompressedIndexedFastaSequenceFile;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import main.Musial;
import org.apache.commons.io.FilenameUtils;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.genome.parsers.gff.Location;
import utility.Bio;
import utility.IO;
import utility.Logger;
import utility.Validation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Parse parameters for the `build` task from a {@link HashMap}.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.1
 */
public final class BuildConfiguration {
    /**
     * Minimum (read) coverage in order to call an SNV.
     */
    public Double minimalCoverage;
    /**
     * Minimum allele frequency (with respect to reads) in order to call an SNV.
     */
    public Double minimalHomozygousFrequency;
    /**
     * Minimum allele frequency (with respect to reads) in order to call heterozygous SNVs.
     */
    public Double minimalHeterozygousFrequency;
    /**
     * Maximum allele frequency (with respect to reads) in order to call heterozygous SNVs.
     */
    public Double maximalHeterozygousFrequency;
    /**
     * Minimum quality, i.e. the value given in the QUAL field of .vcf files, in order to call an SNV.
     */
    public Double minimalQuality;
    /**
     * {@link File} object pointing to the reference genome sequence file.
     */
    public File referenceSequenceFile;
    /**
     * {@link BlockCompressedIndexedFastaSequenceFile} object providing the compressed indexed reference genome sequence.
     */
    public IndexedFastaSequenceFile referenceSequence;
    /**
     * {@link File} object pointing to the reference features file.
     */
    public File referenceFeaturesFile;
    /**
     * {@link FeatureList} object specifying the reference genome features.
     */
    public FeatureList referenceFeatures;
    /**
     * {@link ConcurrentSkipListMap} of {@link String}/{@link Sample} pairs.
     */
    public final ConcurrentSkipListMap<String, Sample> samples = new ConcurrentSkipListMap<>();
    /**
     * {@link ConcurrentSkipListMap} of {@link String}/{@link Feature} pairs.
     */
    public final ConcurrentSkipListMap<String, Feature> features = new ConcurrentSkipListMap<>();
    /**
     * {@link File} object specifying the output directory.
     */
    public File output;
    /**
     * Optional positions to exclude from the analysis.
     */
    public final HashMap<String, TreeSet<Integer>> excludedPositions = new HashMap<>();
    /**
     * Prefix to use for exception messages.
     */
    @SuppressWarnings("FieldCanBeLocal")
    private final String EXCEPTION_PREFIX = "(Task `build` Argument Validation)";

    /**
     * Constructor of the {@link BuildConfiguration} class.
     * <p>
     * Used to parse and validate command line interface input. Parsed arguments are stored - and accessible by other
     * components - via class properties.
     *
     * @param parameters {@link LinkedTreeMap} containing parameter information.
     * @throws MusialException If IO file validation fails; If CLI parameter validation fails.
     */
    public BuildConfiguration(HashMap<String, Object> parameters)
            throws MusialException {
        // Parse minimalCoverage
        if (parameters.containsKey("minimalCoverage")
                && Validation.isPositiveDouble(String.valueOf(parameters.get("minimalCoverage")))) {
            this.minimalCoverage = (Double) parameters.get("minimalCoverage");
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid or missing parameter `minimalCoverage`; expected positive float.");
        }

        // Parse minimalHomozygousFrequency
        if (parameters.containsKey("minimalHomozygousFrequency")
                && Validation.isPercentage(String.valueOf(parameters.get("minimalHomozygousFrequency")))) {
            this.minimalHomozygousFrequency = (Double) parameters.get("minimalHomozygousFrequency");
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid or missing parameter `minimalHomozygousFrequency`; expected float between 0.0 and 1.0.");
        }

        // Parse minimalHeterozygousFrequency
        if (parameters.containsKey("minimalHeterozygousFrequency")
                && Validation.isPercentage(String.valueOf(parameters.get("minimalHeterozygousFrequency")))) {
            this.minimalHeterozygousFrequency = (Double) parameters.get("minimalHeterozygousFrequency");
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid or missing parameter `minimalHeterozygousFrequency`; expected float between 0.0 and 1.0.");
        }

        // Parse maximalHeterozygousFrequency
        if (parameters.containsKey("maximalHeterozygousFrequency")
                && Validation.isPercentage(String.valueOf(parameters.get("maximalHeterozygousFrequency")))) {
            this.maximalHeterozygousFrequency = (Double) parameters.get("maximalHeterozygousFrequency");
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid or missing parameter `maximalHeterozygousFrequency`; expected float between 0.0 and 1.0.");
        }

        // Parse minimalQuality
        if (parameters.containsKey("minimalQuality")
                && Validation.isPositiveDouble(String.valueOf(parameters.get("minimalQuality")))) {
            this.minimalQuality = (Double) parameters.get("minimalQuality");
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid or missing parameter `minimalQuality`; expected positive float.");
        }

        // Parse reference sequence (.fasta)
        if (!parameters.containsKey("referenceSequenceFile")) {
            throw new MusialException(EXCEPTION_PREFIX + " Missing parameter `referenceSequenceFile`; expected path to .fasta file.");
        }
        this.referenceSequenceFile = new File((String) parameters.get("referenceSequenceFile"));
        if (Validation.isFile(referenceSequenceFile)) {
            try {
                try {
                    FastaSequenceIndexCreator.create(referenceSequenceFile.toPath(), false);
                } catch (SAMException e) {
                    // Exception is raised if file already exists. This can be ignored.
                }
                this.referenceSequence = new IndexedFastaSequenceFile(referenceSequenceFile.toPath());
            } catch (IOException e) {
                throw new MusialException(EXCEPTION_PREFIX + " Failed to create or locate index for reference sequence file '" + referenceSequenceFile.getAbsolutePath() + "'; " + e.getMessage());
            }
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + referenceSequenceFile.getAbsolutePath() + "' for parameter `referenceSequenceFile`; failed to access file.");
        }

        // Parse reference features (.gff/.gff3)
        if (!parameters.containsKey("referenceFeaturesFile")) {
            throw new MusialException(EXCEPTION_PREFIX + " Missing parameter `referenceFeaturesFile`; expected path to .gff file.");
        }
        this.referenceFeaturesFile = new File((String) parameters.get("referenceFeaturesFile"));
        if (Validation.isFile(referenceFeaturesFile)) {
            try {
                System.setOut(Musial.EMPTY_STREAM);
                this.referenceFeatures = GFF3Reader.read(referenceFeaturesFile.getCanonicalPath());

            } catch (IOException e) {
                throw new MusialException(EXCEPTION_PREFIX + " Failed to parse reference features from file '" + referenceFeaturesFile.getAbsolutePath() + "'; " + e.getMessage());
            } finally {
                System.setOut(Musial.ORIGINAL_OUT_STREAM);
            }
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + referenceFeaturesFile.getAbsolutePath() + "' for parameter `referenceFeatures`; failed to read file.");
        }

        // Parse output file
        if (!parameters.containsKey("output")) {
            throw new MusialException(EXCEPTION_PREFIX + " Missing parameter `outputFile`; expected path to file.");
        }
        File outputFile = new File((String) parameters.get("output"));
        if (!Validation.isDirectory(new File(outputFile.getParent()))) {
            IO.generateDirectory(outputFile);
        }
        if (Validation.isFile(outputFile)) {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid value '" + outputFile.getAbsolutePath() + "' for parameter `output`; specified file already exists.");
        }
        this.output = outputFile;

        // Parse excluded positions
        if (parameters.containsKey("excludedPositions")) {
            //noinspection rawtypes,unchecked
            LinkedTreeMap<String, ArrayList<Double>> excludedPositions = (LinkedTreeMap) parameters.get("excludedPositions");
            for (String key : excludedPositions.keySet()) {
                this.excludedPositions.put(key, new TreeSet<>());
                for (Double excludedPosition : excludedPositions.get(key)) {
                    this.excludedPositions.get(key).add(excludedPosition.intValue());
                }
            }
        }

        // Parse sample information
        Logger.logStatus("Parse sample information");
        @SuppressWarnings("rawtypes")
        LinkedTreeMap sampleEntry;
        // Check if any samples can be parsed from a passed samplesDirectory key in the configuration.
        if (parameters.containsKey("samplesDirectory")) {
            HashSet<String> sampleFilesFromDirectory;
            try (Stream<Path> stream = Files.list(Paths.get((String) parameters.get("samplesDirectory")))) {
                sampleFilesFromDirectory = (HashSet<String>) stream
                        .filter(file -> !Files.isDirectory(file))
                        .map(Path::getFileName)
                        .map(Path::toString)
                        .collect(Collectors.toSet());
            } catch (IOException e) {
                throw new MusialException(EXCEPTION_PREFIX + " Unable to access directory '" + parameters.get("samplesDirectory") + "' specified for parameter `samplesDirectory`");
            }
            for (String sampleFile : sampleFilesFromDirectory) {
                if (sampleFile.endsWith(".vcf")) {
                    String sampleName = FilenameUtils.removeExtension(sampleFile);
                    File sampleVcfFile = new File(parameters.get("samplesDirectory") + sampleFile);
                    addSample(
                            sampleName,
                            sampleVcfFile,
                            new HashMap<>()
                    );
                }
            }
        }

        // Check if any samples can be parsed explicitly from the samples key in the configuration.
        if (parameters.containsKey("samples")) {
            @SuppressWarnings("rawtypes")
            LinkedTreeMap samples = (LinkedTreeMap) parameters.get("samples");
            for (Object sampleKey : samples.keySet()) {
                //noinspection rawtypes
                sampleEntry = (LinkedTreeMap) samples.get(sampleKey);
                File sampleVcfFile = new File((String) sampleEntry.get("vcfFile"));
                if (sampleEntry.get("vcfFile") != null && !Validation.isFile(sampleVcfFile)) {
                    throw new MusialException(EXCEPTION_PREFIX + " Failed to access .vcf file for sample " + sampleKey);
                }
                Map<String, String> annotations;
                if (sampleEntry.containsKey("annotations")) {
                    //noinspection unchecked
                    annotations = (Map<String, String>) sampleEntry.get("annotations");
                } else {
                    annotations = new HashMap<>();
                }
                addSample(
                        (String) sampleKey,
                        sampleVcfFile,
                        annotations
                );
            }
        }
        if ((!parameters.containsKey("samples") && !parameters.containsKey("samplesDirectory")) || this.samples.keySet().size() == 0) {
            throw new MusialException(EXCEPTION_PREFIX + " Expected at least one sample entry being parsed from `samples` and/or `samplesDirectory` parameter");
        }

        // Parse features
        Logger.logStatus("Parse feature information");
        if (!parameters.containsKey("features")) {
            throw new MusialException(EXCEPTION_PREFIX + " Missing parameter `features`; expected at least one feature entry.");
        }
        //noinspection rawtypes
        LinkedTreeMap features = (LinkedTreeMap) parameters.get("features");
        //noinspection rawtypes
        LinkedTreeMap featureEntry;
        for (Object featureKey : features.keySet()) {
            //noinspection rawtypes
            featureEntry = (LinkedTreeMap) features.get(featureKey);
            File optionalPdbFile;
            boolean cds = false;
            if (featureEntry.containsKey("pdbFile")) {
                optionalPdbFile = new File((String) featureEntry.get("pdbFile"));
                if (featureEntry.get("pdbFile") != null && !Validation.isFile(optionalPdbFile)) {
                    throw new MusialException(EXCEPTION_PREFIX + " Failed to access specified .pdb file for feature " + featureKey);
                }
                cds = true;
            } else {
                optionalPdbFile = null;
                if (featureEntry.containsKey("coding") && ((Boolean) featureEntry.get("coding"))) {
                    cds = ((Boolean) featureEntry.get("coding"));
                }
            }
            //noinspection unchecked
            String matchKey = (String) featureEntry
                    .keySet()
                    .stream()
                    .filter(key -> key.toString().startsWith("match_"))
                    .findFirst()
                    .orElse(null);
            if (matchKey == null || featureEntry.get(matchKey) == null) {
                throw new MusialException(EXCEPTION_PREFIX + " Failed to find .gff attribute key/value pair to match feature " + featureKey);
            }
            Map<String, String> annotations;
            if (featureEntry.containsKey("annotations")) {
                //noinspection unchecked
                annotations = (Map<String, String>) featureEntry.get("annotations");
            } else {
                annotations = new HashMap<>();
            }
            addFeature(
                    (String) featureKey,
                    matchKey.replace("match_", ""),
                    (String) featureEntry.get(matchKey),
                    optionalPdbFile,
                    cds,
                    annotations
            );
        }
    }

    /**
     * Initializes a {@link Sample} object with the specified parameters and adds it to the samples list.
     *
     * @param name        {@link String}; The internal name to use for the sample.
     * @param vcfFile     {@link File} object pointing to a .vcf format file.
     * @param annotations {@link java.util.HashMap} of {@link String} key/pair values; sample meta information.
     * @throws MusialException If the specified .vcf file does not exist or can not be accessed.
     */
    private void addSample(String name, File vcfFile, Map<String, String> annotations) throws MusialException {
        if (vcfFile.isFile()) {
            Sample sample = new Sample(vcfFile, name);
            for (Map.Entry<String, String> annotation : annotations.entrySet()) {
                sample.addAnnotation(annotation.getKey(), annotation.getValue());
            }
            this.samples.put(name, sample);
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Failed to access .vcf file " + vcfFile);
        }
    }

    /**
     * Initializes a {@link Feature} object with the specified parameters and adds it to the features list.
     *
     * @param name        {@link String}; The internal name to use for the feature.
     * @param matchKey    {@link String}; The key of the attribute in the specified .gff format reference annotation to match the feature from.
     * @param matchValue  {@link String}; The value of the attribute in the specified .gff format reference annotation to match the feature from.
     * @param asCds       {@link Boolean}; Whether to consider feature as cds, independent of provided structure.
     * @param pdbFile     {@link File}; Optional object pointing to a .pdb format file yielding a protein structure derived for the (gene) feature.
     * @param annotations {@link java.util.HashMap} of {@link String} key/pair values; feature meta information.
     * @throws MusialException If the initialization of the {@link Feature} fails; If the specified .gff reference annotation or .pdb protein file can not be read; If the specified feature is not found or parsed multiple times from the reference annotation.
     */
    private void addFeature(String name, String matchKey, String matchValue, File pdbFile, boolean asCds, Map<String, String> annotations)
            throws MusialException {
        FeatureList matchedFeatures = this.referenceFeatures.selectByAttribute(matchKey, matchValue);
        if (matchedFeatures.size() == 0) {
            throw new MusialException(EXCEPTION_PREFIX + " Failed to match feature " + name + " with attribute key/value pair " +
                    matchKey + "=" + matchValue);
        } else if (matchedFeatures.size() > 1) {
            throw new MusialException(EXCEPTION_PREFIX + " Feature " + name + " was matched multiple times ");
        } else {
            FeatureI matchedFeature = matchedFeatures.get(0);
            Location featureCoordinates = matchedFeature.location();
            String featureChromosome = matchedFeature.seqname();
            /* FIXME: Starting positions are shifted by minus one, i.e. the returned values do not match with the ones of
                      the `gff` files. This is currently fixed in the Feature.java class.
            */
            Feature feature;
            if (asCds || pdbFile != null) {
                feature = new FeatureCoding(name, featureChromosome, featureCoordinates.getBegin(), featureCoordinates.getEnd(), "coding");
                if (pdbFile != null) {
                    try {
                        ((FeatureCoding) feature).setStructure(pdbFile);
                    } catch (IOException e) {
                        throw new MusialException(EXCEPTION_PREFIX + " Failed to parse structure from file " + pdbFile.getAbsolutePath() + " for feature " + feature.name);
                    }
                }
                ((FeatureCoding) feature).setCodingSequence(
                        Bio.translateNucSequence(
                                referenceSequence.getSubsequenceAt(feature.chromosome, feature.start, feature.end).getBaseString(),
                                true,
                                true,
                                feature.isSense
                        )
                );
            } else {
                feature = new Feature(name, featureChromosome, featureCoordinates.getBegin(), featureCoordinates.getEnd(), "non_coding");
            }
            for (Map.Entry<String, String> annotation : annotations.entrySet()) {
                feature.addAnnotation(annotation.getKey(), annotation.getValue());
            }
            this.features.put(name, feature);
        }
    }

}