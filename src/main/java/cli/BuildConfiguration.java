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
import utility.IO;
import utility.Logger;
import utility.SequenceOperations;
import utility.Validation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Parse parameters for `build` task.
 *
 * @author Simon Hackl
 */
public final class BuildConfiguration {
    /**
     * {@link HashMap} of {@link String}/{@link Sample} pairs.
     */
    public final HashMap<String, Sample> samples = new HashMap<>();
    /**
     * {@link HashMap} of {@link String}/{@link Feature} pairs.
     */
    public final HashMap<String, Feature> features = new HashMap<>();
    /**
     * Optional positions to exclude from the analysis.
     */
    public final HashMap<String, TreeSet<Integer>> excludedPositions = new HashMap<>();
    /**
     * Optional explicit variants to exclude from the analysis.
     */
    public final HashMap<String, HashMap<Integer, HashSet<String>>> excludedVariants = new HashMap<>();
    /**
     * Prefix to use for exception messages.
     */
    @SuppressWarnings("FieldCanBeLocal")
    private final String EXCEPTION_PREFIX = "(Task `build` Argument Validation)";
    /**
     * Minimum (read) coverage in order to accept an allele.
     */
    public Double minimalCoverage;
    /**
     * Minimum allele frequency (with respect to reads) in order to accept an allele.
     */
    public Double minimalFrequency;
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
     * {@link File} object specifying the output directory.
     */
    public File output;

    /**
     * Constructor of the {@link BuildConfiguration} class.
     * <p>
     * Used to parse and validate command line interface input. Parsed arguments are stored - and accessible by other
     * components - as class properties.
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
        // Parse minimalFrequency
        if (parameters.containsKey("minimalFrequency")
                && Validation.isPercentage(String.valueOf(parameters.get("minimalFrequency")))) {
            this.minimalFrequency = (Double) parameters.get("minimalFrequency");
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid or missing parameter `minimalFrequency`; expected float between 0.0 and 1.0.");
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

        // Parse reference features (.gff3); Important: .gff does not work!
        if (!parameters.containsKey("referenceFeaturesFile")) {
            throw new MusialException(EXCEPTION_PREFIX + " Missing parameter `referenceFeaturesFile`; expected path to .gff3 file.");
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
            for (String contig : excludedPositions.keySet()) {
                this.excludedPositions.put(contig, new TreeSet<>());
                for (Double excludedPosition : excludedPositions.get(contig)) {
                    this.excludedPositions.get(contig).add(excludedPosition.intValue());
                }
            }
        }

        // Parse excluded variants
        if (parameters.containsKey("excludedVariants")) {
            //noinspection rawtypes,unchecked
            LinkedTreeMap<String, ArrayList<String>> excludedVariants = (LinkedTreeMap) parameters.get("excludedVariants");
            for (String contig : excludedVariants.keySet()) {
                this.excludedVariants.put(contig, new HashMap<>());
                for (String excludedVariantEntry : excludedVariants.get(contig)) {
                    try {
                        String[] excludedVariantFields = excludedVariantEntry.split(":");
                        Integer position = Integer.parseInt(excludedVariantFields[0]);
                        String excludedVariant = excludedVariantFields[1];
                        this.excludedVariants.get(contig).putIfAbsent(
                                position, new HashSet<>()
                        );
                        this.excludedVariants.get(contig).get(position).add(excludedVariant);
                    } catch (Exception e) {
                        Logger.logWarning("Failed to parse variant to ignore on feature " + contig + ": " + excludedVariantEntry + ".");
                    }
                }
            }
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
            boolean cds = (featureEntry.containsKey("coding") && ((Boolean) featureEntry.get("coding")));
            //noinspection unchecked
            String matchKey = (String) featureEntry
                    .keySet()
                    .stream()
                    .filter(key -> key.toString().startsWith("match_"))
                    .findFirst()
                    .orElse(null);
            if (matchKey == null || featureEntry.get(matchKey) == null) {
                throw new MusialException(EXCEPTION_PREFIX + " Failed to find GFF attribute key/value pair to match feature " + featureKey);
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
                    cds,
                    annotations
            );
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
            Sample sample = new Sample(vcfFile, name, this.features.size());
            for (Map.Entry<String, String> annotation : annotations.entrySet()) {
                sample.addInfo(annotation.getKey(), annotation.getValue());
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
     * @param matchKey    {@link String}; The key of the attribute in the specified GFF format reference annotation to match the feature from.
     * @param matchValue  {@link String}; The value of the attribute in the specified GFF format reference annotation to match the feature from.
     * @param coding      {@link Boolean}; Whether to consider feature as cds, independent of provided structure.
     * @param annotations {@link java.util.HashMap} of {@link String} key/pair values; feature meta information.
     * @throws MusialException If the initialization of the {@link Feature} fails; If the specified GFF reference annotation or .pdb protein file can not be read; If the specified feature is not found or parsed multiple times from the reference annotation.
     */
    private void addFeature(String name, String matchKey, String matchValue, boolean coding, Map<String, String> annotations)
            throws MusialException {
        FeatureList matchedFeatures = this.referenceFeatures.selectByAttribute(matchKey, matchValue);
        if (matchedFeatures.size() == 0) {
            throw new MusialException(EXCEPTION_PREFIX + " Failed to match feature " + name + " with attribute key/value pair " +
                    matchKey + "=" + matchValue);
        } else {
            if (matchedFeatures.size() > 1) {
                boolean isInvalid = false;
                int locationStart = matchedFeatures.get(0).location().bioStart();
                int locationEnd = matchedFeatures.get(0).location().bioEnd();
                for (int i = 1; i < matchedFeatures.size(); i++) {
                    if ((matchedFeatures.get(i).location().bioStart() != locationStart) || (matchedFeatures.get(i).location().bioEnd() != locationEnd)) {
                        isInvalid = true;
                    }
                }
                if (isInvalid) {
                    throw new MusialException(EXCEPTION_PREFIX + " Feature " + name + " was matched multiple times with different locations.");
                }
            }
            FeatureI matchedFeature = matchedFeatures.get(0);
            Location featureLocation = matchedFeature.location();
            String featureChromosome = matchedFeature.seqname();
            Feature feature;
            if (coding) {
                feature = new FeatureCoding(name, featureChromosome, featureLocation.bioStart(), featureLocation.bioEnd(), featureLocation.bioStrand(), "coding");
                ((FeatureCoding) feature).setCodingSequence(
                        SequenceOperations.translateNucSequence(
                                referenceSequence.getSubsequenceAt(feature.contig, feature.start, feature.end).getBaseString(),
                                true,
                                true,
                                feature.isSense
                        )
                );
            } else {
                feature = new Feature(name, featureChromosome, featureLocation.getBegin(), featureLocation.getEnd(), featureLocation.bioStrand(), "non_coding");
            }
            for (Map.Entry<String, String> annotation : annotations.entrySet()) {
                feature.addInfo(annotation.getKey(), annotation.getValue());
            }
            this.features.put(name, feature);
        }
    }

}