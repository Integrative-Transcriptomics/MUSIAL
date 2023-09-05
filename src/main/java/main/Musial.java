package main;

import cli.BuildConfiguration;
import cli.CLIParser;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Sets;
import datastructure.*;
import exceptions.MusialException;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.util.TriConsumer;
import runnables.AlleleAnalyzer;
import runnables.SampleAnalyzer;
import runnables.SnpEffAnnotator;
import utility.Bio;
import utility.IO;
import utility.Logger;
import utility.MusialStorageFactory;

import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * Main class of MUSIAL (MUlti Sample varIant AnaLysis), a tool to calculate SNV, gene, and whole genome alignments,
 * together with other relevant statistics based on vcf files.
 *
 * @author Alexander Seitz
 * @author Simon Hackl
 * @version 2.2
 * @since 1.0
 */
@SuppressWarnings("DuplicatedCode")
public final class Musial {

    /**
     * {@link String} representing the name of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String NAME = "";
    /**
     * {@link String} representing the name of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String VERSION = "";
    /**
     * {@link String} representing author contact of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String CONTACT = "";
    /**
     * {@link String} representing license information of the software; parsed from `/src/main/resources/info.properties`.
     */
    public static String LICENSE = "";
    /**
     * Specifies the task to execute.
     */
    public static String TASK = null;
    /**
     * Original system output stream.
     */
    public static final PrintStream ORIGINAL_OUT_STREAM = System.out;
    /**
     * Original system error stream.
     */
    public static final PrintStream ORIGINAL_ERR_STREAM = System.err;
    /**
     * Alternative output stream to ignore logging.
     */
    public static final PrintStream EMPTY_STREAM = new PrintStream(new OutputStream() {
        public void write(int b) {
        }
    });
    /**
     * The number of threads to use.
     */
    public static int THREADS = 1;
    /**
     * Whether to compress output.
     */
    private static boolean COMPRESS = false;
    /**
     * Project wide formatter to convert decimal numbers to strings.
     */
    public static final DecimalFormat DECIMAL_FORMATTER = new DecimalFormat("#.####", DecimalFormatSymbols.getInstance(Locale.US));
    /**
     * Program wide accessor to command line interface parser.
     */
    @SuppressWarnings("FieldCanBeLocal")
    private static CLIParser CLI_PARSER;
    /**
     * Store start of execution time. Only for internal use.
     */
    @SuppressWarnings({"FieldCanBeLocal", "unused"})
    private static long START_TIME;
    /**
     * Store end of execution time. Only for internal use.
     */
    @SuppressWarnings({"FieldCanBeLocal", "unused"})
    private static long END_TIME;

    /**
     * Main method of MUSIAL; loads meta-information and invokes methods dependent on the command line arguments.
     *
     * @param args {@link String[]} comprising the passed command line arguments.
     */
    public static void main(String[] args) {
        try {
            START_TIME = System.currentTimeMillis();
            loadMetadata();
            if (args.length == 0) {
                System.out.println("No arguments were specified. Call `java -jar " + Musial.NAME + "-" + Musial.VERSION + ".jar [-h|--help]` for more information.");
                System.exit(0);
            } else {
                TASK = args[0];
                CLI_PARSER = new CLIParser(args);
                if (Musial.TASK.equalsIgnoreCase(MusialTasks.BUILD.name())) {
                    // Logger logger = LoggerFactory.getLogger(GFF3Reader.class);
                    Logger.logStatus("Execute task `build`");
                    if (CLI_PARSER.arguments.hasOption("t")) THREADS = Integer.parseInt(CLI_PARSER.arguments.getOptionValue("t"));
                    COMPRESS = CLI_PARSER.arguments.hasOption("k");
                    BuildConfiguration buildConfiguration = new BuildConfiguration(CLI_PARSER.buildConfiguration);
                    build(buildConfiguration);
                } else if (Musial.TASK.equalsIgnoreCase(MusialTasks.VIEW_FEATURES.name())) {
                    Logger.logStatus("Execute task `view_features`");
                    viewFeatures(CLI_PARSER.arguments);
                } else if (Musial.TASK.equalsIgnoreCase(MusialTasks.VIEW_SAMPLES.name())) {
                    Logger.logStatus("Execute task `view_samples`");
                    viewSamples(CLI_PARSER.arguments);
                } else if (Musial.TASK.equalsIgnoreCase(MusialTasks.VIEW_VARIANTS.name())) {
                    Logger.logStatus("Execute task `view_variants`");
                    viewVariants(CLI_PARSER.arguments);
                } else if (Musial.TASK.equalsIgnoreCase(MusialTasks.EXPORT_TABLE.name())) {
                    Logger.logStatus("Execute task `export_table`");
                    exportTable(CLI_PARSER.arguments);
                } else if (Musial.TASK.equalsIgnoreCase(MusialTasks.EXPORT_SEQUENCE.name())) {
                    Logger.logStatus("Execute task `export_sequence`");
                    exportSequence(CLI_PARSER.arguments);
                } else {
                    System.exit(0);
                }
            }
        } catch (Exception e) {
            Logger.logError(e.getMessage());
            e.printStackTrace();
            System.exit(-1);
        } finally {
            END_TIME = System.currentTimeMillis();
            // Logger.logStatus("Total execution time: " + ((END_TIME - START_TIME) / 1000));
        }
    }

    /**
     * Load and store metadata, such as the software title and version from `/src/main/resources/info.properties` and prints the information to stdout.
     *
     * @throws IOException If any metadata can not be loaded.
     */
    private static void loadMetadata() throws IOException {
        Properties properties = new Properties();
        InputStream in = Musial.class.getResourceAsStream("/info.properties");
        properties.load(in);
        Musial.NAME = properties.getProperty("name");
        Musial.VERSION = properties.getProperty("version");
        Musial.CONTACT = properties.getProperty("contact");
        Musial.LICENSE = properties.getProperty("license");
        // Print information to stdout.
        Logger.printSoftwareInfo();
    }

    /**
     * Generates a new or updates an existing MUSIAL storage JSON file based on the specifications parsed from a build configuration JSON file.
     *
     * @param parameters {@link BuildConfiguration} instance yielding parameter specification for the MUSIAL build task.
     * @throws InterruptedException Thrown if any {@link Runnable} implementing class instance is interrupted by the user.
     * @throws IOException          Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialException      If any method fails wrt. biological context, i.e. parsing of unknown symbols; If any method fails wrt. internal logic, i.e. assignment of proteins to genomes; If any input or output file is missing or unable to being generated.
     */
    private static void build(BuildConfiguration parameters) throws InterruptedException, MusialException, IOException {
        // Build new feature variants object.
        MusialStorage musialStorage = MusialStorageFactory.build(parameters);

        // Add feature information.
        for (Feature feature : parameters.features.values()) {
            musialStorage.addFeature(feature);
        }

        // Collect information about sample variants.
        Logger.logStatus("Process variant calls");
        ExecutorService executor = Executors.newFixedThreadPool(Musial.THREADS);
        for (Sample sample : parameters.samples.values()) {
            sample.imputeVcfFileReader();
            musialStorage.addSample(sample);
            executor.execute(new SampleAnalyzer(sample, musialStorage.getFeatureNameIterator(), musialStorage));
        }
        executor.shutdown();
        //noinspection ResultOfMethodCallIgnored
        executor.awaitTermination(30, TimeUnit.MINUTES);
        System.gc();

        // Run SnpEff annotation of variants.
        Logger.logStatus("Annotate variant calls with SnpEff");
        SnpEffAnnotator snpEffAnnotator = new SnpEffAnnotator(
                new File("./temp/"),
                new File("./temp/variants.vcf"),
                parameters.referenceSequenceFile,
                parameters.referenceFeaturesFile,
                musialStorage
        );
        snpEffAnnotator.run();
        System.gc();


        // Infer coding feature variant/proteoform information; For all coding features.
        Logger.logStatus("Infer proteoforms");
        Iterator<String> featureNameIterator = musialStorage.getFeatureNameIterator();
        String featureName;
        Feature feature;
        FeatureCoding featureCoding;
        while (featureNameIterator.hasNext()) {
            featureName = featureNameIterator.next();
            feature = musialStorage.getFeature(featureName);
            if (feature instanceof FeatureCoding) {
                featureCoding = (FeatureCoding) feature;
                executor = Executors.newFixedThreadPool(Musial.THREADS);
                Iterator<String> alleleNameIterator = featureCoding.getAlleleNameIterator();
                while (alleleNameIterator.hasNext()) {
                    executor.execute(new AlleleAnalyzer(featureCoding, featureCoding.getAllele(alleleNameIterator.next()), musialStorage));
                }
                executor.shutdown();
                //noinspection ResultOfMethodCallIgnored
                executor.awaitTermination(30, TimeUnit.MINUTES);
            }
        }
        System.gc();

        // Compute statistics for annotation.
        Logger.logStatus("Annotate statistics");
        float noSamples = musialStorage.getNumberOfSamples();
        float featureSequenceLength;
        long noOccurrences;
        featureNameIterator = musialStorage.getFeatureNameIterator();
        int noSubstitutions;
        int noInsertions;
        int noDeletions;
        while (featureNameIterator.hasNext()) {
            // Feature level statistics.
            featureName = featureNameIterator.next();
            feature = musialStorage.getFeature(featureName);
            featureSequenceLength = Math.abs((feature.end - feature.start) + 1);
            noSubstitutions = 0;
            noInsertions = 0;
            noDeletions = 0;
            int variantLength;
            for (Integer variantPosition : feature.getNucleotideVariantPositions()) {
                for (Map.Entry<String, VariantAnnotation> variantEntry : feature.getNucleotideVariants(variantPosition).entrySet()) {
                    noOccurrences = variantEntry.getValue().getPropertyKeys().stream().filter(s -> s.contains(MusialConstants.VARIANT_OCCURRENCE_SAMPLE_PREFIX)).count();
                    variantLength = variantEntry.getKey().length();
                    variantEntry.getValue().addProperty(
                            MusialConstants.FREQUENCY,
                            DECIMAL_FORMATTER.format((noOccurrences / noSamples) * 100)
                    );
                    if (variantEntry.getKey().contains("-")) {
                        noDeletions += variantLength;
                    } else if (variantEntry.getValue().getProperty(MusialConstants.REFERENCE_CONTENT).contains("-")) {
                        noInsertions += variantLength;
                    } else {
                        noSubstitutions += variantLength;
                    }
                }
            }
            feature.addAnnotation(MusialConstants.NUMBER_SUBSTITUTIONS, String.valueOf(noSubstitutions));
            feature.addAnnotation(MusialConstants.NUMBER_INSERTIONS, String.valueOf(noInsertions));
            feature.addAnnotation(MusialConstants.NUMBER_DELETIONS, String.valueOf(noDeletions));
            feature.addAnnotation(MusialConstants.VARIABLE_POSITIONS, DECIMAL_FORMATTER.format(
                    ((noSubstitutions + noInsertions + noDeletions) / featureSequenceLength) * 100
            ));
            Iterator<String> formNameIterator = feature.getAlleleNameIterator();
            String formName;
            Form form;
            while (formNameIterator.hasNext()) {
                // Allele level statistics.
                formName = formNameIterator.next();
                form = feature.getAllele(formName);
                noSubstitutions = 0;
                noInsertions = 0;
                noDeletions = 0;
                form.addAnnotation(
                        MusialConstants.FREQUENCY,
                        DECIMAL_FORMATTER.format((form.getOccurrence().size() / noSamples) * 100)
                );
                try {
                    for (Map.Entry<Integer, String> variantEntry : feature.getNucleotideVariants(formName).entrySet()) {
                        variantLength = variantEntry.getValue().length();
                        if (variantEntry.getValue().contains("-")) {
                            noDeletions += variantLength;
                        } else if (feature.getNucleotideVariantAnnotation(variantEntry.getKey(), variantEntry.getValue()).getProperty(MusialConstants.REFERENCE_CONTENT).contains("-")) {
                            noInsertions += variantLength;
                        } else {
                            noSubstitutions += variantLength;
                        }
                    }
                } catch (NullPointerException e) {
                    e.printStackTrace();
                    System.exit(1);
                }
                form.addAnnotation(MusialConstants.NUMBER_SUBSTITUTIONS, String.valueOf(noSubstitutions));
                form.addAnnotation(MusialConstants.NUMBER_INSERTIONS, String.valueOf(noInsertions));
                form.addAnnotation(MusialConstants.NUMBER_DELETIONS, String.valueOf(noDeletions));
                form.addAnnotation(MusialConstants.VARIABLE_POSITIONS, DECIMAL_FORMATTER.format(
                        ((noSubstitutions + noInsertions + noDeletions) / featureSequenceLength) * 100
                ));
            }
            if (feature.isCoding()) {
                featureCoding = (FeatureCoding) feature;
                formNameIterator = featureCoding.getProteoformNameIterator();
                while (formNameIterator.hasNext()) {
                    // Proteoform level statistics.
                    formName = formNameIterator.next();
                    form = featureCoding.getProteoform(formName);
                    noOccurrences = 0;
                    noSubstitutions = 0;
                    noInsertions = 0;
                    noDeletions = 0;
                    for (String alleleName : form.getOccurrence()) {
                        noOccurrences += featureCoding.getAllele(alleleName).getOccurrence().size();
                    }
                    form.addAnnotation(
                            MusialConstants.FREQUENCY,
                            DECIMAL_FORMATTER.format((noOccurrences / noSamples) * 100)
                    );
                    if (!form.hasAnnotation(MusialConstants.PROTEOFORM_DIFFERENTIAL_SEQUENCE)) {
                        for (Map.Entry<Integer, String> variantEntry : featureCoding.getAminoacidVariants(formName).entrySet()) {
                            variantLength = variantEntry.getValue().length();
                            if (variantEntry.getValue().contains("-")) {
                                noDeletions += variantLength;
                            } else if (featureCoding.getAminoacidVariantAnnotation(variantEntry.getKey(), variantEntry.getValue()).getProperty(MusialConstants.REFERENCE_CONTENT).contains("-")) {
                                noInsertions += variantLength;
                            } else {
                                noSubstitutions += variantLength;
                            }
                        }
                        form.addAnnotation(MusialConstants.NUMBER_SUBSTITUTIONS, String.valueOf(noSubstitutions));
                        form.addAnnotation(MusialConstants.NUMBER_INSERTIONS, String.valueOf(noInsertions));
                        form.addAnnotation(MusialConstants.NUMBER_DELETIONS, String.valueOf(noDeletions));
                        form.addAnnotation(MusialConstants.VARIABLE_POSITIONS, DECIMAL_FORMATTER.format(
                                ((noSubstitutions + noInsertions + noDeletions) / featureCoding.getCodingSequence().length()) * 100L
                        ));
                    }
                }
                // Variant level statistics.
                for (Integer variantPosition : featureCoding.getAminoacidVariantPositions()) {
                    for (Map.Entry<String, VariantAnnotation> variantEntry : featureCoding.getAminoacidVariants(variantPosition).entrySet()) {
                        noOccurrences = variantEntry.getValue().getPropertyKeys().stream().filter(s -> s.contains(MusialConstants.VARIANT_OCCURRENCE_SAMPLE_PREFIX)).count();
                        variantEntry.getValue().addProperty(
                                MusialConstants.FREQUENCY,
                                DECIMAL_FORMATTER.format((noOccurrences / noSamples) * 100)
                        );
                    }
                }
            }
        }

        // Sample level statistics.
        Iterator<String> sampleNameIterator = musialStorage.getSampleNameIterator();
        Sample sample;
        while (sampleNameIterator.hasNext()) {
            sample = musialStorage.getSample(sampleNameIterator.next());
            noSubstitutions = 0;
            noInsertions = 0;
            noDeletions = 0;
            featureNameIterator = musialStorage.getFeatureNameIterator();
            String alleleName;
            while (featureNameIterator.hasNext()) {
                feature = musialStorage.getFeature(featureNameIterator.next());
                alleleName = sample.getAnnotation(MusialConstants.SAMPLE_ANNOTATION_ALLELE_PREFIX + feature.name);
                noSubstitutions += Integer.parseInt(feature.getAllele(alleleName).getAnnotation(MusialConstants.NUMBER_SUBSTITUTIONS));
                noInsertions += Integer.parseInt(feature.getAllele(alleleName).getAnnotation(MusialConstants.NUMBER_INSERTIONS));
                noDeletions += Integer.parseInt(feature.getAllele(alleleName).getAnnotation(MusialConstants.NUMBER_DELETIONS));
            }
            sample.addAnnotation(MusialConstants.NUMBER_SUBSTITUTIONS, String.valueOf(noSubstitutions));
            sample.addAnnotation(MusialConstants.NUMBER_INSERTIONS, String.valueOf(noInsertions));
            sample.addAnnotation(MusialConstants.NUMBER_DELETIONS, String.valueOf(noDeletions));
        }
        System.gc();

        // Write updated storage to file.
        Logger.logStatus("Dump to file `" + parameters.output.getAbsolutePath() + "`");
        String outputPath = parameters.output.getAbsolutePath();
        File outputFile = new File(outputPath);
        if (!outputFile.exists()) {
            IO.generateFile(outputFile);
        }
        musialStorage.dump(outputFile, COMPRESS);
    }

    /**
     * Views all features of an existing MUSIAL storage as tsv formatted table.
     *
     * @param parameters {@link CommandLine} instance yielding parameter specification for the MUSIAL view_features task.
     * @throws IOException     If reading of the input MUSIAL storage fails.
     * @throws MusialException If printing the output fails.
     */
    private static void viewFeatures(CommandLine parameters) throws IOException, MusialException {
        // Parse feature variants object.
        MusialStorage musialStorage = MusialStorageFactory.load(new File(parameters.getOptionValue("I")));

        // Remove all features that should not be viewed.
        HashSet<String> featuresToView;
        if (parameters.getOptionValues("f") == null) {
            featuresToView = new HashSet<>(Sets.newHashSet(musialStorage.getFeatureNameIterator()));
        } else {
            featuresToView = (HashSet<String>) Arrays.stream(parameters.getOptionValues("f")).collect(Collectors.toSet());
        }

        // Collect feature information.
        Iterator<String> featureNameIterator = musialStorage.getFeatureNameIterator();
        LinkedHashMap<String, ArrayList<String>> outputMap = new LinkedHashMap<>() {{
            put("name", new ArrayList<>(featuresToView.size()));
            put("chromosome", new ArrayList<>(featuresToView.size()));
            put("start", new ArrayList<>(featuresToView.size()));
            put("end", new ArrayList<>(featuresToView.size()));
            put("strand", new ArrayList<>(featuresToView.size()));
            put("number_of_alleles", new ArrayList<>(featuresToView.size()));
            put("number_of_proteoforms", new ArrayList<>(featuresToView.size()));
        }};
        Feature feature;
        int featureIndex = 0;
        while (featureNameIterator.hasNext()) {
            feature = musialStorage.getFeature(featureNameIterator.next());
            if (!featuresToView.contains(feature.name))
                continue;
            outputMap.get("name").add(feature.name);
            outputMap.get("chromosome").add(feature.chromosome);
            outputMap.get("start").add(String.valueOf(feature.start));
            outputMap.get("end").add(String.valueOf(feature.end));
            outputMap.get("strand").add(feature.isSense ? "+" : "-");
            outputMap.get("number_of_alleles").add(String.valueOf(feature.getAlleleCount()));
            outputMap.get("number_of_proteoforms").add(feature.isCoding() ? String.valueOf(((FeatureCoding) feature).getProteoformCount()) : "null");
            HashSet<String> featureAnnotationKeys = new HashSet<>() {{
                add("name");
                add("chromosome");
                add("start");
                add("end");
                add("strand");
                add("number_of_alleles");
                add("number_of_proteoforms");
            }};
            for (Map.Entry<String, String> annotation : feature.getAnnotations()) {
                String annotationKey = annotation.getKey();
                transferAnnotation(outputMap, featureIndex, featureAnnotationKeys, annotationKey, annotation.getValue());
            }
            for (String globalAnnotationKey : outputMap.keySet()) {
                if (!featureAnnotationKeys.contains(globalAnnotationKey)) {
                    outputMap.get(globalAnnotationKey).add("null");
                }
            }
            featureIndex += 1;
        }

        // Construct and output string content from collected information
        outputView(outputMap, featureIndex, parameters.hasOption("o") ? parameters.getOptionValue("o") : null);
    }

    /**
     * Views all samples of an existing MUSIAL storage as tsv formatted table.
     *
     * @param parameters {@link CommandLine} instance yielding parameter specification for the MUSIAL view_samples task.
     * @throws IOException     If reading of the input MUSIAL storage fails.
     * @throws MusialException If printing the output fails.
     */
    private static void viewSamples(CommandLine parameters) throws IOException, MusialException {
        // Parse feature variants object.
        MusialStorage musialStorage = MusialStorageFactory.load(new File(parameters.getOptionValue("I")));

        // Remove all samples that should not be viewed.
        HashSet<String> samplesToView;
        if (parameters.getOptionValues("s") == null) {
            samplesToView = new HashSet<>(Sets.newHashSet(musialStorage.getSampleNameIterator()));
        } else {
            samplesToView = (HashSet<String>) Arrays.stream(parameters.getOptionValues("s")).collect(Collectors.toSet());
        }

        // Collect sample information.
        Iterator<String> sampleNameIterator = musialStorage.getSampleNameIterator();
        LinkedHashMap<String, ArrayList<String>> outputMap = new LinkedHashMap<>() {{
            put("name", new ArrayList<>(samplesToView.size()));
            put(MusialConstants.NUMBER_SUBSTITUTIONS, new ArrayList<>(samplesToView.size()));
            put(MusialConstants.NUMBER_INSERTIONS, new ArrayList<>(samplesToView.size()));
            put(MusialConstants.NUMBER_DELETIONS, new ArrayList<>(samplesToView.size()));
        }};
        Sample sample;
        int sampleIndex = 0;
        while (sampleNameIterator.hasNext()) {
            sample = musialStorage.getSample(sampleNameIterator.next());
            if (!(samplesToView.contains(sample.name)))
                continue;
            outputMap.get("name").add(sample.name);
            // outputMap.get("count.substitutions").add(String.valueOf(sample.getAnnotation(MusialConstants.NUMBER_SUBSTITUTIONS)));
            // outputMap.get("count.insertions").add(String.valueOf(sample.getAnnotation(MusialConstants.NUMBER_INSERTIONS)));
            // outputMap.get("count.deletions").add(String.valueOf(sample.getAnnotation(MusialConstants.NUMBER_DELETIONS)));
            HashSet<String> sampleAnnotationKeys = new HashSet<>() {{
                add("name");
            }};
            for (Map.Entry<String, String> annotation : sample.getAnnotations()) {
                String annotationKey = annotation.getKey();
                transferAnnotation(outputMap, sampleIndex, sampleAnnotationKeys, annotationKey, annotation.getValue());
            }
            for (String globalAnnotationKey : outputMap.keySet()) {
                if (!sampleAnnotationKeys.contains(globalAnnotationKey)) {
                    outputMap.get(globalAnnotationKey).add("null");
                }
            }
            sampleIndex += 1;
        }

        // Construct and output string content from collected information
        outputView(outputMap, sampleIndex, parameters.hasOption("o") ? parameters.getOptionValue("o") : null);
    }

    /**
     * Views all variants of an existing MUSIAL storage as tsv formatted table.
     *
     * @param parameters {@link CommandLine} instance yielding parameter specification for the MUSIAL view_variants task.
     * @throws IOException     If reading of the input MUSIAL storage fails.
     * @throws MusialException If printing the output fails.
     */
    private static void viewVariants(CommandLine parameters) throws IOException, MusialException {
        // Parse feature variants object.
        MusialStorage musialStorage = MusialStorageFactory.load(new File(parameters.getOptionValue("I")));

        // Remove all features that should not be viewed.
        HashSet<String> featuresToView;
        if (parameters.getOptionValues("f") == null) {
            featuresToView = new HashSet<>(Sets.newHashSet(musialStorage.getFeatureNameIterator()));
        } else {
            featuresToView = (HashSet<String>) Arrays.stream(parameters.getOptionValues("f")).collect(Collectors.toSet());
        }

        // Collect all samples that should be viewed.
        HashSet<String> samplesToView;
        if (parameters.getOptionValues("s") == null) {
            samplesToView = new HashSet<>(Sets.newHashSet(musialStorage.getSampleNameIterator()));
        } else {
            samplesToView = (HashSet<String>) Arrays.stream(parameters.getOptionValues("s")).collect(Collectors.toSet());
        }

        // Iterate over features to collect variants.
        Iterator<String> featureNameIterator = musialStorage.getFeatureNameIterator();
        Feature feature;
        VariantAnnotation variantAnnotation;
        String contentMode = parameters.hasOption("c") ? parameters.getOptionValue("c") : "nucleotide";
        String variantPropertyValue;
        int variantsCount = featuresToView
                .stream()
                .map(featureName -> {
                    if (contentMode.equals("nucleotide")) {
                        return musialStorage.getFeature(featureName).getNucleotideVariantPositions().size();
                    } else {
                        return ((FeatureCoding) musialStorage.getFeature(featureName)).getAminoacidVariantPositions().size();
                    }
                })
                .reduce(0, Integer::sum);
        LinkedHashMap<String, ArrayList<String>> outputMap = new LinkedHashMap<>() {{
            put("position", new ArrayList<>(variantsCount));
            put("type", new ArrayList<>(variantsCount));
            put("reference_content", new ArrayList<>(variantsCount));
            put("alternate_content", new ArrayList<>(variantsCount));
            put("frequency", new ArrayList<>(variantsCount));
            put("feature", new ArrayList<>(variantsCount));
            put("occurrence", new ArrayList<>(variantsCount));
        }};
        int variantIndex = 0;
        while (featureNameIterator.hasNext()) {
            feature = musialStorage.getFeature(featureNameIterator.next());
            if (!(featuresToView.contains(feature.name)))
                continue;
            NavigableSet<Integer> variantPositions = new TreeSet<>();
            ConcurrentSkipListMap<String, VariantAnnotation> variants = new ConcurrentSkipListMap<>();
            if (contentMode.equals("aminoacid") && !(feature instanceof FeatureCoding))
                continue;
            switch (contentMode) {
                case "nucleotide" -> variantPositions = feature.getNucleotideVariantPositions();
                case "aminoacid" -> variantPositions = ((FeatureCoding) feature).getAminoacidVariantPositions();
            }
            for (Integer variantPosition : variantPositions) {
                switch (contentMode) {
                    case "nucleotide" -> variants = feature.getNucleotideVariants(variantPosition);
                    case "aminoacid" -> variants = ((FeatureCoding) feature).getAminoacidVariants(variantPosition);
                }
                for (Map.Entry<String, VariantAnnotation> variant : variants.entrySet()) {
                    // Check if variant has to be skipped wrt. samples.
                    variantAnnotation = variant.getValue();
                    boolean skip = false;
                    for (String sampleToView : samplesToView) {
                        if (variantAnnotation.hasProperty(MusialConstants.VARIANT_OCCURRENCE_SAMPLE_PREFIX + sampleToView)) {
                            skip = true;
                            break;
                        }
                    }
                    if (!skip)
                        continue;
                    ArrayList<String> occurrence = new ArrayList<>();
                    HashSet<String> variantAnnotationKeys = new HashSet<>() {{
                        add("position");
                        add("type");
                        add("reference_content");
                        add("alternate_content");
                        add("frequency");
                        add("feature");
                        add("occurrence");
                    }};
                    transferAnnotation(outputMap, variantIndex, variantAnnotationKeys, "position", String.valueOf(variantPosition));
                    String type = "substitution";
                    if (variant.getKey().contains(String.valueOf(Bio.DELETION_AA1))) {
                        type = "deletion";
                    } else if (variantAnnotation.getProperty(MusialConstants.REFERENCE_CONTENT).contains(String.valueOf(Bio.DELETION_AA1))) {
                        type = "insertion";
                    }
                    transferAnnotation(outputMap, variantIndex, variantAnnotationKeys, "type", type);
                    transferAnnotation(outputMap, variantIndex, variantAnnotationKeys, "alternate_content", variant.getKey());
                    transferAnnotation(outputMap, variantIndex, variantAnnotationKeys, "feature", feature.name);
                    for (String variantPropertyKey : variantAnnotation.getPropertyKeys()) {
                        variantPropertyValue = variantAnnotation.getProperty(variantPropertyKey);
                        if (variantPropertyKey.equals(MusialConstants.REFERENCE_CONTENT)) {
                            transferAnnotation(outputMap, variantIndex, variantAnnotationKeys, "reference_content", variantPropertyValue);
                        } else if (variantPropertyKey.equals(MusialConstants.FREQUENCY)) {
                            transferAnnotation(outputMap, variantIndex, variantAnnotationKeys, "frequency", variantPropertyValue);
                        } else if (variantPropertyKey.startsWith(MusialConstants.VARIANT_OCCURRENCE_SAMPLE_PREFIX)) {
                            switch (contentMode) {
                                case "nucleotide" ->
                                        occurrence.add(variantPropertyKey.replace(MusialConstants.VARIANT_OCCURRENCE_SAMPLE_PREFIX, "") + ":" + variantPropertyValue);
                                case "aminoacid" ->
                                        occurrence.add(variantPropertyKey.replace(MusialConstants.VARIANT_OCCURRENCE_SAMPLE_PREFIX, ""));
                            }
                        } else {
                            if (variantPropertyKey.startsWith(MusialConstants.SNP_EFF_PROPERTY_PREFIX))
                                variantPropertyKey = variantPropertyKey.replace(MusialConstants.SNP_EFF_PROPERTY_PREFIX, "SnpEff.");
                            transferAnnotation(outputMap, variantIndex, variantAnnotationKeys, variantPropertyKey, variantPropertyValue);
                        }
                    }
                    transferAnnotation(outputMap, variantIndex, variantAnnotationKeys, "occurrence", String.join(",", occurrence));
                    for (String globalAnnotationKey : outputMap.keySet()) {
                        if (!variantAnnotationKeys.contains(globalAnnotationKey)) {
                            outputMap.get(globalAnnotationKey).add("null");
                        }
                    }
                    variantIndex += 1;
                }
            }
        }
        // Construct and output string content from collected information
        outputView(outputMap, outputMap.get("position").size(), parameters.hasOption("o") ? parameters.getOptionValue("o") : null);
    }

    /**
     * Exports a variant table in tsv format wrt. a single feature.
     *
     * @param parameters {@link CommandLine} instance yielding parameter specification for the MUSIAL export_table task.
     * @throws IOException     If reading of the input MUSIAL storage fails.
     * @throws MusialException If any computation fails.
     */
    private static void exportTable(CommandLine parameters) throws IOException, MusialException {
        // Parse feature variants object.
        MusialStorage musialStorage = MusialStorageFactory.load(new File(parameters.getOptionValue("I")));

        // Store content mode, grouping and rejected output options.
        String contentMode = parameters.hasOption("c") ? parameters.getOptionValue("c") : "nucleotide";
        boolean group = parameters.hasOption("g");
        boolean rejectedAsReference = parameters.hasOption("x");
        boolean conservedSites = parameters.hasOption("k");

        // Check if specified feature exists and supports the chosen content mode.
        Feature feature = musialStorage.getFeature(parameters.getOptionValue("F"));
        if (feature == null) {
            throw new MusialException("The specified feature " + parameters.getOptionValue("F") + " could not be found in the provided MUSIAL storage file; Available features are " + String.join(", ", ImmutableList.copyOf(musialStorage.getFeatureNameIterator())));
        }
        if (contentMode.equals("aminoacid") && !(feature.isCoding())) {
            throw new MusialException("The specified feature " + parameters.getOptionValue("F") + " is no coding feature, but `aminoacid` content mode was specified.");
        }

        // Definition of helper functions.
        Function<String, String> sampleNameToFormName = s -> {
            if (contentMode.equals("nucleotide")) {
                return musialStorage.getSample(s).getAnnotation(MusialConstants.SAMPLE_ANNOTATION_ALLELE_PREFIX + feature.name);
            } else {
                return musialStorage.getSample(s).getAnnotation(MusialConstants.SAMPLE_ANNOTATION_PROTEOFORM_PREFIX + feature.name);
            }
        };
        Function<String, String> sampleNameToRepresentativeName = s -> {
            if (contentMode.equals("nucleotide")) {
                return feature.getAllele(sampleNameToFormName.apply(s)).getOccurrence().first();
            } else {
                return feature.getAllele(((FeatureCoding) feature).getProteoform(sampleNameToFormName.apply(s)).getOccurrence().first()).getOccurrence().first();
            }
        };
        Function<String, LinkedHashSet<String>> getFormOccurrence = s -> {
            LinkedHashSet<String> occurrence = new LinkedHashSet<>();
            if (contentMode.equals("nucleotide")) {
                occurrence.addAll(feature.getAllele(s).getOccurrence());
            } else {
                for (String alleleName : ((FeatureCoding) feature).getProteoform(s).getOccurrence()) {
                    occurrence.addAll(feature.getAllele(alleleName).getOccurrence());
                }
            }
            return occurrence;
        };

        // Parse reference sequence information, if provided.
        char[] referenceSequenceChars = null;
        if (conservedSites && contentMode.equals("nucleotide") && parameters.hasOption("r")) {
            File referenceSequenceFile = new File(parameters.getOptionValue("r"));
            IndexedFastaSequenceFile referenceSequence;
            try {
                FastaSequenceIndexCreator.create(referenceSequenceFile.toPath(), false);
            } catch (SAMException e) {
                // Exception is raised if file already exists. This can be ignored.
            }
            referenceSequence = new IndexedFastaSequenceFile(referenceSequenceFile.toPath());
            musialStorage.setReference(referenceSequence);
            referenceSequenceChars = musialStorage.getReferenceSequenceOfFeature(feature.name).toCharArray();
        } else if (conservedSites && contentMode.equals("aminoacid")) {
            referenceSequenceChars = ((FeatureCoding) feature).getCodingSequence().toCharArray();
        }

        // Collect all samples that should be exported.
        HashSet<String> sampleNames;
        if (parameters.getOptionValues("s") == null) {
            sampleNames = new HashSet<>(Sets.newHashSet(musialStorage.getSampleNameIterator()));
        } else {
            sampleNames = (HashSet<String>) Arrays.stream(parameters.getOptionValues("s")).collect(Collectors.toSet());
        }
        HashSet<String> representativeSampleNames = new HashSet<>();
        for (String sampleName : sampleNames) {
            representativeSampleNames.add(sampleNameToRepresentativeName.apply(sampleName));
        }

        // Generate variants table for feature.
        TreeMap<String, LinkedHashMap<String, String>> variantsTable = constructVariantsTable(
                musialStorage,
                feature.name,
                representativeSampleNames,
                contentMode,
                rejectedAsReference,
                conservedSites
        );
        System.gc();
        // Build output.
        try (BufferedWriter outputWriter = new BufferedWriter(new FileWriter(parameters.getOptionValue("O")))) {
            // Generate header.
            if (group) {
                outputWriter.write("Position\tReference");
                for (String representativeName : representativeSampleNames) {
                    outputWriter.write("\t" + sampleNameToFormName.apply(representativeName));
                }
                outputWriter.newLine();
            } else {
                outputWriter.write("\nPosition\tReference");
                for (String representativeName : representativeSampleNames) {
                    outputWriter.write(
                            "\t" + String.join(
                                    "\t",
                                    getFormOccurrence.apply(sampleNameToFormName.apply(representativeName)).stream().filter(sampleNames::contains).collect(Collectors.toSet())
                            )
                    );
                }
                outputWriter.newLine();
            }

            // Generate variants rows.
            Iterator<String> positionsIterator = variantsTable.navigableKeySet().iterator();
            Iterator<LinkedHashMap<String, String>> contentIterator = variantsTable.values().iterator();
            LinkedHashMap<String, String> variantContents;
            String variantContent;
            String position;
            do {
                position = positionsIterator.next();
                // Position.
                outputWriter.write(position + "\t");
                // Reference content.
                if (conservedSites) {
                    //noinspection ConstantConditions
                    outputWriter
                            .write(referenceSequenceChars[Integer.parseInt(position.split("\\.")[0]) - (contentMode.equals("nucleotide") ? feature.start : 1)] + "\t");
                } else {
                    outputWriter
                            .write(variantsTable.get(position).get(MusialConstants.REFERENCE_ID) + "\t");
                }
                // Variant content.
                variantContents = contentIterator.next();
                if (group) {
                    for (String representativeName : representativeSampleNames) {
                        if (variantContents == null) {
                            variantContent = ".";
                        } else {
                            variantContent = variantContents.getOrDefault(representativeName, ".");
                        }
                        outputWriter.write(variantContent + "\t");
                    }
                } else {
                    for (String representativeName : representativeSampleNames) {
                        if (variantContents == null) {
                            variantContent = ".";
                        } else {
                            variantContent = variantContents.getOrDefault(representativeName, ".");
                        }
                        for (String ignored : getFormOccurrence.apply(sampleNameToFormName.apply(representativeName)).stream().filter(sampleNames::contains).collect(Collectors.toSet())) {
                            outputWriter.write(variantContent + "\t");
                        }
                    }
                }
                outputWriter.newLine();
                outputWriter.flush();
            } while (positionsIterator.hasNext());
        }
        Logger.logStatus("Done writing table to file `" + new File(parameters.getOptionValue("O")).getAbsolutePath() + "`");
    }

    /**
     * Exports sequences in fasta format wrt. a single feature.
     *
     * @param parameters {@link CommandLine} instance yielding parameter specification for the MUSIAL export_sequence task.
     * @throws IOException     If reading of the input MUSIAL storage fails.
     * @throws MusialException If any computation fails.
     */
    private static void exportSequence(CommandLine parameters) throws IOException, MusialException {
        // Parse feature variants object.
        MusialStorage musialStorage = MusialStorageFactory.load(new File(parameters.getOptionValue("I")));

        // Store content mode, grouping and rejected output options.
        String contentMode = parameters.hasOption("c") ? parameters.getOptionValue("c") : "nucleotide";
        boolean group = parameters.hasOption("g");
        boolean rejectedAsReference = parameters.hasOption("x");
        boolean conservedSites = parameters.hasOption("k");
        boolean aligned = parameters.hasOption("a");

        // Check if specified feature exists and supports the chosen content mode.
        Feature feature = musialStorage.getFeature(parameters.getOptionValue("F"));
        if (feature == null) {
            throw new MusialException("The specified feature " + parameters.getOptionValue("F") + " could not be found in the provided MUSIAL storage file; Available features are " + String.join(", ", ImmutableList.copyOf(musialStorage.getFeatureNameIterator())));
        }
        if (contentMode.equals("aminoacid") && !(feature.isCoding())) {
            throw new MusialException("The specified feature " + parameters.getOptionValue("F") + " is no coding feature, but `aminoacid` content mode was specified.");
        }

        // Definition of helper functions.
        Function<String, String> sampleNameToFormName = s -> {
            if (s.equals(MusialConstants.REFERENCE_ID))
                return MusialConstants.REFERENCE_ID;
            if (contentMode.equals("nucleotide")) {
                return musialStorage.getSample(s).getAnnotation(MusialConstants.SAMPLE_ANNOTATION_ALLELE_PREFIX + feature.name);
            } else {
                return musialStorage.getSample(s).getAnnotation(MusialConstants.SAMPLE_ANNOTATION_PROTEOFORM_PREFIX + feature.name);
            }
        };
        Function<String, String> sampleNameToRepresentativeName = s -> {
            if (contentMode.equals("nucleotide")) {
                return feature.getAllele(sampleNameToFormName.apply(s)).getOccurrence().first();
            } else {
                return feature.getAllele(((FeatureCoding) feature).getProteoform(sampleNameToFormName.apply(s)).getOccurrence().first()).getOccurrence().first();
            }
        };
        Function<String, LinkedHashSet<String>> getFormOccurrence = s -> {
            LinkedHashSet<String> occurrence = new LinkedHashSet<>();
            if (contentMode.equals("nucleotide")) {
                if ( !feature.hasAllele(s) )
                    return occurrence;
                occurrence.addAll(feature.getAllele(s).getOccurrence());
            } else {
                if ( !((FeatureCoding) feature).hasProteoform(s) )
                    return occurrence;
                for (String alleleName : ((FeatureCoding) feature).getProteoform(s).getOccurrence()) {
                    if ( !feature.hasAllele(alleleName) )
                        continue;
                    occurrence.addAll(feature.getAllele(alleleName).getOccurrence());
                }
            }
            return occurrence;
        };
        TriConsumer<BufferedWriter, String, String> writeSequence = (w, s, n) -> {
            String sequence;
            if (aligned)
                sequence = String.join("\n", Splitter.fixedLength(80).split(s));
            else
                sequence = String.join("\n", Splitter.fixedLength(80).split(s.replace("-", "")));
            try {
                w.write(">" + n);
                w.newLine();
                w.write(sequence);
                w.newLine();
            } catch (IOException e) {
                Logger.logWarning("Failed to write sequence of entry " + n + ";" + e.getMessage());
            }
        };

        // Parse reference sequence information, if provided.
        char[] referenceSequenceChars = null;
        if (conservedSites && contentMode.equals("nucleotide") && parameters.hasOption("r")) {
            File referenceSequenceFile = new File(parameters.getOptionValue("r"));
            IndexedFastaSequenceFile referenceSequence;
            try {
                FastaSequenceIndexCreator.create(referenceSequenceFile.toPath(), false);
            } catch (SAMException e) {
                // Exception is raised if file already exists. This can be ignored.
            }
            referenceSequence = new IndexedFastaSequenceFile(referenceSequenceFile.toPath());
            musialStorage.setReference(referenceSequence);
            referenceSequenceChars = musialStorage.getReferenceSequenceOfFeature(feature.name).toCharArray();
        } else if (conservedSites && contentMode.equals("aminoacid")) {
            referenceSequenceChars = ((FeatureCoding) feature).getCodingSequence().toCharArray();
        }

        // Collect all samples that should be exported.
        HashSet<String> sampleNames;
        if (parameters.getOptionValues("s") == null) {
            sampleNames = new HashSet<>(Sets.newHashSet(musialStorage.getSampleNameIterator()));
        } else {
            sampleNames = (HashSet<String>) Arrays.stream(parameters.getOptionValues("s")).collect(Collectors.toSet());
        }
        HashSet<String> representativeSampleNames = new HashSet<>();
        representativeSampleNames.add(MusialConstants.REFERENCE_ID);
        for (String sampleName : sampleNames) {
            if (sampleNameToFormName.apply(sampleName).equals(MusialConstants.REFERENCE_ID))
                continue;
            representativeSampleNames.add(sampleNameToRepresentativeName.apply(sampleName));
        }

        // Generate variants table for feature.
        TreeMap<String, LinkedHashMap<String, String>> variantsTable = constructVariantsTable(
                musialStorage,
                feature.name,
                representativeSampleNames,
                contentMode,
                rejectedAsReference,
                conservedSites
        );
        System.gc();

        // Build output.
        Iterator<String> representativeSampleIterator = representativeSampleNames.iterator();
        try (BufferedWriter outputWriter = new BufferedWriter(new FileWriter(parameters.getOptionValue("O")))) {
            do {
                Iterator<String> positionsIterator = variantsTable.navigableKeySet().iterator();
                Iterator<LinkedHashMap<String, String>> contentIterator = variantsTable.values().iterator();
                LinkedHashMap<String, String> variantContents;
                StringBuilder sequenceBuilder = new StringBuilder((feature.end - feature.start) + 1);
                String representativeName = representativeSampleIterator.next();
                String variantContent;
                String position;
                do {
                    position = positionsIterator.next();
                    variantContents = contentIterator.next();
                    if (variantContents == null) {
                        //noinspection ConstantConditions
                        variantContent = String.valueOf(referenceSequenceChars[Integer.parseInt(position.split("\\.")[0]) - (contentMode.equals("nucleotide") ? feature.start : 1)]);
                    } else {
                        variantContent = variantContents.getOrDefault(representativeName, variantContents.get(MusialConstants.REFERENCE_ID));
                    }
                    sequenceBuilder.append(variantContent);
                } while (positionsIterator.hasNext());
                if (group) {
                    writeSequence.accept(outputWriter, sequenceBuilder.toString(), sampleNameToFormName.apply(representativeName));
                } else {
                    if (representativeName.equals(MusialConstants.REFERENCE_ID)) {
                        writeSequence.accept(outputWriter, sequenceBuilder.toString(), representativeName);
                        for (String sampleName : getFormOccurrence.apply(MusialConstants.REFERENCE_ID).stream().filter(sampleNames::contains).collect(Collectors.toSet())) {
                            writeSequence.accept(outputWriter, sequenceBuilder.toString(), sampleName);
                        }
                    } else {
                        for (String sampleName : getFormOccurrence.apply(sampleNameToFormName.apply(representativeName)).stream().filter(sampleNames::contains).collect(Collectors.toSet())) {
                            writeSequence.accept(outputWriter, sequenceBuilder.toString(), sampleName);
                        }
                    }
                }
                outputWriter.flush();
            } while (representativeSampleIterator.hasNext());
        }
        Logger.logStatus("Done writing sequences to file `" + new File(parameters.getOptionValue("O")).getAbsolutePath() + "`");
    }

    /**
     * Internal method to output view_x tasks of MUSIAL.
     *
     * @param content    The content for which a tsv format string is to be built.
     * @param entryIndex The number of entries to iterate over.
     * @param output     {@link String} representing a file path to write the output to; If set to null, output will be printed to console.
     * @throws MusialException If I/O operation fails in case of a provided output file.
     */
    private static void outputView(LinkedHashMap<String, ArrayList<String>> content, int entryIndex, String output) throws MusialException {
        Logger.logStatus("Write output to " + (output == null ? "stdout" : output));
        // Build string content from collected information.
        StringBuilder outputBuilder = new StringBuilder(entryIndex);
        ArrayList<String> perEntryAnnotation;
        outputBuilder.append(String.join("\t", content.keySet())).append("\n");
        for (int i = 0; i < entryIndex; i++) {
            perEntryAnnotation = new ArrayList<>();
            for (ArrayList<String> annotationValues : content.values()) {
                perEntryAnnotation.add(annotationValues.get(i));
            }
            outputBuilder.append(String.join("\t", perEntryAnnotation)).append("\n");
        }

        // Write output to file or console.
        if (output != null) {
            IO.writeFile(new File(output), outputBuilder.toString());
            Logger.logStatus("Done writing output to file `" + new File(output).getAbsolutePath() + "`");
        } else {
            System.out.println(outputBuilder);
        }
    }

    /**
     * Internal method to transfers an annotation value to a container map. This is used to supplement missing annotation values in any object to consider with null values.
     *
     * @param annotations         The container to add annotations to.
     * @param entryIndex          The number of existing entries in the container.
     * @param entryAnnotationKeys Annotation keys of the entry to consider.
     * @param annotationKey       The key of the annotation.
     * @param annotationValue     The value of the annotation.
     */
    private static void transferAnnotation(LinkedHashMap<String, ArrayList<String>> annotations, int entryIndex, HashSet<String> entryAnnotationKeys, String annotationKey, String annotationValue) {
        entryAnnotationKeys.add(annotationKey);
        if (annotations.containsKey(annotationKey)) {
            annotations.get(annotationKey).add(annotationValue);
        } else {
            annotations.put(annotationKey, new ArrayList<>(entryIndex));
            for (int i = 0; i < entryIndex; i++) {
                annotations.get(annotationKey).add("null");
            }
            annotations.get(annotationKey).add(annotationValue);
        }
    }

    /**
     * Internal method to construct a two-layer map structure used for export of tables and sequences.
     * <br>
     * The first layer stores positions in the format x+y, where x represent true genomic positions and y represent insertion positions.
     * <br>
     * The second layer stores sample names mapping to variant contents.
     *
     * @param musialStorage       {@link MusialStorage} to build from.
     * @param featureName         Name of the feature to build the table of; Should be accepted for {@link MusialStorage#getFeature(String)}.
     * @param sampleNames         Set of sample names to restrict the construction to; Should be accepted for {@link MusialStorage#getSample(String)}.
     * @param contentMode         Either `aminoacid` or `nucleotide`; Specifies the variants to query.
     * @param rejectedAsReference If rejected variants should be added as reference or ambiguous base symbol.
     * @param conservedSites      If non-variant sites should be considered.
     * @return Two-layer map structure.
     */
    private static TreeMap<String, LinkedHashMap<String, String>> constructVariantsTable(MusialStorage musialStorage, String featureName, HashSet<String> sampleNames, String contentMode, boolean rejectedAsReference, boolean conservedSites) {
        TreeMap<String, LinkedHashMap<String, String>> variantsTable = new TreeMap<>((p1, p2) -> {
            String[] f1 = p1.split("\\.");
            String[] f2 = p2.split("\\.");
            int x1 = Integer.parseInt(f1[0]);
            int y1 = Integer.parseInt(f1[1]);
            int x2 = Integer.parseInt(f2[0]);
            int y2 = Integer.parseInt(f2[1]);
            if (x1 == x2) {
                return Integer.compare(y1, y2);
            } else {
                return Integer.compare(x1, x2);
            }
        });
        Supplier<String> getAmbiguousSymbol = ( ) -> switch (contentMode) {
            case "nucleotide" -> Bio.ANY_NUC;
            case "aminoacid" -> String.valueOf(Bio.ANY_AA1);
            default -> ".";
        };
        Feature feature = musialStorage.getFeature(featureName);
        NavigableSet<Integer> variantPositions = new TreeSet<>();
        ConcurrentSkipListMap<String, VariantAnnotation> variants = new ConcurrentSkipListMap<>();
        VariantAnnotation variantAnnotation;
        char[] variantChars;
        char[] referenceChars;
        String[] perSampleVariantProperties;
        boolean rejected;
        TriConsumer<String, String, String> insertToVariantsTable = (position, sampleName, content) -> {
            if (!variantsTable.containsKey(position)) {
                variantsTable.put(position, new LinkedHashMap<>(sampleNames.size()));
            }
            variantsTable.get(position).put(sampleName, content);
        };
        switch (contentMode) {
            case "nucleotide" -> variantPositions.addAll(feature.getNucleotideVariantPositions());
            case "aminoacid" -> variantPositions.addAll(((FeatureCoding) feature).getAminoacidVariantPositions());
        }
        if (conservedSites) {
            int positionStart = contentMode.equals("nucleotide") ? feature.start : 1;
            int positionEnd = contentMode.equals("nucleotide") ? feature.end : ((FeatureCoding) feature).getCodingSequence().length();
            for (int position = positionStart; position <= positionEnd; position++) {
                variantPositions.add(position);
            }
        }
        for (int variantPosition : variantPositions) {
            switch (contentMode) {
                case "nucleotide" -> variants = feature.getNucleotideVariants(variantPosition);
                case "aminoacid" -> variants = ((FeatureCoding) feature).getAminoacidVariants(variantPosition);
            }
            if (variants == null) {
                variantsTable.put(variantPosition + ".0", null);
                continue;
            }
            boolean anyVariantOccurrence = false; // Flag, if any variant at this position occurs in the sample set.
            boolean variantOccurrence; // Flag, if specific variant at this position occurs in the sample set.
            for (Map.Entry<String, VariantAnnotation> variant : variants.entrySet()) {
                variantOccurrence = false;
                variantChars = variant.getKey().toCharArray();
                variantAnnotation = variant.getValue();
                referenceChars = variantAnnotation.getProperty(MusialConstants.REFERENCE_CONTENT).toCharArray();
                String type = "substitution";
                if (variant.getKey().contains(String.valueOf(Bio.DELETION_AA1))) {
                    type = "deletion";
                } else if (variantAnnotation.getProperty(MusialConstants.REFERENCE_CONTENT).contains(String.valueOf(Bio.DELETION_AA1))) {
                    type = "insertion";
                }
                for (String sampleName : sampleNames) {
                    if (variantAnnotation.hasProperty(MusialConstants.VARIANT_OCCURRENCE_SAMPLE_PREFIX + sampleName)) {
                        if (!anyVariantOccurrence)
                            anyVariantOccurrence = true;
                        if (!variantOccurrence)
                            variantOccurrence = true;
                        perSampleVariantProperties = variantAnnotation.getProperty(MusialConstants.VARIANT_OCCURRENCE_SAMPLE_PREFIX + sampleName).split(":");
                        rejected = contentMode.equals("nucleotide") && Objects.equals(perSampleVariantProperties[1], "true");
                        switch (type) {
                            case "substitution":
                                if (rejected) {
                                    insertToVariantsTable.accept(variantPosition + ".0", sampleName, rejectedAsReference ? String.valueOf(referenceChars[0]) : getAmbiguousSymbol.get());
                                } else {
                                    insertToVariantsTable.accept(variantPosition + ".0", sampleName, variant.getKey());
                                }
                                break;
                            case "deletion":
                                for (int i = 0; i < variantChars.length; i++) {
                                    if (rejected) {
                                        insertToVariantsTable.accept((variantPosition + i + 1) + ".0", sampleName, rejectedAsReference ? String.valueOf(referenceChars[i]) : getAmbiguousSymbol.get());
                                    } else {
                                        insertToVariantsTable.accept((variantPosition + i + 1) + ".0", sampleName, String.valueOf(variantChars[i]));
                                    }
                                }
                                break;
                            case "insertion":
                                for (int i = 0; i < variantChars.length; i++) {
                                    if (rejected) {
                                        insertToVariantsTable.accept((variantPosition) + "." + (i + 1), sampleName, rejectedAsReference ? String.valueOf(referenceChars[i]) : getAmbiguousSymbol.get());
                                    } else {
                                        insertToVariantsTable.accept((variantPosition) + "." + (i + 1), sampleName, String.valueOf(variantChars[i]));
                                    }
                                }
                                break;
                        }
                    }
                }
                if (variantOccurrence) { // If this specific variant occurs in the sample set, add reference content to table.
                    switch (type) {
                        case "substitution":
                            insertToVariantsTable.accept(variantPosition + ".0", MusialConstants.REFERENCE_ID, String.valueOf(referenceChars[0]));
                            break;
                        case "deletion":
                            for (int i = 0; i < referenceChars.length; i++) {
                                insertToVariantsTable.accept((variantPosition + i + 1) + ".0", MusialConstants.REFERENCE_ID, String.valueOf(referenceChars[i]));
                            }
                            break;
                        case "insertion":
                            for (int i = 0; i < referenceChars.length; i++) {
                                insertToVariantsTable.accept((variantPosition) + "." + (i + 1), MusialConstants.REFERENCE_ID, String.valueOf(referenceChars[i]));
                            }
                            break;
                    }
                }
            }

            if (!anyVariantOccurrence && conservedSites) // If no variant is contained in the sample set, add null content.
                variantsTable.put(variantPosition + ".0", null);
        }
        return variantsTable;
    }
}