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
import org.apache.commons.compress.utils.Lists;
import org.apache.commons.lang3.SerializationUtils;
import org.apache.logging.log4j.util.TriConsumer;
import utility.*;

import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Main class of MUSIAL (MUlti Sample varIant AnaLysis), a tool to calculate SNV, gene, and whole genome alignments,
 * together with other relevant statistics based on vcf files.
 *
 * @author Alexander Seitz
 * @author Simon Hackl
 * @version 2.3
 * @since 1.0
 */
@SuppressWarnings("DuplicatedCode")
public final class Musial {

    /**
     * Original system output stream.
     */
    public static final PrintStream ORIGINAL_OUT_STREAM = System.out;
    /**
     * Alternative output stream to ignore logging.
     */
    public static final PrintStream EMPTY_STREAM = new PrintStream(new OutputStream() {
        public void write(int b) {
        }
    });
    /**
     * Project wide formatter to convert decimal numbers to strings.
     */
    public static final DecimalFormat DECIMAL_FORMATTER = new DecimalFormat("#.####", DecimalFormatSymbols.getInstance(Locale.US));
    /**
     * Whether to compress output. This can not be set by the user and is for debugging purposes only.
     */
    private static final boolean COMPRESS = true;
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
                if (Musial.TASK.equalsIgnoreCase(Tasks.BUILD.name())) {
                    // Logger logger = LoggerFactory.getLogger(GFF3Reader.class);
                    Logger.logStatus("Execute task `build`");
                    BuildConfiguration buildConfiguration = new BuildConfiguration(CLI_PARSER.buildConfiguration);
                    build(buildConfiguration);
                } else if (Musial.TASK.equalsIgnoreCase(Tasks.VIEW_FEATURES.name())) {
                    Logger.logStatus("Execute task `view_features`");
                    viewFeatures(CLI_PARSER.arguments);
                } else if (Musial.TASK.equalsIgnoreCase(Tasks.VIEW_SAMPLES.name())) {
                    Logger.logStatus("Execute task `view_samples`");
                    viewSamples(CLI_PARSER.arguments);
                } else if (Musial.TASK.equalsIgnoreCase(Tasks.VIEW_VARIANTS.name())) {
                    Logger.logStatus("Execute task `view_variants`");
                    viewVariants(CLI_PARSER.arguments);
                } else if (Musial.TASK.equalsIgnoreCase(Tasks.EXPORT_TABLE.name())) {
                    Logger.logStatus("Execute task `export_table`");
                    exportTable(CLI_PARSER.arguments);
                } else if (Musial.TASK.equalsIgnoreCase(Tasks.EXPORT_SEQUENCE.name())) {
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
            Logger.logStatus("Execution time: " + ((END_TIME - START_TIME) / 1000));
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
     * @throws IOException     Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialException If any method fails wrt. biological context, i.e. parsing of unknown symbols; If any method fails wrt. internal logic, i.e. assignment of proteins to genomes; If any input or output file is missing or unable to being generated.
     */
    private static void build(BuildConfiguration parameters) throws MusialException, IOException {
        // Build new feature variants object.
        MusialStorage musialStorage = MusialStorageFactory.build(parameters);
        Iterator<String> iterator;

        // Collect information about sample variants.
        Logger.logStatus("Process variant calls");
        iterator = musialStorage.getSampleNameIterator();
        Sample sample;
        do {
            sample = musialStorage.getSample(iterator.next());
            sample.imputeVcfFileReader();
            SampleAnalyzer.run(sample, musialStorage.getFeatureNameIterator(), musialStorage);
        } while (iterator.hasNext());
        System.gc();

        // Run SnpEff annotation of variants.
        Logger.logStatus("Annotate variants with SnpEff");
        String temporaryWorkingDir = "./tmp/";
        if (CLI_PARSER.arguments.hasOption("w"))
            temporaryWorkingDir = CLI_PARSER.arguments.getOptionValue("w");
        if (!temporaryWorkingDir.endsWith("/"))
            temporaryWorkingDir += "/";
        SnpEffAnnotator snpEffAnnotator = new SnpEffAnnotator(
                new File(temporaryWorkingDir),
                new File(temporaryWorkingDir + "variants.vcf"),
                parameters.referenceSequenceFile,
                parameters.referenceFeaturesFile,
                musialStorage
        );
        snpEffAnnotator.run();
        System.gc();

        // Infer coding feature variant/proteoform information; For all coding features.
        Logger.logStatus("Infer proteoforms (if any)");
        iterator = musialStorage.getFeatureNameIterator();
        Feature feature;
        FeatureCoding featureCoding;
        do {
            feature = musialStorage.getFeature(iterator.next());
            if (feature instanceof FeatureCoding) {
                featureCoding = (FeatureCoding) feature;
                AlleleAnalyzer.run(featureCoding, musialStorage);
            }
        } while (iterator.hasNext());
        System.gc();

        // Update info statistics.
        updateStatistics(musialStorage);

        // Write updated storage to file.
        musialStorage.dump(parameters.output.getAbsolutePath(), COMPRESS);
    }

    /**
     * Views all variants of an existing MUSIAL storage as tsv formatted table.
     *
     * @param parameters {@link CommandLine} instance yielding parameter specification for the MUSIAL view_variants task.
     * @throws IOException     If reading of the input MUSIAL storage fails.
     * @throws MusialException If printing the output fails.
     */
    private static void viewVariants(CommandLine parameters) throws IOException, MusialException {
        // Load MUSIAL storage and parameters.
        MusialStorage musialStorage = MusialStorageFactory.load(new File(parameters.getOptionValue("I")));
        Feature feature;
        VariantInformation variantInformation;
        String value;
        String variantContent;
        String contentType = parameters.hasOption("c") ? parameters.getOptionValue("c") : Constants.CONTENT_MODE_NUCLEOTIDE;
        int variantsCount;
        int rowCount = 0;
        HashSet<String> featureNames;
        HashSet<String> sampleNames;
        HashSet<String> naNames;
        //noinspection MismatchedQueryAndUpdateOfCollection
        LinkedHashSet<String> defaultKeys;
        Iterator<String> featureNameIterator;
        NavigableSet<Integer> variantPositions;
        LinkedHashMap<String, VariantInformation> variants;
        LinkedHashMap<String, ArrayList<String>> outputContainer;

        // Parse features to view.
        if (parameters.getOptionValues("f") == null) {
            featureNames = new HashSet<>(Sets.newHashSet(musialStorage.getFeatureNameIterator()));
        } else {
            featureNames = (HashSet<String>) Arrays.stream(parameters.getOptionValues("f")).collect(Collectors.toSet());
        }
        naNames = new HashSet<>();
        for (String featureName : featureNames) {
            if (!musialStorage.hasFeature(featureName)) {
                naNames.add(featureName);
                Logger.logWarning("Feature with name " + featureName + " is not stored in specified input.");
            }
            if (contentType.equals(Constants.CONTENT_MODE_AMINOACID) && !musialStorage.getFeature(featureName).isCoding())
                naNames.add(featureName);
        }
        featureNames.removeAll(naNames);

        // Parse samples to view.
        if (parameters.getOptionValues("s") == null) {
            sampleNames = new HashSet<>(Sets.newHashSet(musialStorage.getSampleNameIterator()));
        } else {
            sampleNames = (HashSet<String>) Arrays.stream(parameters.getOptionValues("s")).collect(Collectors.toSet());
        }
        naNames = new HashSet<>();
        for (String sampleName : sampleNames) {
            if (!musialStorage.hasSample(sampleName)) {
                naNames.add(sampleName);
                Logger.logWarning("Sample with name " + sampleName + " is not stored in specified input.");
            }
        }
        sampleNames.removeAll(naNames);

        // Iterate over features to collect variants.
        featureNameIterator = featureNames.iterator();
        variantsCount = featureNames
                .stream()
                .map(featureName -> {
                    if (contentType.equals(Constants.CONTENT_MODE_NUCLEOTIDE))
                        return musialStorage.getFeature(featureName).getNucleotideVariantPositions().size();
                    else
                        return ((FeatureCoding) musialStorage.getFeature(featureName)).getAminoacidVariantPositions().size();
                })
                .reduce(0, Integer::sum);
        defaultKeys = new LinkedHashSet<>() {{
            add("position");
            add("reference_content");
            add("alternate_content");
            add("feature");
            add("occurrence");
        }};
        outputContainer = new LinkedHashMap<>(defaultKeys.size());
        for (String defaultKey : defaultKeys)
            outputContainer.put(defaultKey, new ArrayList<>(variantsCount));
        while (featureNameIterator.hasNext()) {
            feature = musialStorage.getFeature(featureNameIterator.next());
            variants = new LinkedHashMap<>(variantsCount);
            if (contentType.equals(Constants.CONTENT_MODE_AMINOACID) && !feature.isCoding()) {
                Logger.logWarning("Cannot view aminoacid variants for non coding feature " + feature.name + ".");
                continue;
            }
            variantPositions = Constants.CONTENT_MODE_AMINOACID.equals(contentType) ? ((FeatureCoding) feature).getAminoacidVariantPositions() : feature.getNucleotideVariantPositions();
            for (Integer variantPosition : variantPositions) {
                switch (contentType) {
                    case Constants.CONTENT_MODE_NUCLEOTIDE -> variants = feature.getNucleotideVariantsAt(variantPosition);
                    case Constants.CONTENT_MODE_AMINOACID -> variants = ((FeatureCoding) feature).getAminoacidVariantsAt(variantPosition);
                }
                for (Map.Entry<String, VariantInformation> variant : variants.entrySet()) {
                    // Check if variant has to be skipped wrt. samples.
                    variantContent = variant.getKey();
                    variantInformation = variant.getValue();
                    if (sampleNames.stream().noneMatch(variantInformation::hasOccurrence))
                        continue;
                    transferContent(outputContainer, "position", String.valueOf(variantPosition), rowCount);
                    transferContent(outputContainer, "reference_content", variantInformation.referenceContent, rowCount);
                    transferContent(outputContainer, "alternate_content", variantContent, rowCount);
                    transferContent(outputContainer, "feature", feature.name, rowCount);
                    transferContent(outputContainer, "occurrence", String.join(",", variantInformation.getOccurrenceSet()), rowCount);
                    for (String key : variantInformation.getInfoKeys()) {
                        value = variantInformation.getInfo(key);
                        transferContent(outputContainer, key, value, rowCount);
                    }
                    for (String key : outputContainer.keySet()) {
                        if (!defaultKeys.contains(key) && !variantInformation.hasInfoKey(key)) {
                            outputContainer.get(key).add("null");
                        }
                    }
                    rowCount += 1;
                }
            }
        }
        // Construct and output string content from collected information
        outputView(outputContainer, rowCount, parameters.hasOption("o") ? parameters.getOptionValue("o") : null);
    }

    /**
     * Views all samples of an existing MUSIAL storage as tsv formatted table.
     *
     * @param parameters {@link CommandLine} instance yielding parameter specification for the MUSIAL view_samples task.
     * @throws IOException     If reading of the input MUSIAL storage fails.
     * @throws MusialException If printing the output fails.
     */
    private static void viewSamples(CommandLine parameters) throws IOException, MusialException {
        MusialStorage musialStorage = MusialStorageFactory.load(new File(parameters.getOptionValue("I")));
        Sample sample;
        String featureName;
        String key;
        String value;
        int samplesCount;
        int rowCount = 0;
        HashSet<String> sampleNames;
        HashSet<String> naNames;
        //noinspection MismatchedQueryAndUpdateOfCollection
        LinkedHashSet<String> defaultKeys;
        Iterator<String> sampleNameIterator;
        Iterator<String> featureNameIterator;
        LinkedHashMap<String, ArrayList<String>> outputContainer;

        // Parse samples to view.
        if (parameters.getOptionValues("s") == null) {
            sampleNames = new HashSet<>(Sets.newHashSet(musialStorage.getSampleNameIterator()));
        } else {
            sampleNames = (HashSet<String>) Arrays.stream(parameters.getOptionValues("s")).collect(Collectors.toSet());
        }
        naNames = new HashSet<>();
        for (String sampleName : sampleNames) {
            if (!musialStorage.hasSample(sampleName)) {
                naNames.add(sampleName);
                Logger.logWarning("Sample with name " + sampleName + " is not stored in specified input.");
            }
        }
        sampleNames.removeAll(naNames);
        samplesCount = sampleNames.size();

        // Collect sample information.
        sampleNameIterator = sampleNames.iterator();

        defaultKeys = new LinkedHashSet<>() {{
            add("name");
            add(Constants.NUMBER_OF_SUBSTITUTIONS);
            add(Constants.NUMBER_OF_INSERTIONS);
            add(Constants.NUMBER_OF_DELETIONS);
            add(Constants.NUMBER_OF_AMBIGUOUS);
        }};
        outputContainer = new LinkedHashMap<>(defaultKeys.size());
        for (String defaultKey : defaultKeys)
            outputContainer.put(defaultKey, new ArrayList<>(samplesCount));
        while (sampleNameIterator.hasNext()) {
            sample = musialStorage.getSample(sampleNameIterator.next());
            outputContainer.get("name").add(sample.name);
            featureNameIterator = musialStorage.getFeatureNameIterator();
            do {
                featureName = featureNameIterator.next();
                transferContent(outputContainer, "allele_" + featureName, Objects.isNull(sample.getAllele(featureName)) ? "null" : sample.getAllele(featureName), rowCount);
                transferContent(outputContainer, "proteoform_" + featureName, Objects.isNull(sample.getProteoform(featureName)) ? "null" : sample.getProteoform(featureName), rowCount);
            } while (featureNameIterator.hasNext());
            for (Map.Entry<String, String> info : sample.getInfoSet()) {
                key = info.getKey();
                value = info.getValue();
                transferContent(outputContainer, key, value, rowCount);
            }
            for (String k : outputContainer.keySet()) {
                if (!defaultKeys.contains(k) && !sample.hasInfoKey(k) && !k.startsWith("allele_") && !k.startsWith("proteoform_")) {
                    outputContainer.get(k).add("null");
                }
            }
            rowCount += 1;
        }

        // Construct and output string content from collected information
        outputView(outputContainer, rowCount, parameters.hasOption("o") ? parameters.getOptionValue("o") : null);
    }

    /**
     * Views all features of an existing MUSIAL storage as tsv formatted table.
     *
     * @param parameters {@link CommandLine} instance yielding parameter specification for the MUSIAL view_features task.
     * @throws IOException     If reading of the input MUSIAL storage fails.
     * @throws MusialException If printing the output fails.
     */
    private static void viewFeatures(CommandLine parameters) throws IOException, MusialException {
        MusialStorage musialStorage = MusialStorageFactory.load(new File(parameters.getOptionValue("I")));
        Feature feature;
        String key;
        String value;
        int featuresCount;
        int rowCount = 0;
        HashSet<String> featureNames;
        HashSet<String> naNames;
        //noinspection MismatchedQueryAndUpdateOfCollection
        LinkedHashSet<String> defaultKeys;
        Iterator<String> featureNameIterator;
        LinkedHashMap<String, ArrayList<String>> outputContainer;

        // Parse features to view.
        if (parameters.getOptionValues("f") == null) {
            featureNames = new HashSet<>(Sets.newHashSet(musialStorage.getFeatureNameIterator()));
        } else {
            featureNames = (HashSet<String>) Arrays.stream(parameters.getOptionValues("f")).collect(Collectors.toSet());
        }
        naNames = new HashSet<>();
        for (String featureName : featureNames) {
            if (!musialStorage.hasFeature(featureName)) {
                naNames.add(featureName);
                Logger.logWarning("Feature with name " + featureName + " is not stored in specified input.");
            }
        }
        featureNames.removeAll(naNames);
        featuresCount = featureNames.size();

        // Collect feature information.
        featureNameIterator = featureNames.iterator();
        defaultKeys = new LinkedHashSet<>() {{
            add("name");
            add("location");
            add("start");
            add("end");
            add("strand");
            add("number_of_alleles");
            add("number_of_proteoforms");
        }};
        outputContainer = new LinkedHashMap<>(defaultKeys.size());
        for (String defaultKey : defaultKeys)
            outputContainer.put(defaultKey, new ArrayList<>(featuresCount));
        while (featureNameIterator.hasNext()) {
            feature = musialStorage.getFeature(featureNameIterator.next());
            outputContainer.get("name").add(feature.name);
            outputContainer.get("location").add(feature.contig);
            outputContainer.get("start").add(String.valueOf(feature.start));
            outputContainer.get("end").add(String.valueOf(feature.end));
            outputContainer.get("strand").add(feature.isSense ? "+" : "-");
            outputContainer.get("number_of_alleles").add(String.valueOf(feature.getAlleleCount()));
            outputContainer.get("number_of_proteoforms").add(feature.isCoding() ? String.valueOf(((FeatureCoding) feature).getProteoformCount()) : "null");
            for (Map.Entry<String, String> info : feature.getInfoSet()) {
                key = info.getKey();
                value = info.getValue();
                transferContent(outputContainer, key, value, rowCount);
            }
            for (String k : outputContainer.keySet()) {
                if (!defaultKeys.contains(k) && !feature.hasInfoKey(k)) {
                    outputContainer.get(k).add("null");
                }
            }
            rowCount += 1;
        }

        // Construct and output string content from collected information
        outputView(outputContainer, rowCount, parameters.hasOption("o") ? parameters.getOptionValue("o") : null);
    }

    /**
     * Exports a variant table in tsv format wrt. a single feature.
     *
     * @param parameters {@link CommandLine} instance yielding parameter specification for the MUSIAL export_table task.
     * @throws IOException     If reading of the input MUSIAL storage fails.
     * @throws MusialException If any computation fails.
     */
    private static void exportTable(CommandLine parameters) throws IOException, MusialException {
        // Load MUSIAL storage.
        MusialStorage musialStorage = MusialStorageFactory.load(new File(parameters.getOptionValue("I")));
        // Store content mode.
        String contentMode = parameters.hasOption("c") ? parameters.getOptionValue("c") : Constants.CONTENT_MODE_NUCLEOTIDE;
        // Check if specified feature exists and supports the chosen content mode.
        Feature feature = musialStorage.getFeature(parameters.getOptionValue("F"));
        if (feature == null) {
            throw new MusialException("The specified feature " + parameters.getOptionValue("F") + " is not stored in the provided MUSIAL storage file; Available features are " + String.join(", ", ImmutableList.copyOf(musialStorage.getFeatureNameIterator())));
        }
        if (contentMode.equals(Constants.CONTENT_MODE_AMINOACID) && !(feature.isCoding())) {
            throw new MusialException("The specified feature " + parameters.getOptionValue("F") + " is no coding feature, but `aminoacid` content was specified.");
        }
        // Parse samples to view.
        ArrayList<String> sampleNames;
        if (parameters.getOptionValues("s") == null) {
            sampleNames = Lists.newArrayList(musialStorage.getSampleNameIterator());
        } else {
            sampleNames = Lists.newArrayList(Arrays.stream(parameters.getOptionValues("s")).iterator());
        }
        HashSet<String> naNames = new HashSet<>();
        for (String sampleName : sampleNames) {
            if (!musialStorage.hasSample(sampleName)) {
                naNames.add(sampleName);
                Logger.logWarning("Sample with name " + sampleName + " is not stored in specified input.");
            }
        }
        sampleNames.removeAll(naNames);
        int noSamples = sampleNames.size();
        Map<String, Integer> sampleNameIndex = IntStream.range(0, noSamples).boxed().collect(Collectors.toMap(sampleNames::get, Function.identity()));
        // Build output.
        try (BufferedWriter outputWriter = new BufferedWriter(new FileWriter(parameters.getOptionValue("O")), 32768)) {
            StringBuilder line = new StringBuilder();
            StringBuilder sampleVariantContentBuilder = new StringBuilder();
            String[] sampleContents;
            char ref;
            Runnable writeLine = () -> {
                try {
                    outputWriter.write(line.toString());
                    outputWriter.newLine();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                line.setLength(0);
            }; // Internal helper function to dump current line content to file.
            // Write header.
            line.append("Position").append("\t").append("Reference").append("\t").append(String.join("\t", sampleNames));
            writeLine.run();
            NavigableSet<Integer> variantPositions = contentMode.equals(Constants.CONTENT_MODE_AMINOACID) ? ((FeatureCoding) feature).getAminoacidVariantPositions() : feature.getNucleotideVariantPositions();
            Set<Map.Entry<String, VariantInformation>> variants;
            for (Integer position : variantPositions) {
                variants = (contentMode.equals(Constants.CONTENT_MODE_AMINOACID) ? ((FeatureCoding) feature).getAminoacidVariantsAt(position) : feature.getNucleotideVariantsAt(position)).entrySet();
                sampleContents = new String[noSamples];
                //noinspection OptionalGetWithoutIsPresent
                ref = variants.stream().findFirst().get().getValue().referenceContent.charAt(0);
                line.append(position).append("\t").append(ref).append("\t");
                for (String sampleName : sampleNames) {
                    if (contentMode.equals(Constants.CONTENT_MODE_NUCLEOTIDE))
                        sampleContents[sampleNameIndex.get(sampleName)] = musialStorage.getSample(sampleName).getCallVariantsStringAt(feature.name, position);
                    else {
                        sampleVariantContentBuilder.setLength(0);
                        for (Map.Entry<String, VariantInformation> variant : variants) {
                            if (variant.getValue().hasOccurrence(sampleName)) {
                                if (variant.getValue().getInfo(Constants.TYPE).contains(Constants.TYPE_SUBSTITUTION)) {
                                    if (!sampleVariantContentBuilder.isEmpty())
                                        sampleVariantContentBuilder.replace(0, 1, variant.getKey());
                                    else
                                        sampleVariantContentBuilder.append(variant.getKey());
                                } else if (variant.getValue().getInfo(Constants.TYPE).contains(Constants.TYPE_INSERTION)) {
                                    if (!sampleVariantContentBuilder.isEmpty())
                                        sampleVariantContentBuilder.insert(1, variant.getKey().substring(1));
                                    else
                                        sampleVariantContentBuilder.append(variant.getKey());
                                } else if (variant.getValue().getInfo(Constants.TYPE).contains(Constants.TYPE_DELETION)) {
                                    if (!sampleVariantContentBuilder.isEmpty())
                                        sampleVariantContentBuilder.append(variant.getKey().substring(1));
                                    else
                                        sampleVariantContentBuilder.append(variant.getKey());
                                } else {
                                    Logger.logWarning("Unable to process variant " + position + " " + variant.getKey() + " of type " + variant.getValue().getInfo(Constants.TYPE) + ".");
                                }
                            }
                        }
                        if (sampleVariantContentBuilder.isEmpty())
                            sampleVariantContentBuilder.append(Constants.CALL_INFO_NO_VARIANT);
                        sampleContents[sampleNameIndex.get(sampleName)] = sampleVariantContentBuilder.toString();
                    }
                }
                line.append(String.join("\t", sampleContents));
                writeLine.run();
            }
            outputWriter.flush();
        }
        Logger.logStatus("Done writing output to file `" + new File(parameters.getOptionValue("O")).getAbsolutePath() + "`");
    }

    /**
     * Exports sequences in fasta format wrt. a single feature.
     *
     * @param parameters {@link CommandLine} instance yielding parameter specification for the MUSIAL export_sequence task.
     * @throws IOException     If reading of the input MUSIAL storage fails.
     * @throws MusialException If any computation fails.
     */
    private static void exportSequence(CommandLine parameters) throws IOException, MusialException {
        // Load MUSIAL storage.
        MusialStorage musialStorage = MusialStorageFactory.load(new File(parameters.getOptionValue("I")));
        HashSet<String> naNames;
        // Store parameters.
        String contentMode = parameters.hasOption("c") ? parameters.getOptionValue("c") : Constants.CONTENT_MODE_NUCLEOTIDE;
        boolean merge = parameters.hasOption("m");
        boolean conserved = parameters.hasOption("k");
        boolean aligned = parameters.hasOption("a");
        // Check if specified feature exists and supports the chosen content mode.
        Feature feature = musialStorage.getFeature(parameters.getOptionValue("F"));
        if (feature == null)
            throw new MusialException("The specified feature " + parameters.getOptionValue("F") + " is not stored in the provided MUSIAL storage file; Available features are " + String.join(", ", ImmutableList.copyOf(musialStorage.getFeatureNameIterator())));
        if (contentMode.equals(Constants.CONTENT_MODE_AMINOACID) && !(feature.isCoding()))
            throw new MusialException("The specified feature " + parameters.getOptionValue("F") + " is no coding feature, but `aminoacid` content mode was specified.");
        // Parse samples to view.
        ArrayList<String> sampleNames;
        HashSet<String> sampleNamesSet;
        if (parameters.getOptionValues("s") == null)
            sampleNames = Lists.newArrayList(musialStorage.getSampleNameIterator());
        else
            sampleNames = Lists.newArrayList(Arrays.stream(parameters.getOptionValues("s")).iterator());
        naNames = new HashSet<>();
        for (String sampleName : sampleNames) {
            if (!musialStorage.hasSample(sampleName)) {
                naNames.add(sampleName);
                Logger.logWarning("Sample with name " + sampleName + " is not stored in specified input.");
            }
        }
        sampleNames.removeAll(naNames);
        sampleNamesSet = new HashSet<>(sampleNames);
        Function<VariantInformation, Boolean> variantOccurs = (variantInformation -> variantInformation.getOccurrenceSet().stream().anyMatch(sampleNamesSet::contains));

        // Construct reference table.
        LinkedHashMap<Integer, String> referenceTable;
        if (conserved) {
            String[] perPositionReferenceContent;
            NavigableSet<Integer> variantPositions;
            Function<Integer, Collection<VariantInformation>> getVariantInformation;
            if (contentMode.equals(Constants.CONTENT_MODE_NUCLEOTIDE) && parameters.hasOption("r")) {
                File referenceSequenceFile = new File(parameters.getOptionValue("r"));
                IndexedFastaSequenceFile referenceSequence;
                try {
                    FastaSequenceIndexCreator.create(referenceSequenceFile.toPath(), false);
                } catch (SAMException e) {
                    // Exception is raised if file already exists. This can be ignored.
                }
                referenceSequence = new IndexedFastaSequenceFile(referenceSequenceFile.toPath());
                musialStorage.setReference(referenceSequence);
                perPositionReferenceContent = musialStorage.getReferenceSequenceOfFeature(feature.name).split("");
                referenceTable = new LinkedHashMap<>(perPositionReferenceContent.length, 10);
                variantPositions = feature.getNucleotideVariantPositions();
                getVariantInformation = (position) -> feature.getNucleotideVariantsAt(position).values();
            } else if (contentMode.equals(Constants.CONTENT_MODE_AMINOACID)) {
                FeatureCoding featureCoding = (FeatureCoding) feature;
                perPositionReferenceContent = featureCoding.getCodingSequence().split("");
                referenceTable = new LinkedHashMap<>(perPositionReferenceContent.length, 10);
                variantPositions = featureCoding.getAminoacidVariantPositions();
                getVariantInformation = (position) -> featureCoding.getAminoacidVariantsAt(position).values();
            } else {
                throw new MusialException("Failed to export sequences with parameters --content=" + contentMode + " and --reference=" + parameters.getOptionValue("r") + ".");
            }
            for (int i = 0; i < perPositionReferenceContent.length; i++)
                if (contentMode.equals(Constants.CONTENT_MODE_NUCLEOTIDE))
                    referenceTable.put(feature.start + i, perPositionReferenceContent[i]);
                else
                    referenceTable.put(1 + i, perPositionReferenceContent[i]);
            for (Integer position : variantPositions)
                for (VariantInformation variantInformation : getVariantInformation.apply(position))
                    // Replaces reference content with maximal insertion length wrt. samples to use.
                    if (variantInformation.getInfo(Constants.TYPE).contains(Constants.TYPE_INSERTION) && variantOccurs.apply(variantInformation))
                        if (variantInformation.referenceContent.length() > referenceTable.getOrDefault(position, "").length())
                            referenceTable.put(position, variantInformation.referenceContent);
        } else {
            NavigableSet<Integer> variantPositions;
            Function<Integer, Collection<Map.Entry<String, VariantInformation>>> getVariantEntries;
            String alt;
            VariantInformation variantInformation;
            char[] chars;
            if (contentMode.equals(Constants.CONTENT_MODE_NUCLEOTIDE)) {
                variantPositions = feature.getNucleotideVariantPositions();
                referenceTable = new LinkedHashMap<>(variantPositions.size(), 10);
                getVariantEntries = (position) -> feature.getNucleotideVariantsAt(position).entrySet();
            } else if (contentMode.equals(Constants.CONTENT_MODE_AMINOACID)) {
                FeatureCoding featureCoding = (FeatureCoding) feature;
                variantPositions = featureCoding.getAminoacidVariantPositions();
                referenceTable = new LinkedHashMap<>(variantPositions.size(), 10);
                getVariantEntries = (position) -> featureCoding.getAminoacidVariantsAt(position).entrySet();
            } else {
                throw new MusialException("Failed to export sequences with parameter --content=" + contentMode + ".");
            }
            for (Integer position : variantPositions) {
                for (Map.Entry<String, VariantInformation> entry : getVariantEntries.apply(position)) {
                    alt = entry.getKey();
                    variantInformation = entry.getValue();
                    if (variantOccurs.apply(variantInformation)) {
                        if (variantInformation.getInfo(Constants.TYPE).contains(Constants.TYPE_INSERTION)) {
                            if (variantInformation.referenceContent.length() > referenceTable.getOrDefault(position, "").length())
                                referenceTable.put(position, variantInformation.referenceContent);
                        } else if (variantInformation.getInfo(Constants.TYPE).contains(Constants.TYPE_DELETION)) {
                            // Insert position wise reference content of deletion that occurs in any sample.
                            chars = variantInformation.referenceContent.toCharArray();
                            for (int i = 0; i < chars.length; i++)
                                referenceTable.putIfAbsent(position + i, String.valueOf(chars[i]));
                        } else if (variantInformation.getInfo(Constants.TYPE).contains(Constants.TYPE_SUBSTITUTION)) {
                            // Insert reference content of substitution that occurs in any sample.
                            referenceTable.putIfAbsent(position, variantInformation.referenceContent);
                        } else {
                            Logger.logWarning("Unable to process variant " + position + " " + alt + " of type " + variantInformation.getInfo(Constants.TYPE) + ".");
                        }
                    }
                }
            }
        }

        // Build output.
        try (BufferedWriter outputWriter = new BufferedWriter(new FileWriter(parameters.getOptionValue("O")), 32768)) {
            TriConsumer<String, String, Collection<String>> writeEntry = (n, k, C) -> {
                try {
                    outputWriter.write(">" + n);
                    outputWriter.newLine();
                    if (Objects.nonNull(k)) {
                        outputWriter.write(";" + k);
                        outputWriter.newLine();
                    }
                    String sequence = String.join("", C);
                    if (!aligned)
                        sequence = sequence.replaceAll(Constants.DELETION_OR_GAP_STRING, "");
                    outputWriter.write(String.join(IO.LINE_SEPARATOR, Splitter.fixedLength(80).split(sequence)));
                    outputWriter.newLine();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }; // Internal helper function to dump content to file.
            BiFunction<String, Integer, String> substringOrEmpty = (s, i) -> {
                if (i > s.length())
                    return "";
                else
                    return s.substring(i);
            }; // Internal helper function to truncate strings.
            writeEntry.accept(Constants.REFERENCE_FORM_NAME, null, referenceTable.values());
            Form form;
            Iterator<String> formNameIterator;
            LinkedHashMap<Integer, String> formVariants;
            int position;
            String alt;
            String type;
            String[] variantFields;
            if (contentMode.equals(Constants.CONTENT_MODE_NUCLEOTIDE)) {
                formNameIterator = feature.getAlleleNameIterator();
                do {
                    form = feature.getAllele(formNameIterator.next());
                    formVariants = SerializationUtils.clone(referenceTable);
                    if (!form.name.equals(Constants.REFERENCE_FORM_NAME)) {
                        for (String variant : form.variants.split(Constants.FIELD_SEPARATOR_2)) {
                            variantFields = variant.split(Constants.FIELD_SEPARATOR_1);
                            position = Integer.parseInt(variantFields[0]);
                            alt = variantFields[1];
                            type = feature.getNucleotideVariant(position, alt).getInfo(Constants.TYPE);
                            if (type.contains(Constants.TYPE_DELETION)) {
                                // Call is deletion.
                                for (int i = 1; i < alt.length(); i++)
                                    formVariants.put(
                                            position + i,
                                            alt.charAt(i) + formVariants.get(position + i).substring(1)
                                    );
                            } else if (type.contains(Constants.TYPE_SUBSTITUTION)) {
                                // Call is substitution.
                                formVariants.put(
                                        position,
                                        alt + formVariants.get(position).substring(1)
                                );
                            } else if (type.contains(Constants.TYPE_INSERTION)) {
                                // Call is insertion.
                                formVariants.put(
                                        position,
                                        formVariants.get(position).charAt(0) + alt.substring(1) + substringOrEmpty.apply(formVariants.get(position), alt.length())
                                );
                            }
                        }
                        if (merge)
                            writeEntry.accept(form.name, null, formVariants.values());
                        else
                            for (String sampleName : form.getOccurrenceSet())
                                writeEntry.accept(sampleName, null, formVariants.values());
                    }
                    outputWriter.flush();
                } while (formNameIterator.hasNext());
            } else {
                FeatureCoding featureCoding = ((FeatureCoding) feature);
                formNameIterator = featureCoding.getProteoformNameIterator();
                do {
                    form = featureCoding.getProteoform(formNameIterator.next());
                    formVariants = SerializationUtils.clone(referenceTable);
                    if (!form.name.equals(Constants.REFERENCE_FORM_NAME)) {
                        for (String variant : form.variants.split(Constants.FIELD_SEPARATOR_2)) {
                            variantFields = variant.split(Constants.FIELD_SEPARATOR_1);
                            position = Integer.parseInt(variantFields[0]);
                            alt = variantFields[1];
                            type = featureCoding.getAminoacidVariant(position, alt).getInfo(Constants.TYPE);
                            if (type.contains(Constants.TYPE_DELETION)) {
                                // Call is deletion.
                                for (int i = 1; i < alt.length(); i++)
                                    formVariants.put(
                                            position + i,
                                            alt.charAt(i) + formVariants.get(position).substring(1)
                                    );
                            } else if (type.contains(Constants.TYPE_SUBSTITUTION)) {
                                // Call is substitution.
                                formVariants.put(
                                        position,
                                        alt + formVariants.get(position).substring(1)
                                );
                            } else if (type.contains(Constants.TYPE_INSERTION)) {
                                // Call is insertion.
                                formVariants.put(
                                        position,
                                        formVariants.get(position).charAt(0) + alt.substring(1) + substringOrEmpty.apply(formVariants.get(position), alt.length())
                                );
                            }
                        }
                        if (merge)
                            writeEntry.accept(form.name, null, formVariants.values());
                        else
                            for (String sampleName : form.getOccurrenceSet())
                                writeEntry.accept(sampleName, null, formVariants.values());
                    }
                    outputWriter.flush();
                } while (formNameIterator.hasNext());
            }
        }
        Logger.logStatus("Done writing output to file `" + new File(parameters.getOptionValue("O")).getAbsolutePath() + "`");
    }

    /**
     * Internal method to output view_x tasks of MUSIAL.
     *
     * @param content   The content for which a tsv format string is to be built.
     * @param noEntries The number of entries to iterate over.
     * @param output    {@link String} representing a file path to write the output to; If set to null, output will be printed to console.
     * @throws MusialException If I/O operation fails in case of a provided output file.
     */
    private static void outputView(LinkedHashMap<String, ArrayList<String>> content, int noEntries, String output) throws MusialException {
        Logger.logStatus("Write output to " + (output == null ? "stdout" : output));
        // Build string content from collected information.
        StringBuilder outputBuilder = new StringBuilder(noEntries);
        ArrayList<String> entry;
        outputBuilder.append(String.join("\t", content.keySet())).append("\n");
        for (int i = 0; i < noEntries; i++) {
            entry = new ArrayList<>();
            for (ArrayList<String> annotationValues : content.values()) {
                entry.add(annotationValues.get(i));
            }
            outputBuilder.append(String.join("\t", entry)).append("\n");
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
     * Internal method to transfers a value to a container that represents column vectors of a table like object to output.
     * <p>
     * This supplements missing values in any object to consider with JSON serializable 'null' entries.
     *
     * @param container The container to add an entry to.
     * @param noEntries The number of existing entries in the container.
     * @param key       The key of the entry.
     * @param value     The value of the entry.
     */
    private static void transferContent(LinkedHashMap<String, ArrayList<String>> container, String key, String value, int noEntries) {
        if (container.containsKey(key)) {
            container.get(key).add(value);
        } else {
            container.put(key, new ArrayList<>(noEntries));
            for (int i = 0; i < noEntries; i++) {
                container.get(key).add("null");
            }
            container.get(key).add(value);
        }
    }

    /**
     * Updates the statistics related to genetic variants and allele frequencies based on the data stored in the provided MusialStorage.
     * This method iterates over features, samples, and proteoforms, calculates statistics such as substitution count, insertion count, deletion count,
     * ambiguous count, and allele frequencies, and updates the information accordingly.
     *
     * @param musialStorage The MusialStorage containing the genetic variant data.
     */
    private static void updateStatistics(MusialStorage musialStorage) {
        Logger.logStatus("Update info statistics");
        float noSamples = musialStorage.getNumberOfSamples();
        Iterator<String> iterator;
        // Update feature statistics.
        Feature feature;
        FeatureCoding featureCoding;
        Sample sample;
        float length;
        iterator = musialStorage.getFeatureNameIterator();
        int noSubstitutions;
        int noInsertions;
        int noDeletions;
        int noAmbiguous;
        String alt;
        String ref;
        String type;
        VariantInformation variantInformation;
        while (iterator.hasNext()) {
            // Feature level statistics.
            feature = musialStorage.getFeature(iterator.next());
            length = Math.abs((feature.end - feature.start) + 1);
            noSubstitutions = 0;
            noInsertions = 0;
            noDeletions = 0;
            noAmbiguous = 0;
            for (Integer p : feature.getNucleotideVariantPositions()) {
                for (Map.Entry<String, VariantInformation> entry : feature.getNucleotideVariantsAt(p).entrySet()) {
                    alt = entry.getKey();
                    variantInformation = entry.getValue();
                    variantInformation.addInfo(
                            Constants.FREQUENCY,
                            DECIMAL_FORMATTER.format((variantInformation.getOccurrenceCount() / noSamples) * 100)
                    );
                    ref = variantInformation.referenceContent;
                    type = variantInformation.getInfo(Constants.TYPE);
                    switch (type) {
                        case (Constants.TYPE_AMBIGUOUS_PREFIX + Constants.TYPE_SUBSTITUTION) -> noAmbiguous += 1;
                        case (Constants.TYPE_AMBIGUOUS_PREFIX + Constants.TYPE_INSERTION) ->
                                noAmbiguous += ref.chars().filter(c -> c == '-').count();
                        case (Constants.TYPE_AMBIGUOUS_PREFIX + Constants.TYPE_DELETION) ->
                                noAmbiguous += alt.chars().filter(c -> c == 'N').count() - 1;
                        case (Constants.TYPE_SUBSTITUTION) -> noSubstitutions += 1;
                        case (Constants.TYPE_INSERTION) -> noInsertions += ref.chars().filter(c -> c == '-').count();
                        case (Constants.TYPE_DELETION) -> noDeletions += alt.chars().filter(c -> c == '-').count();
                    }
                }
            }
            if (feature.isCoding()) {
                featureCoding = (FeatureCoding) feature;
                for (Integer p : featureCoding.getAminoacidVariantPositions())
                    for (VariantInformation vI : featureCoding.getAminoacidVariantsAt(p).values())
                        vI.addInfo(
                                Constants.FREQUENCY,
                                DECIMAL_FORMATTER.format((vI.getOccurrenceCount() / noSamples) * 100)
                        );
            }
            feature.addInfo(Constants.NUMBER_OF_SUBSTITUTIONS, String.valueOf(noSubstitutions));
            feature.addInfo(Constants.NUMBER_OF_INSERTIONS, String.valueOf(noInsertions));
            feature.addInfo(Constants.NUMBER_OF_DELETIONS, String.valueOf(noDeletions));
            feature.addInfo(Constants.NUMBER_OF_AMBIGUOUS, String.valueOf(noAmbiguous));
            feature.addInfo(Constants.VARIABLE_POSITIONS_PERCENTAGE, DECIMAL_FORMATTER.format(
                    ((noSubstitutions + noInsertions + noDeletions) / length) * 100
            ));
            int pos;
            Iterator<String> iterator2 = feature.getAlleleNameIterator();
            Form form;
            String[] formVariants;
            String[] formVariantFields;
            while (iterator2.hasNext()) {
                // Allele level statistics.
                form = feature.getAllele(iterator2.next());
                noSubstitutions = 0;
                noInsertions = 0;
                noDeletions = 0;
                noAmbiguous = 0;
                form.addInfo(
                        Constants.FREQUENCY,
                        DECIMAL_FORMATTER.format((form.getOccurrenceCount() / noSamples) * 100)
                );
                if (!form.name.equals(Constants.REFERENCE_FORM_NAME))
                    try {
                        formVariants = form.variants.split(Constants.FIELD_SEPARATOR_2);
                        for (String formVariant : formVariants) {
                            formVariantFields = formVariant.split(Constants.FIELD_SEPARATOR_1);
                            pos = Integer.parseInt(formVariantFields[0]);
                            alt = formVariantFields[1];
                            variantInformation = feature.getNucleotideVariant(pos, alt);
                            ref = variantInformation.referenceContent;
                            type = variantInformation.getInfo(Constants.TYPE);
                            switch (type) {
                                case (Constants.TYPE_AMBIGUOUS_PREFIX + Constants.TYPE_SUBSTITUTION) -> noAmbiguous += 1;
                                case (Constants.TYPE_AMBIGUOUS_PREFIX + Constants.TYPE_INSERTION) ->
                                        noAmbiguous += ref.chars().filter(c -> c == '-').count();
                                case (Constants.TYPE_AMBIGUOUS_PREFIX + Constants.TYPE_DELETION) ->
                                        noAmbiguous += alt.chars().filter(c -> c == 'N').count() - 1;
                                case (Constants.TYPE_SUBSTITUTION) -> noSubstitutions += 1;
                                case (Constants.TYPE_INSERTION) -> noInsertions += ref.chars().filter(c -> c == '-').count();
                                case (Constants.TYPE_DELETION) -> noDeletions += alt.chars().filter(c -> c == '-').count();
                            }
                        }
                    } catch (NullPointerException e) {
                        e.printStackTrace();
                        System.exit(1);
                    }
                form.addInfo(Constants.NUMBER_OF_SUBSTITUTIONS, String.valueOf(noSubstitutions));
                form.addInfo(Constants.NUMBER_OF_INSERTIONS, String.valueOf(noInsertions));
                form.addInfo(Constants.NUMBER_OF_DELETIONS, String.valueOf(noDeletions));
                form.addInfo(Constants.NUMBER_OF_AMBIGUOUS, String.valueOf(noAmbiguous));
                form.addInfo(Constants.VARIABLE_POSITIONS_PERCENTAGE, DECIMAL_FORMATTER.format(
                        ((noSubstitutions + noInsertions + noDeletions) / length) * 100
                ));
            }
            // Proteoform level statistics.
            if (feature.isCoding()) {
                featureCoding = (FeatureCoding) feature;
                iterator2 = featureCoding.getProteoformNameIterator();
                while (iterator2.hasNext()) {
                    form = featureCoding.getProteoform(iterator2.next());
                    noSubstitutions = 0;
                    noInsertions = 0;
                    noDeletions = 0;
                    noAmbiguous = 0;
                    form.addInfo(
                            Constants.FREQUENCY,
                            DECIMAL_FORMATTER.format((form.getOccurrenceCount() / noSamples) * 100)
                    );
                    if (!form.name.equals(Constants.REFERENCE_FORM_NAME))
                        try {
                            formVariants = form.variants.split(Constants.FIELD_SEPARATOR_2);
                            for (String formVariant : formVariants) {
                                formVariantFields = formVariant.split(Constants.FIELD_SEPARATOR_1);
                                pos = Integer.parseInt(formVariantFields[0]);
                                alt = formVariantFields[1];
                                variantInformation = featureCoding.getAminoacidVariant(pos, alt);
                                ref = variantInformation.referenceContent;
                                type = variantInformation.getInfo(Constants.TYPE);
                                switch (type) {
                                    case (Constants.TYPE_AMBIGUOUS_PREFIX + Constants.TYPE_SUBSTITUTION) -> noAmbiguous += 1;
                                    case (Constants.TYPE_AMBIGUOUS_PREFIX + Constants.TYPE_INSERTION) ->
                                            noAmbiguous += ref.chars().filter(c -> c == '-').count();
                                    case (Constants.TYPE_AMBIGUOUS_PREFIX + Constants.TYPE_DELETION) ->
                                            noAmbiguous += alt.chars().filter(c -> c == 'X').count() - 1;
                                    case (Constants.TYPE_SUBSTITUTION) -> noSubstitutions += 1;
                                    case (Constants.TYPE_INSERTION) -> noInsertions += ref.chars().filter(c -> c == '-').count();
                                    case (Constants.TYPE_DELETION) -> noDeletions += alt.chars().filter(c -> c == '-').count();
                                }
                            }
                        } catch (NullPointerException e) {
                            e.printStackTrace();
                            System.exit(1);
                        }
                    form.addInfo(Constants.NUMBER_OF_SUBSTITUTIONS, String.valueOf(noSubstitutions));
                    form.addInfo(Constants.NUMBER_OF_INSERTIONS, String.valueOf(noInsertions));
                    form.addInfo(Constants.NUMBER_OF_DELETIONS, String.valueOf(noDeletions));
                    form.addInfo(Constants.NUMBER_OF_AMBIGUOUS, String.valueOf(noAmbiguous));
                    form.addInfo(Constants.VARIABLE_POSITIONS_PERCENTAGE, DECIMAL_FORMATTER.format(
                            ((noSubstitutions + noInsertions + noDeletions) / length) * 100
                    ));
                }
            }
        }
        // Sample level statistics.
        iterator = musialStorage.getSampleNameIterator();
        while (iterator.hasNext()) {
            sample = musialStorage.getSample(iterator.next());
            noSubstitutions = 0;
            noInsertions = 0;
            noDeletions = 0;
            noAmbiguous = 0;
            Iterator<String> iterator2 = musialStorage.getFeatureNameIterator();
            while (iterator2.hasNext()) {
                feature = musialStorage.getFeature(iterator2.next());
                noSubstitutions += Integer.parseInt(feature.getAllele(sample.getAllele(feature.name)).getInfo(Constants.NUMBER_OF_SUBSTITUTIONS));
                noInsertions += Integer.parseInt(feature.getAllele(sample.getAllele(feature.name)).getInfo(Constants.NUMBER_OF_INSERTIONS));
                noDeletions += Integer.parseInt(feature.getAllele(sample.getAllele(feature.name)).getInfo(Constants.NUMBER_OF_DELETIONS));
                noAmbiguous += Integer.parseInt(feature.getAllele(sample.getAllele(feature.name)).getInfo(Constants.NUMBER_OF_AMBIGUOUS));
            }
            sample.addInfo(Constants.NUMBER_OF_SUBSTITUTIONS, String.valueOf(noSubstitutions));
            sample.addInfo(Constants.NUMBER_OF_INSERTIONS, String.valueOf(noInsertions));
            sample.addInfo(Constants.NUMBER_OF_DELETIONS, String.valueOf(noDeletions));
            sample.addInfo(Constants.NUMBER_OF_AMBIGUOUS, String.valueOf(noAmbiguous));
        }
        System.gc();
    }
}