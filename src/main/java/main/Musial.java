package main;

import cli.*;
import com.google.gson.internal.LinkedTreeMap;
import components.*;
import datastructure.*;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import runnables.SampleAnalyzerRunnable;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Properties;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import static datastructure.VariantsDictionary.FIELD_SEPARATOR_1;

/**
 * Main class of MUSIAL (MUlti Sample varIant AnaLysis), a tool to calculate SNV, gene, and whole genome alignments,
 * together with other relevant statistics based on vcf files.
 *
 * @author Alexander Seitz
 * @author Simon Hackl
 * @version 2.1
 * @since 1.0
 */
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
     * Specifies the module to execute during runtime.
     */
    public static MusialModules MODULE = null;
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
     * Whether to compress output files.
     */
    public static boolean COMPRESS = false;
    /**
     * Whether to output cli information messages.
     */
    public static boolean SILENT = false;
    /**
     * Whether to print stack trace.
     */
    private static final boolean DEBUG = true;
    /**
     * Project wide formatter to convert decimal numbers to strings.
     */
    public static final DecimalFormat DECIMAL_FORMATTER = new DecimalFormat("#.##");

    /**
     * Main method of MUSIAL; loads meta-information and invokes methods dependent on the run configuration.
     *
     * @param args {@link String[]} comprising the arguments parsed from command line.
     */
    @SuppressWarnings("unchecked")
    public static void main(String[] args) {
        try {
            loadMetadata();
            if (args.length == 0) {
                Logging.logError("No arguments were specified");
                printInfo();
            } else {
                CLIParser cliParser = new CLIParser(args);
                if (cliParser.configuration.size() == 0) {
                    Logging.logError("No modules were specified in " + Logging.colorParameter(cliParser.arguments.getOptionValue("c")));
                    printInfo();
                }
                for (Object o : cliParser.configuration) {
                    MODULE = null;
                    String moduleId = String.valueOf(((LinkedTreeMap<?, ?>) o).get("MODULE"));
                    // Match specified modules with implemented modules.
                    for (MusialModules mm : MusialModules.values()) {
                        if (mm.name().equals(moduleId)) {
                            MODULE = mm;
                            break;
                        }
                    }
                    if (MODULE == null) {
                        Logging.logWarning("Skip unknown module " + Logging.getPurpleTag(moduleId));
                    } else {
                        switch (MODULE) {
                            case BUILD -> executeBUILD(
                                    new ModuleBuildParameters(
                                            (LinkedTreeMap<Object, Object>) o
                                    )
                            );
                            case EXTRACT -> executeEXTRACT(
                                    new ModuleExtractParameters(
                                            (LinkedTreeMap<Object, Object>) o
                                    )
                            );
                        }
                    }
                }
            }
        } catch (Exception e) {
            Logging.logError(e.getMessage());
            if (DEBUG) e.printStackTrace();
            System.exit(-1);
        }
    }

    /**
     * Prints information about available modules to the user and exits.
     */
    public static void printInfo() {
        System.out.println("Please type " + CLIColors.BLACK_BOLD + "java -jar " + Musial.NAME + "-" + Musial.VERSION + ".jar [-h|--help]" + CLIColors.RESET + " for more information.");
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
        Logging.logSoftwareInfo();
    }

    /**
     * Generates a new or updates an existing variants dictionary JSON file based on the specifications parsed from a variants dictionary specification JSON file.
     *
     * @param parameters {@link ModuleBuildParameters} instance yielding parameter specification for the MUSIAL BUILD module.
     * @throws InterruptedException Thrown if any {@link Runnable} implementing class instance is interrupted by the user.
     * @throws IOException          Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialException      If any method fails wrt. biological context, i.e. parsing of unknown symbols; If any method fails wrt. internal logic, i.e. assignment of proteins to genomes; If any input or output file is missing or unable to being generated.
     */
    private static void executeBUILD(ModuleBuildParameters parameters)
            throws InterruptedException, MusialException, IOException {
        Logging.logStatus("Execute module " + Logging.getPurpleTag("BUILD"));

        // Build new empty variants dictionary.
        VariantsDictionary variantsDictionary = VariantsDictionaryFactory.build(parameters);

        // Add feature information.
        Logging.logStatus(Logging.getStartTag() + " Add feature information");
        // Fetch fasta container of feature parent sequence.
        FastaContainer featureContig;
        for (String featureId : parameters.features.keySet()) {
            featureContig = null;
            for (FastaContainer fastaContainer : IO.readFastaToSet(parameters.referenceFASTA)) {
                String contig = fastaContainer.getHeader().split(" ")[0].trim();
                if (parameters.features.get(featureId).chromosome.equals(contig)) {
                    featureContig = fastaContainer;
                    break;
                }
            }
            if (featureContig == null) {
                throw new MusialException("Failed to match contig " + Logging.colorParameter(parameters.features.get(featureId).chromosome) + " from `fasta` file " + Logging.colorParameter(parameters.referenceFASTA.getAbsolutePath()));
            } else {
                parameters.features.get(featureId).imputeNucleotideSequence(featureContig);
                if (parameters.features.get(featureId).isCodingSequence) {
                    parameters.features.get(featureId).imputeProteinInformation();
                }
                variantsDictionary.features.put(featureId, parameters.features.get(featureId));
            }
        }
        // Assign contigs as features, if genome analysis is specified.
        for (FastaContainer fastaContainer : IO.readFastaToSet(parameters.referenceFASTA)) {
            String contig = fastaContainer.getHeader().split(" ")[0].trim();
            if (parameters.excludedPositions.containsKey(contig)) {
                variantsDictionary.excludedPositions.put(contig, parameters.excludedPositions.get(contig));
            }
            if (parameters.genomeAnalysis) {
                variantsDictionary.features.put(contig, new FeatureEntry(
                        contig,
                        contig,
                        0,
                        fastaContainer.getSequence().length(),
                        false
                ));
                variantsDictionary.features.get(contig).imputeNucleotideSequence(fastaContainer);
            }
        }
        Logging.logStatus(Logging.getDoneTag());

        // Collect information about sample variants.
        Logging.logStatus(Logging.getStartTag() + " Add sample information ");
        ExecutorService executor = Executors.newFixedThreadPool(Musial.THREADS);
        for (String sampleId : parameters.samples.keySet()) {
            SampleEntry sampleEntry = parameters.samples.get(sampleId);
            SampleEntry.imputeVCFFileReader(sampleEntry);
            variantsDictionary.samples.put(sampleEntry.name, sampleEntry);
            for (String featureId : variantsDictionary.features.keySet()) {
                executor.execute(
                        new SampleAnalyzerRunnable(sampleEntry, variantsDictionary.features.get(featureId), variantsDictionary)
                );
            }
        }
        executor.shutdown();
        //noinspection ResultOfMethodCallIgnored
        executor.awaitTermination(30, TimeUnit.MINUTES);
        Logging.logStatus(Logging.getDoneTag());

        // Run SnpEff annotation for variants.
        Logging.logStatus(Logging.getStartTag() + " Annotate variants with SnpEff");
        File tmpDirectory = new File("./tmp/");
        try {
            IO.generateDirectory(tmpDirectory);
            File variantsFile = new File("./tmp/variants.vcf");
            IO.generateFile(variantsFile);
            IO.writeVcf(variantsFile, variantsDictionary);
            SnpEffAnnotator.runSnpEff(
                    tmpDirectory,
                    variantsFile,
                    parameters.referenceFASTA,
                    parameters.referenceGFF
            );
            List<String> variantAnnotations = IO.readLinesFromFile("./tmp/annotated_variants.vcf").stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());
            String[] splitAnnotation;
            String position;
            String refContent;
            StringBuilder altContent;
            String[] annotationFields;
            for (String variantAnnotation : variantAnnotations) {
                splitAnnotation = variantAnnotation.split("\t");
                position = splitAnnotation[1];
                refContent = splitAnnotation[3];
                altContent = new StringBuilder(splitAnnotation[4]);
                if (refContent.length() > altContent.length()) {
                    altContent.append(String.valueOf(Bio.DELETION_AA1).repeat(refContent.length() - altContent.length()));
                }
                if (splitAnnotation[7].equals(".")) {
                    continue;
                } else {
                    annotationFields = splitAnnotation[7].replace("ANN=", "").split("\\|");
                }
                for (int i = 1; i < annotationFields.length; i++) {
                    if (i >= 12) {
                        break;
                    }
                    variantsDictionary.nucleotideVariants.get(Integer.valueOf(position)).get(altContent.toString()).annotations.put(
                            NucleotideVariantEntry.PROPERTY_PREFIX_SNP_EF + NucleotideVariantEntry.PROPERTY_LIST_SNP_EFF.get(i - 1),
                            annotationFields[i]
                    );
                }
                altContent.setLength(0);
            }
        } finally {
            IO.deleteDirectory(tmpDirectory);
        }
        Logging.logStatus(Logging.getDoneTag());

        // Infer feature alleles from collected variant information.
        Logging.logStatus(Logging.getStartTag() + " Infer alleles");
        for (FeatureEntry featureEntry : variantsDictionary.features.values()) {
            featureEntry.inferAlleleInformation(variantsDictionary);
            // Compute allele frequency.
            for (AlleleEntry alleleEntry : featureEntry.alleles.values()) {
                alleleEntry.annotations.put(
                        AlleleEntry.PROPERTY_NAME_FREQUENCY,
                        Musial.DECIMAL_FORMATTER.format(100.0 * (alleleEntry.samples.size() / (float) variantsDictionary.samples.size())).replace(",", ".")
                );
            }
        }
        Logging.logStatus(Logging.getDoneTag());

        // Infer coding feature variant/proteoform information; For all coding features.
        Logging.logStatus(Logging.getStartTag() + " Infer proteoforms");
        ConcurrentSkipListMap<String, Tuple<String, String>> aminoacidVariants;
        List<String> featureIdsWithProteinInformation = variantsDictionary.features.keySet().stream().filter(
                feKey -> variantsDictionary.features.get(feKey).isCodingSequence
        ).collect(Collectors.toList());
        for (String featureId : featureIdsWithProteinInformation) {
            for (String sampleId : variantsDictionary.samples.keySet()) {
                aminoacidVariants = Bio.inferAminoacidVariants(variantsDictionary, featureId, sampleId);
                variantsDictionary.features.get(featureId).addProteoform(sampleId, aminoacidVariants, variantsDictionary);
            }
            for (ProteoformEntry proteoformEntry : variantsDictionary.features.get(featureId).proteoforms.values()) {
                proteoformEntry.annotations.put(
                        ProteoformEntry.PROPERTY_NAME_FREQUENCY,
                        Musial.DECIMAL_FORMATTER.format(100.0 * (proteoformEntry.samples.size() / (float) variantsDictionary.samples.size())).replace(",", ".")
                );
            }
        }
        Logging.logStatus(Logging.getDoneTag());

        // Infer variant frequencies.
        Logging.logStatus(Logging.getStartTag() + " Compute statistics");
        float noSamples = variantsDictionary.samples.size();
        float noOccurrences;
        // Compute frequencies of nucleotide variants.
        for (ConcurrentSkipListMap<String, NucleotideVariantEntry> perPositionNucleotideVariants : variantsDictionary.nucleotideVariants.values()) {
            for (NucleotideVariantEntry nucleotideVariant : perPositionNucleotideVariants.values()) {
                noOccurrences = 0;
                for (String occurrence : nucleotideVariant.occurrence) {
                    String featureId = occurrence.split(FIELD_SEPARATOR_1)[0];
                    String alleleId = occurrence.split(FIELD_SEPARATOR_1)[1];
                    noOccurrences += variantsDictionary.features.get(featureId).alleles.get(alleleId).samples.size();
                }
                nucleotideVariant.annotations.put(
                        NucleotideVariantEntry.PROPERTY_NAME_FREQUENCY,
                        Musial.DECIMAL_FORMATTER.format(100L * (noOccurrences / noSamples)).replace(",", ".")
                );
            }
        }
        // Compute frequencies of aminoacid variants.
        for (FeatureEntry feature : variantsDictionary.features.values()) {
            for (ConcurrentSkipListMap<String, AminoacidVariantEntry> perPositionAminoacidVariants : feature.aminoacidVariants.values()) {
                for (AminoacidVariantEntry aminoacidVariant : perPositionAminoacidVariants.values()) {
                    noOccurrences = 0;
                    for (String occurrence : aminoacidVariant.occurrence) {
                        noOccurrences += feature.proteoforms.get(occurrence).samples.size();
                    }
                    aminoacidVariant.annotations.put(
                            AminoacidVariantEntry.PROPERTY_NAME_FREQUENCY,
                            Musial.DECIMAL_FORMATTER.format(100 * (noOccurrences / noSamples)).replace(",", ".")
                    );
                }
            }
        }
        Logging.logStatus(Logging.getDoneTag());

        // Write updated database to file.
        Logging.logStatus(Logging.getStartTag() + " Dump to file");
        String outputFile = parameters.outputFile.getAbsolutePath();
        File vDictOutfile = new File(outputFile);
        if (!vDictOutfile.exists()) {
            IO.generateFile(vDictOutfile);
        }
        variantsDictionary.dump(vDictOutfile);
        Logging.logStatus(Logging.getDoneTag());
    }

    /**
     * TODO
     *
     * @param parameters {@link ModuleExtractParameters} instance yielding parameter specification for MUSIAL EXTRACT modules.
     * @throws IOException     Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialException If any method fails wrt. biological context or internal logic, e.g. only samples were specified that are not stored in the specified variants dictionary.
     */
    @SuppressWarnings("DuplicatedCode")
    private static void executeEXTRACT(ModuleExtractParameters parameters) throws IOException, MusialException {
        Logging.logStatus(
                "Execute module "
                        + Logging.getPurpleTag("EXTRACT")
                        + " " + Logging.getPurpleTag(String.valueOf(parameters.contentMode))
                        + " " + Logging.getPurpleTag(String.valueOf(parameters.outputMode))
                        + " " + Logging.getPurpleTag("Excl. indels: " + parameters.excludeIndels)
                        + " " + Logging.getPurpleTag("Incl. non-variant positions: " + parameters.excludeConservedPositions)
                        + " " + Logging.getPurpleTag("Group samples: " + parameters.grouped)
        );

        // Load variants dictionary.
        VariantsDictionary variantsDictionary = VariantsDictionaryFactory.load(parameters.inputFile);
        if (parameters.samples.size() == 0) {
            parameters.samples = new ArrayList<>(variantsDictionary.samples.keySet());
        } else {
            for (String sampleIdentifier : parameters.samples) {
                if (!variantsDictionary.samples.containsKey(sampleIdentifier)) {
                    parameters.samples.remove(sampleIdentifier);
                    Logging.logWarning("Removed sample " + sampleIdentifier + " from analysis: No data is stored under the key.");
                }
            }
        }
        if (parameters.features.size() == 0) {
            parameters.features = new ArrayList<>(variantsDictionary.features.keySet());
        } else {
            for (String featureIdentifer : parameters.features) {
                if (!variantsDictionary.features.containsKey(featureIdentifer)) {
                    parameters.features.remove(featureIdentifer);
                    Logging.logWarning("Removed sample " + featureIdentifer + " from analysis: No data is stored under the key.");
                }
            }
        }
        // Process each feature.
        for (String featureIdentifier : parameters.features) {
            switch (parameters.outputMode) {
                case TABLE -> {
                    if (parameters.contentMode.equals(ModuleExtractContentModes.AMINOACID)
                            && !variantsDictionary.features.get(featureIdentifier).isCodingSequence) {
                        Logging.logWarning(
                                "Skipping non-coding feature "
                                        + featureIdentifier
                                        + " (incompatible with mode "
                                        + parameters.contentMode
                                        + ")."
                        );
                        continue;
                    }
                    Logging.logStatus(
                            Logging.getStartTag() +
                                    " Construct variants table ("
                                    + parameters.contentMode
                                    + ") of feature "
                                    + featureIdentifier
                                    + " and "
                                    + parameters.samples.size()
                                    + " samples"
                    );
                    VariantsTable variantsTable = new VariantsTable(
                            variantsDictionary,
                            parameters.samples,
                            featureIdentifier,
                            parameters.contentMode,
                            parameters.excludeIndels,
                            !parameters.excludeConservedPositions,
                            parameters.keepOnlyVariantsWith,
                            false
                    );
                    Logging.logStatus(Logging.getDoneTag());
                    Logging.logStatus(Logging.getStartTag() + " Write table output");
                    variantsTable.writeTableToFile(parameters.outputDirectory + "/" + featureIdentifier + ".tsv", parameters.grouped);
                    if (Musial.COMPRESS) {
                        IO.compressFile(new File(parameters.outputDirectory + "/" + featureIdentifier + ".tsv"));
                    }
                    Logging.logStatus(Logging.getDoneTag());
                }
                case SEQUENCE_ALIGNED -> {
                    if (parameters.contentMode.equals(ModuleExtractContentModes.AMINOACID)
                            && !variantsDictionary.features.get(featureIdentifier).isCodingSequence) {
                        Logging.logWarning(
                                "Skipping non-coding feature "
                                        + featureIdentifier
                                        + " (incompatible with mode "
                                        + parameters.contentMode
                                        + ")."
                        );
                        continue;
                    }
                    Logging.logStatus(
                            Logging.getStartTag() +
                                    " Construct variants table ("
                                    + parameters.contentMode
                                    + ") of feature "
                                    + featureIdentifier
                                    + " and "
                                    + parameters.samples.size()
                                    + " samples"
                    );
                    VariantsTable variantsTable = new VariantsTable(
                            variantsDictionary,
                            parameters.samples,
                            featureIdentifier,
                            parameters.contentMode,
                            parameters.excludeIndels,
                            !parameters.excludeConservedPositions,
                            parameters.keepOnlyVariantsWith,
                            true
                    );
                    Logging.logStatus(Logging.getDoneTag());
                    // Generate fasta entries.
                    Logging.logStatus(Logging.getStartTag() + " Write fasta output");
                    ArrayList<Tuple<String, String>> fastaEntries = new ArrayList<>();
                    if (parameters.grouped) {
                        HashSet<String> accessors = new HashSet<>();
                        for (String sample : parameters.samples) {
                            switch (parameters.contentMode) {
                                case AMINOACID -> accessors.add(
                                        variantsDictionary.samples.get(sample).annotations.get("PF" + FIELD_SEPARATOR_1 + featureIdentifier)
                                );
                                case NUCLEOTIDE -> accessors.add(
                                        variantsDictionary.samples.get(sample).annotations.get("AL" + FIELD_SEPARATOR_1 + featureIdentifier)
                                );
                            }
                        }
                        for (String accessor : accessors) {
                            fastaEntries.add(variantsTable.getFastaEntry(accessor, false));
                        }
                    } else {
                        for (String sample : parameters.samples) {
                            fastaEntries.add(variantsTable.getFastaEntry(sample, false));
                        }
                    }
                    switch (parameters.contentMode) {
                        case AMINOACID -> fastaEntries.add(
                                variantsTable.getFastaEntry(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID, false)
                        );
                        case NUCLEOTIDE -> fastaEntries.add(
                                variantsTable.getFastaEntry(AlleleEntry.PROPERTY_NAME_REFERENCE_ID, false)
                        );
                    }
                    IO.writeFasta(
                            new File(parameters.outputDirectory + "/" + featureIdentifier + ".fasta"),
                            fastaEntries
                    );
                    if (Musial.COMPRESS) {
                        IO.compressFile(new File(parameters.outputDirectory + "/" + featureIdentifier + ".fasta"));
                    }
                    Logging.logStatus(Logging.getDoneTag());
                }
                case SEQUENCE -> {
                    if (parameters.contentMode.equals(ModuleExtractContentModes.AMINOACID)
                            && !variantsDictionary.features.get(featureIdentifier).isCodingSequence) {
                        Logging.logWarning(
                                "Skipping non-coding feature "
                                        + featureIdentifier
                                        + " (incompatible with mode "
                                        + parameters.contentMode
                                        + ")."
                        );
                        continue;
                    }
                    Logging.logStatus(
                            Logging.getStartTag() +
                                    " Construct variants table ("
                                    + parameters.contentMode
                                    + ") of feature "
                                    + featureIdentifier
                                    + " and "
                                    + parameters.samples.size()
                                    + " samples"
                    );
                    VariantsTable variantsTable = new VariantsTable(
                            variantsDictionary,
                            parameters.samples,
                            featureIdentifier,
                            parameters.contentMode,
                            parameters.excludeIndels,
                            !parameters.excludeConservedPositions,
                            parameters.keepOnlyVariantsWith,
                            true
                    );
                    Logging.logStatus(Logging.getDoneTag());
                    // Generate fasta entries.
                    Logging.logStatus(Logging.getStartTag() + " Write fasta output");
                    ArrayList<Tuple<String, String>> fastaEntries = new ArrayList<>();
                    if (parameters.grouped) {
                        HashSet<String> accessors = new HashSet<>();
                        for (String sample : parameters.samples) {
                            switch (parameters.contentMode) {
                                case AMINOACID -> accessors.add(
                                        variantsDictionary.samples.get(sample).annotations.get("PF" + FIELD_SEPARATOR_1 + featureIdentifier)
                                );
                                case NUCLEOTIDE -> accessors.add(
                                        variantsDictionary.samples.get(sample).annotations.get("AL" + FIELD_SEPARATOR_1 + featureIdentifier)
                                );
                            }
                        }
                        for (String accessor : accessors) {
                            fastaEntries.add(variantsTable.getFastaEntry(accessor, true));
                        }
                    } else {
                        for (String sample : parameters.samples) {
                            fastaEntries.add(variantsTable.getFastaEntry(sample, true));
                        }
                    }
                    IO.writeFasta(
                            new File(parameters.outputDirectory + "/" + featureIdentifier + ".fasta"),
                            fastaEntries
                    );
                    if (Musial.COMPRESS) {
                        IO.compressFile(new File(parameters.outputDirectory + "/" + featureIdentifier + ".fasta"));
                    }
                    Logging.logStatus(Logging.getDoneTag());
                }
            }
        }
    }
}