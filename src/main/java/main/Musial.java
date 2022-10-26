package main;

import cli.CLIColors;
import cli.CLIParser;
import cli.ModuleBuildParameters;
import cli.ModuleExtractParameters;
import com.aayushatharva.brotli4j.Brotli4jLoader;
import com.aayushatharva.brotli4j.encoder.Encoder;
import com.google.gson.internal.LinkedTreeMap;
import com.sun.source.tree.Tree;
import components.*;
import datastructure.*;
import exceptions.MusialException;
import runnables.SampleAnalyzerRunnable;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.*;
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
                            case EXTRACT_NUC_VARIANTS_TABLE -> executeEXTRACT_NUC_VARIANTS_TABLE(
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
            e.printStackTrace();
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
     * @param parameters {@link ModuleBuildParameters} instance yielding parameter specification for the MUSIAL update variants dictionary module.
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
        Set<String> specifiedFeatures = parameters.features.keySet();
        // Fetch fasta container of feature parent sequence.
        FastaContainer featureContig;
        for (String featureId : specifiedFeatures) {
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
                if (parameters.features.get(featureId).pdbFile != null) {
                    parameters.features.get(featureId).imputeProteinInformation();
                }
                variantsDictionary.features.put(featureId, parameters.features.get(featureId));
            }
        }
        Logging.logStatus(Logging.getDoneTag() + " Add feature information");

        // Collect information about sample variants.
        Logging.logStatus(Logging.getStartTag() + " Add sample information ");
        Set<String> sampleIdsUpdate = parameters.samples.keySet();
        HashSet<String> sampleIdsAll = new HashSet<>();
        sampleIdsAll.addAll(sampleIdsUpdate);
        sampleIdsAll.addAll(variantsDictionary.samples.keySet());
        ExecutorService executor = Executors.newFixedThreadPool(Musial.THREADS);
        for (String sampleId : sampleIdsAll) {
            SampleEntry sampleEntry = parameters.samples.get(sampleId);
            SampleEntry.imputeVCFFileReader(sampleEntry);
            variantsDictionary.samples.put(sampleEntry.name, sampleEntry);
            for (String featureId : specifiedFeatures) {
                executor.execute(
                        new SampleAnalyzerRunnable(sampleEntry, parameters.features.get(featureId), variantsDictionary)
                );
            }
        }
        executor.shutdown();
        //noinspection ResultOfMethodCallIgnored
        executor.awaitTermination(30, TimeUnit.MINUTES);
        Logging.logStatus(Logging.getDoneTag() + " Add sample information");

        // Run SnpEff annotation for variants.
        Logging.logStatus(Logging.getStartTag() + " Annotate novel variants with SnpEff");
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
                variantsDictionary.nucleotideVariants.get(Integer.valueOf(position)).get(altContent.toString()).annotations.put(
                        NucleotideVariantEntry.PROPERTY_NAME_SNP_EFF_TYPE,
                        annotationFields[1]
                );
                variantsDictionary.nucleotideVariants.get(Integer.valueOf(position)).get(altContent.toString()).annotations.put(
                        NucleotideVariantEntry.PROPERTY_NAME_SNP_EFF_IMPACT,
                        annotationFields[2]
                );
                altContent.setLength(0);
            }
        } finally {
            IO.deleteDirectory(tmpDirectory);
        }
        Logging.logStatus(Logging.getDoneTag() + " Annotate novel variants with SnpEff");

        // Infer feature alleles from collected variant information.
        Logging.logStatus(Logging.getStartTag() + " Infer alleles");
        for (FeatureEntry featureEntry : variantsDictionary.features.values()) {
            featureEntry.inferAlleleInformation(variantsDictionary);
        }
        Logging.logStatus(Logging.getDoneTag() + " Infer alleles");

        // Infer coding feature variant/proteoform information; For all coding features.
        Logging.logStatus(Logging.getStartTag() + " Infer protein variants");
        ConcurrentSkipListMap<String, String> aminoacidVariants;
        List<String> featureIdsWithProteinInformation = variantsDictionary.features.keySet().stream().filter(
                feKey -> variantsDictionary.features.get(feKey).isCodingSequence
        ).collect(Collectors.toList());
        for (String featureId : featureIdsWithProteinInformation) {
            for (String sampleId : variantsDictionary.samples.keySet()) {
                aminoacidVariants = Bio.inferAminoacidVariants(variantsDictionary, featureId, sampleId);
                variantsDictionary.features.get(featureId).addProteoform(sampleId, aminoacidVariants, variantsDictionary);
            }
        }
        Logging.logStatus(Logging.getDoneTag() + " Infer protein variants");

        // Infer variant frequencies.
        Logging.logStatus(Logging.getStartTag() + " Infer variant statistics");
        /** FIXME: Currently unused code snippet.
         * int totalNoVariants = 0;
         * for (ConcurrentSkipListMap<String, NucleotideVariantEntry> perPositionNucleotideVariants : variantsDictionary.nucleotideVariants.values()) {
         *     totalNoVariants += perPositionNucleotideVariants.size();
         * }
         * for (FeatureEntry feature : variantsDictionary.features.values()) {
         *     totalNoVariants += feature.aminoacidVariants.size();
         * }
         */
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
        Logging.logStatus(Logging.getDoneTag() + " Infer variant statistics");

        // Write updated database to file.
        Logging.logStatus(Logging.getStartTag() + " Dump to file");
        String outputFile = parameters.outputFile.getAbsolutePath();
        File vDictOutfile = new File(outputFile);
        if (!vDictOutfile.exists()) {
            IO.generateFile(vDictOutfile);
        }
        variantsDictionary.dump(vDictOutfile);
        Logging.logStatus(Logging.getDoneTag() + " Dump to file");
    }

    /**
     *
     */
    private static void executeEXTRACT_NUC_VARIANTS_TABLE(ModuleExtractParameters parameters) throws IOException, MusialException {
        Logging.logStatus("Execute module " + Logging.getPurpleTag("EXTRACT_NUC_VARIANTS_TABLE"));

        // Load variants dictionary.
        VariantsDictionary variantsDictionary = VariantsDictionaryFactory.load(parameters.inputFile);

        // Validate parameters, i.e., ensure that specified samples and features are actually stored.
        Logging.logStatus(Logging.getStartTag() + " Validate parameters");
        List<String> sampleIds;
        List<String> featureIds;
        // Validate sample information.
        if (parameters.samples.size() != 0) {
            sampleIds = parameters.samples.stream().filter(variantsDictionary.samples::containsKey).collect(Collectors.toList());
        } else {
            sampleIds = new ArrayList<>(variantsDictionary.samples.keySet());
        }
        if (sampleIds.size() == 0) {
            throw new MusialException("The number of specified valid samples is zero; " + String.join(" ,", sampleIds));
        }
        // Validate feature information.
        if (parameters.features.size() != 0) {
            featureIds = parameters.features.stream().filter(variantsDictionary.features::containsKey).collect(Collectors.toList());
        } else {
            featureIds = new ArrayList<>(variantsDictionary.features.keySet());
        }
        if (featureIds.size() == 0) {
            throw new MusialException("The number of specified valid features is zero; " + String.join(" ,", featureIds));
        }
        Logging.logStatus(Logging.getDoneTag() + " Validate parameters");

        // Collect variant information.
        Logging.logStatus(Logging.getStartTag() + " Write variant table per feature");
        // HashMap mapping position -> allele -> nuc. content.
        TreeMap<String, HashMap<String, String>> mergedVariants;
        // HashMap containing the variants wrt. one sample/allele and one feature.
        HashMap<Integer, String> singularVariants;
        HashSet<String> processedAlleles;
        String sampleAlleleId;
        int variantStart;
        String variantPosition;
        char[] variantContents;
        char[] referenceContents;
        String variantContent;
        boolean isInsertion;
        boolean isDeletion;
        Writer outputWriter;
        StringBuilder outputContentBuilder;
        for (String featureId : featureIds) {
            // Collect variant information for feature.
            mergedVariants = new TreeMap<>((s1, s2) -> {
                int p1 = Integer.parseInt(s1.split("\\+")[0]);
                int p2 = Integer.parseInt(s2.split("\\+")[0]);
                if (p1 != p2) {
                    return Integer.compare(p1, p2);
                } else {
                    int i1 = s1.contains("\\+") ? Integer.parseInt(s1.split("\\+")[1]) : 0;
                    int i2 = s2.contains("\\+") ? Integer.parseInt(s2.split("\\+")[1]) : 0;
                    return Integer.compare(i1, i2);
                }
            });
            processedAlleles = new HashSet<>();
            for (String sampleId : sampleIds) {
                sampleAlleleId = variantsDictionary.samples.get(sampleId).annotations.get("AL" + FIELD_SEPARATOR_1 + featureId);
                if (sampleAlleleId.equals(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)) {
                    continue;
                }
                if (!processedAlleles.contains(sampleAlleleId)) {
                    singularVariants = variantsDictionary.getNucleotideVariants(featureId, sampleId, false);
                    for (Map.Entry<Integer, String> variantEntry : singularVariants.entrySet()) {
                        variantStart = variantEntry.getKey();
                        variantContents = variantEntry.getValue().toCharArray();
                        referenceContents = variantsDictionary.nucleotideVariants.get(variantStart).get(variantEntry.getValue()).annotations.get(NucleotideVariantEntry.PROPERTY_NAME_REFERENCE_CONTENT).toCharArray();
                        isInsertion = variantContents.length > referenceContents.length;
                        isDeletion = variantEntry.getValue().contains(String.valueOf(Bio.DELETION_AA1));
                        if ((isInsertion || isDeletion) && parameters.SNVOnly) {
                            continue;
                        }
                        if (isInsertion) {
                            referenceContents = variantContents;
                            for (int i = 1; i < referenceContents.length; i++) {
                                referenceContents[i] = Bio.DELETION_AA1;
                            }
                        }
                        for (int i = 0; i < variantContents.length; i++) {
                            variantContent = String.valueOf(variantContents[i]);
                            if (isInsertion && i > 0) {
                                variantPosition = variantStart + "+" + i;
                            } else {
                                variantPosition = String.valueOf(variantStart + i);
                            }
                            if (!mergedVariants.containsKey(variantPosition)) {
                                mergedVariants.put(variantPosition, new HashMap<>());
                            }
                            if (!mergedVariants.get(variantPosition).containsKey("Reference")) {
                                mergedVariants.get(variantPosition).put("Reference", String.valueOf(referenceContents[i]));
                            }
                            mergedVariants.get(variantPosition).put(sampleAlleleId, variantContent);
                        }
                    }
                    processedAlleles.add(sampleAlleleId);
                }
            }
            // Fill in reference content for non-variant positions per allele.
            for (String processedVariantPosition : mergedVariants.keySet()) {
                for (String processedAllele : processedAlleles) {
                    if (!mergedVariants.get(processedVariantPosition).containsKey(processedAllele)) {
                        mergedVariants.get(processedVariantPosition).put(processedAllele, ".");
                    }
                }
            }
            // Construct output string.
            // Write info and header lines.
            outputContentBuilder = new StringBuilder();
            outputContentBuilder
                    .append("#")
                    .append(NAME)
                    .append(VERSION)
                    .append(";")
                    .append(Logging.getTimestampDayTime())
                    .append(";")
                    .append("NUC_VARIANTS_TABLE")
                    .append(IO.LINE_SEPARATOR);
            outputContentBuilder
                    .append("#")
                    .append("FEATURE=")
                    .append(featureId)
                    .append(";START=")
                    .append(variantsDictionary.features.get(featureId).start)
                    .append(";END=")
                    .append(variantsDictionary.features.get(featureId).end)
                    .append(";CHR=")
                    .append(variantsDictionary.features.get(featureId).chromosome)
                    .append(";SENSE=")
                    .append(variantsDictionary.features.get(featureId).isSense)
                    .append(IO.LINE_SEPARATOR);
            outputContentBuilder
                    .append("Position\tReference");
            if (variantsDictionary.features.get(featureId).alleles.containsKey(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)) {
                for (String sampleId : variantsDictionary.features.get(featureId).alleles.get(AlleleEntry.PROPERTY_NAME_REFERENCE_ID).samples) {
                    outputContentBuilder.append("\t").append(sampleId);
                }
            }
            for (String processedAllele : processedAlleles) {
                for (String sampleId : variantsDictionary.features.get(featureId).alleles.get(processedAllele).samples) {
                    outputContentBuilder.append("\t").append(sampleId);
                }
            }
            outputContentBuilder.append(IO.LINE_SEPARATOR);
            for (Map.Entry<String, HashMap<String, String>> processedVariant : mergedVariants.entrySet()) {
                outputContentBuilder
                        .append(processedVariant.getKey())
                        .append("\t")
                        .append(processedVariant.getValue().get("Reference"));
                if (variantsDictionary.features.get(featureId).alleles.containsKey(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)) {
                    for (String ignored : variantsDictionary.features.get(featureId).alleles.get(AlleleEntry.PROPERTY_NAME_REFERENCE_ID).samples) {
                        outputContentBuilder
                                .append("\t")
                                .append(".");
                    }
                }
                for (String processedAllele : processedAlleles) {
                    for (String ignored : variantsDictionary.features.get(featureId).alleles.get(processedAllele).samples) {
                        outputContentBuilder
                                .append("\t")
                                .append(processedVariant.getValue().get(processedAllele));
                    }
                }
                outputContentBuilder.append(IO.LINE_SEPARATOR);
            }
            //Write built content string to file.
            if (Musial.COMPRESS) {
                // Write compressed (brotli) output.
                Brotli4jLoader.ensureAvailability();
                byte[] compressed = Encoder.compress(outputContentBuilder.toString().getBytes());
                Files.write(Paths.get(parameters.outputDirectory + "/" + featureId + ".tsv.br"), compressed);
            } else {
                // Write plain output.
                outputWriter = new BufferedWriter(
                        new FileWriter(parameters.outputDirectory + "/" + featureId + ".tsv")
                );
                outputWriter.write(outputContentBuilder.toString());
                outputWriter.close();
            }
        }
        Logging.logStatus(Logging.getDoneTag() + " Write variant table per feature");
    }

    /**
     * Extract genotype and/or proteoform sequences as alignment or unaligned in .fasta format from an existing variants dictionary.
     * <p>
     * Sequences are written in coding direction.
     *
     * @param cliarguments {@link CLIParametersInferSequences} instance yielding parameter specification for the MUSIAL infer sequences module.
     * @throws MusialIOException  If the generation of output directories or files fails.
     * @throws MusialBioException If sequence extraction procedures from the specified {@link VariantsDictionary} instance fail.
     *
    private static void run (CLIParametersInferSequences cliarguments) throws MusialIOException, MusialBioException {
    // TODO: Generalize methods; Currently code is heavily duplicated for debugging.
    HashMap<String, ArrayList<String>> sequences; // Stores sequences as keys pointing to list of proteoforms/samples.
    String sequence;
    String VSWAB;
    // Write genotype sequences.
    if (cliarguments.GS) {
    Logging.logStatus("Write genotype sequences.");
    try (ProgressBar pb = buildProgress()) {
    pb.maxHint(cliarguments.inputVDict.features.size());
    int gtIndex;
    for (FeatureEntry featureEntry : cliarguments.inputVDict.features.values()) {
    sequences = new HashMap<>();
    // Access and store wild-type nucleotide sequence.
    String featureWildTypeNucleotideSequence = featureEntry.nucleotideSequence;
    sequences.put(featureWildTypeNucleotideSequence, new ArrayList<>());
    sequences.get(featureWildTypeNucleotideSequence).add(VariantsDictionary.WILD_TYPE_SAMPLE_ID);
    for (SampleEntry sampleEntry : cliarguments.inputVDict.samples.values()) {
    VSWAB = sampleEntry.annotations.get(featureEntry.name + VariantsDictionary.ATTRIBUTE_VARIANT_SWAB_NAME);
    if (VSWAB == null) {
    sequences.get(featureWildTypeNucleotideSequence).add(sampleEntry.name);
    } else {
    sequence = cliarguments.inputVDict.getSampleNucleotideSequence(featureEntry.name, sampleEntry.name);
    if (sequences.containsKey(sequence)) {
    sequences.get(sequence).add(sampleEntry.name);
    } else {
    sequences.put(sequence, new ArrayList<>());
    sequences.get(sequence).add(sampleEntry.name);
    }
    }
    }
    gtIndex = 1;
    for (String s : sequences.keySet()) {
    sequences.get(s).add(0, "GT" + gtIndex);
    gtIndex += 1;
    }
    File gtsOutDir = new File(cliarguments.outputDirectory.getAbsolutePath() + "/GenotypeSequences/");
    if (!Validation.isDirectory(gtsOutDir)) {
    IO.generateDirectory(gtsOutDir);
    }
    IO.writeFasta(
    new File(cliarguments.outputDirectory.getAbsolutePath() + "/GenotypeSequences/" + featureEntry.chromosome + "_" + featureEntry.name + "_" + cliarguments.samples.size() + "_genotypeSequences.fasta"),
    sequences
    );
    pb.step();
    }
    }
    }
    // Write proteoform sequences.
    if (cliarguments.PS) {
    Logging.logStatus("Write proteoform sequences.");
    try (ProgressBar pb = buildProgress()) {
    pb.maxHint(cliarguments.inputVDict.features.size());
    for (FeatureEntry featureEntry : cliarguments.inputVDict.features.values()) {
    sequences = new HashMap<>();
    for (String pfId : cliarguments.inputVDict.features.get(featureEntry.name).allocatedProtein.proteoforms.keySet()) {
    sequence = cliarguments.inputVDict.getProteoformSequence(featureEntry.name, pfId);
    sequences.put(sequence, new ArrayList<>());
    sequences.get(sequence).add(pfId);
    sequences.get(sequence).add(
    "VP=" + featureEntry.allocatedProtein.proteoforms.get(pfId).annotations.get("VP")
    + ";PT=" + featureEntry.allocatedProtein.proteoforms.get(pfId).annotations.get("PT")
    );
    sequences.get(sequence).addAll(featureEntry.allocatedProtein.proteoforms.get(pfId).samples);
    }
    File pfsOutDir = new File(cliarguments.outputDirectory.getAbsolutePath() + "/ProteoformSequences/");
    if (!Validation.isDirectory(pfsOutDir)) {
    IO.generateDirectory(pfsOutDir);
    }
    IO.writeFasta(
    new File(cliarguments.outputDirectory.getAbsolutePath() + "/ProteoformSequences/" + featureEntry.chromosome + "_" + featureEntry.name + "_" + cliarguments.samples.size() + "_proteoformSequences.fasta"),
    sequences
    );
    pb.step();
    }
    }
    }
    // Write genotype sequences MSA.
    ConcurrentSkipListMap<String, HashMap<String, Character>> perPositionContents =
    new ConcurrentSkipListMap<>((k1, k2) -> {
    int p1 = Integer.parseInt(k1.split("\\+")[0]);
    int p2 = Integer.parseInt(k2.split("\\+")[0]);
    if (p1 == p2) {
    p1 = Integer.parseInt(k1.split("\\+")[1]);
    p2 = Integer.parseInt(k2.split("\\+")[1]);
    }
    return Integer.compare(p1, p2);
    });
    char[] sequenceChars;
    String positionContent;
    String positionStr;
    int positionInt;
    int entryIndex;
    int VSWABHashCode;
    HashMap<String, ArrayList<String>> faSequences;
    ArrayList<String> faHeadersList;
    StringBuilder faSequenceBuilder = new StringBuilder();
    if (cliarguments.GSMSA) {
    Logging.logStatus("Write genotype sequences alignment.");
    try (ProgressBar pb = buildProgress()) {
    pb.maxHint(cliarguments.inputVDict.features.size());
    HashMap<Integer, String> variants;
    HashMap<Integer, String> VSWABHashToGTIdentifier = new HashMap<>();
    HashMap<String, TreeSet<String>> GTIdentifierToSamples = new HashMap<>();
    for (FeatureEntry featureEntry : cliarguments.inputVDict.features.values()) {
    perPositionContents.clear();
    VSWABHashToGTIdentifier.clear();
    GTIdentifierToSamples.clear();
    entryIndex = 1;
    // Extract reference sequence information.
    sequenceChars = featureEntry.isSense ? featureEntry.nucleotideSequence.toCharArray() : Bio.reverseComplement(featureEntry.nucleotideSequence).toCharArray();
    for (int i = 0; i < sequenceChars.length; i++) {
    if (!perPositionContents.containsKey((i + 1) + "+0")) {
    perPositionContents.put((i + 1) + "+0", new HashMap<>());
    }
    perPositionContents.get((i + 1) + "+0").put("GT0", sequenceChars[i]);
    }
    VSWABHashToGTIdentifier.put(0, "GT0");
    GTIdentifierToSamples.put("GT0", new TreeSet<>(cliarguments.inputVDict.samples.keySet()));
    // Infer per sample variant content information.
    for (SampleEntry sampleEntry : cliarguments.inputVDict.samples.values()) {
    if (sampleEntry.annotations.containsKey(featureEntry.name + VariantsDictionary.ATTRIBUTE_VARIANT_SWAB_NAME)) {
    VSWAB = sampleEntry.annotations.get(featureEntry.name + VariantsDictionary.ATTRIBUTE_VARIANT_SWAB_NAME);
    VSWABHashCode = VSWAB.hashCode();
    if (!VSWABHashToGTIdentifier.containsKey(VSWABHashCode)) {
    VSWABHashToGTIdentifier.put(VSWABHashCode, "GT" + entryIndex);
    GTIdentifierToSamples.put("GT" + entryIndex, new TreeSet<>());
    // Add variants from VSWAB:
    variants = cliarguments.inputVDict.getSampleNucleotideVariants(featureEntry.name, sampleEntry.name);
    for (Map.Entry<Integer, String> variantEntry : variants.entrySet()) {
    positionInt = variantEntry.getKey();
    positionContent = variantEntry.getValue();
    sequenceChars = positionContent.toCharArray();
    if (positionContent.contains(String.valueOf(Bio.DELETION_AA1))) {
    for (int i = 0; i < sequenceChars.length; i++) {
    if (!perPositionContents.containsKey(positionInt + i + "+0")) {
    perPositionContents.put(positionInt + i + "+0", new HashMap<>());
    }
    perPositionContents.get(positionInt + i + "+0").put("GT" + entryIndex, sequenceChars[i]);
    }
    } else {
    for (int i = 0; i < sequenceChars.length; i++) {
    if (!perPositionContents.containsKey(positionInt + "+" + i)) {
    perPositionContents.put(positionInt + "+" + i, new HashMap<>());
    }
    perPositionContents.get(positionInt + "+" + i).put("GT" + entryIndex, sequenceChars[i]);
    }
    }
    }
    entryIndex += 1;
    }
    GTIdentifierToSamples.get(VSWABHashToGTIdentifier.get(VSWABHashCode)).add(sampleEntry.name);
    GTIdentifierToSamples.get("GT0").remove(sampleEntry.name);
    }
    }
    faSequences = new HashMap<>();
    for (String gtIdentifier : GTIdentifierToSamples.keySet()) {
    faHeadersList = new ArrayList<>();
    faHeadersList.add(gtIdentifier);
    faHeadersList.addAll(GTIdentifierToSamples.get(gtIdentifier));
    faSequenceBuilder.setLength(0);
    for (String p : perPositionContents.keySet()) {
    if (!(cliarguments.skipConserved && perPositionContents.get(p).keySet().size() == 1)) {
    if (perPositionContents.get(p).containsKey(gtIdentifier)) {
    faSequenceBuilder.append(perPositionContents.get(p).get(gtIdentifier));
    } else {
    if (perPositionContents.get(p).containsKey("GT0")) {
    faSequenceBuilder.append(perPositionContents.get(p).get("GT0"));
    } else {
    faSequenceBuilder.append(Bio.GAP);
    }
    }
    }
    }
    faSequences.put(faSequenceBuilder.toString(), faHeadersList);
    }
    IO.writeFasta(
    new File(cliarguments.outputDirectory.getAbsolutePath() + "/GenotypeSequences/" + featureEntry.chromosome + "_" + featureEntry.name + "_" + cliarguments.samples.size() + "_genotypeSequencesAlignment.fasta"),
    faSequences
    );
    pb.step();
    }
    }
    }
    // Write proteoform sequences MSA.
    TreeSet<String> faEntryIds;
    if (cliarguments.PSMSA) {
    Logging.logStatus("Write proteoform sequences alignment.");
    try (ProgressBar pb = buildProgress()) {
    pb.maxHint(cliarguments.inputVDict.features.size());
    for (FeatureEntry featureEntry : cliarguments.inputVDict.features.values()) {
    perPositionContents.clear();
    HashMap<String, String> variants;
    faEntryIds = new TreeSet<>();
    for (ProteoformEntry proteoform : featureEntry.allocatedProtein.proteoforms.values()) {
    if (proteoform.name.equals(AllocatedProteinEntry.WILD_TYPE_PROTEOFORM_ID)) {
    sequenceChars = cliarguments.inputVDict.getProteoformSequence(featureEntry.name, proteoform.name).toCharArray();
    for (int i = 0; i < sequenceChars.length; i++) {
    if (!perPositionContents.containsKey((i + 1) + "+0")) {
    perPositionContents.put((i + 1) + "+0", new HashMap<>());
    }
    perPositionContents.get((i + 1) + "+0").put(proteoform.name, sequenceChars[i]);
    }
    faEntryIds.add(proteoform.name);
    } else {
    variants = cliarguments.inputVDict.getProteoformAminoacidVariants(featureEntry.name, proteoform.name);
    for (Map.Entry<String, String> variantEntry : variants.entrySet()) {
    positionStr = variantEntry.getKey();
    positionContent = variantEntry.getValue();
    if (!perPositionContents.containsKey(positionStr)) {
    perPositionContents.put(positionStr, new HashMap<>());
    }
    perPositionContents.get(positionStr).put(proteoform.name, positionContent.charAt(0));
    }
    faEntryIds.add(proteoform.name);
    }
    }
    faSequences = new HashMap<>();
    for (String faEntryId : faEntryIds) {
    faSequenceBuilder.setLength(0);
    faHeadersList = new ArrayList<>();
    for (String p : perPositionContents.keySet()) {
    if (!(cliarguments.skipConserved && perPositionContents.get(p).keySet().size() == 1)) {
    if (perPositionContents.get(p).containsKey(faEntryId)) {
    faSequenceBuilder.append(perPositionContents.get(p).get(faEntryId));
    } else {
    if (perPositionContents.get(p).containsKey(VariantsDictionary.WILD_TYPE_SAMPLE_ID)) {
    faSequenceBuilder.append(perPositionContents.get(p).get(VariantsDictionary.WILD_TYPE_SAMPLE_ID));
    } else {
    faSequenceBuilder.append(Bio.GAP);
    }
    }
    }
    }
    faHeadersList.add(faEntryId);
    faHeadersList.add(
    "VP=" + featureEntry.allocatedProtein.proteoforms.get(faEntryId).annotations.get("VP")
    + ";PT=" + featureEntry.allocatedProtein.proteoforms.get(faEntryId).annotations.get("PT")
    );
    faHeadersList.addAll(
    featureEntry.allocatedProtein.proteoforms.get(faEntryId).samples
    );
    faSequences.put(faSequenceBuilder.toString(), faHeadersList);
    }
    IO.writeFasta(
    new File(cliarguments.outputDirectory.getAbsolutePath() + "/ProteoformSequences/" + featureEntry.chromosome + "_" + featureEntry.name + "_" + cliarguments.samples.size() + "_proteoformSequencesAlignment.fasta"),
    faSequences
    );
    pb.step();
    }
    }
    }
    }
     */

}