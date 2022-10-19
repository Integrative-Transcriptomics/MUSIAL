package main;

import cli.CLIParser;
import cli.ModuleParametersBuild;
import components.*;
import datastructure.*;
import exceptions.MusialException;
import org.json.simple.JSONObject;
import runnables.SampleAnalyzerRunnable;

import java.io.*;
import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

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
     * Original output stream.
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
     * Number of threads to use.
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

    public static final DecimalFormat DECIMAL_FORMATTER = new DecimalFormat("#.#");

    /**
     * Main method of MUSIAL; loads meta-information and invokes methods dependent on the user specified module.
     *
     * @param args {@link String[]} comprising the arguments parsed from the command line.
     */
    public static void main(String[] args) {
        try {
            loadMetadata();
            if (args.length == 0) {
                Logging.logError("No arguments were specified");
                printInfo();
            } else {
                CLIParser cliParser = new CLIParser(args);
                if (cliParser.configuration.keySet().size() == 0) {
                    Logging.logError("No modules were specified in " + Logging.colorParameter(cliParser.arguments.getOptionValue("c")));
                    printInfo();
                }
                for (Object o : cliParser.configuration.keySet()) {
                    MODULE = null;
                    String moduleId = String.valueOf(o);
                    // Match specified modules with implemented modules.
                    for (MusialModules mm : MusialModules.values()) {
                        if (mm.name().equals(moduleId)) {
                            MODULE = mm;
                            break;
                        }
                    }
                    if (MODULE == null) {
                        Logging.logWarning("Skip unknown module " + Logging.colorParameter(moduleId));
                    } else {
                        JSONObject parameters = (JSONObject) cliParser.configuration.get(moduleId);
                        switch (MODULE) {
                            case BUILD -> executeBUILD(new ModuleParametersBuild(parameters));
                            //case inferSequences -> runInferSequences((CLIParametersInferSequences) arguments);
                            //case statistics -> runStatistics((CLIParametersStatistics) arguments);
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
        System.out.println("| available MUSIAL modules:");
        System.out.println("\tBUILD : Generates a variants dictionary JSON file.");
        //System.out.println("\tinferSequences : Infer sequence information from an existing variants dictionary.");
        //System.out.println("\tstatistics : Run SnpEff annotation and compute various statistics.");
        System.out.println("| specify -h for any module to obtain more information.");
        System.exit(0);
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
     * @param cliarguments {@link ModuleParametersBuild} instance yielding parameter specification for the MUSIAL update variants dictionary module.
     * @throws InterruptedException Thrown if any {@link Runnable} implementing class instance is interrupted by the user.
     * @throws IOException          Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialException      If any method fails wrt. biological context, i.e. parsing of unknown symbols; If any method fails wrt. internal logic, i.e. assignment of proteins to genomes; If any input or output file is missing or unable to being generated.
     */
    private static void executeBUILD(ModuleParametersBuild cliarguments)
            throws InterruptedException, MusialException, IOException {
        Logging.logStatus("Execute module " + Logging.getCustomTag("UPDATE"));

        // Read-in existing variants dictionary or build new one.
        VariantsDictionary variantsDictionary = VariantsDictionaryFactory.build(cliarguments);

        // Update feature information.
        Logging.logStatus(Logging.getStartTag() + " Update feature information");
        Set<String> specifiedFeatures = cliarguments.features.keySet();
        HashSet<String> featureIdsAll = new HashSet<>();
        featureIdsAll.addAll(specifiedFeatures);
        featureIdsAll.addAll(variantsDictionary.features.keySet());
        //Logging.logStatus("Match chromosome " + Logging.colorParameter(variantsDictionary.chromosome) + " from " + Logging.colorParameter(cliarguments.referenceFASTA.getAbsolutePath()));

        FastaContainer referenceChromosomeFastaContainer = null;
        for (FastaContainer fastaContainer : IO.readFastaToSet(cliarguments.referenceFASTA)) {
            String chromosome = fastaContainer.getHeader().split(" ")[0].trim();
            if (variantsDictionary.chromosome.equals(chromosome)) {
                referenceChromosomeFastaContainer = fastaContainer;
                break;
            }
        }
        if (referenceChromosomeFastaContainer == null) {
            throw new MusialException("Failed to match chromosome " + Logging.colorParameter(variantsDictionary.chromosome) + " from " + Logging.colorParameter(cliarguments.referenceFASTA.getAbsolutePath()));
        }

        for (String featureId : featureIdsAll) {
            if (!specifiedFeatures.contains(featureId) && variantsDictionary.features.containsKey(featureId)) {
                // TODO: Remove feature information.
            } else {
                if (!variantsDictionary.features.containsKey(featureId)) {
                    cliarguments.features.get(featureId).imputeNucleotideSequence(referenceChromosomeFastaContainer);
                    if (cliarguments.features.get(featureId).pdbFile != null) {
                        cliarguments.features.get(featureId).imputeProteinInformation();
                    }
                    variantsDictionary.features.put(featureId, cliarguments.features.get(featureId));
                }
            }
        }
        Logging.logStatus(Logging.getDoneTag() + " Update feature information");

        // Collect information about sample variants.
        Logging.logStatus("Collect sample information " + Logging.getStartTag());
        Set<String> sampleIdsUpdate = cliarguments.samples.keySet();
        HashSet<String> sampleIdsAll = new HashSet<>();
        sampleIdsAll.addAll(sampleIdsUpdate);
        sampleIdsAll.addAll(variantsDictionary.samples.keySet());
        ExecutorService executor = Executors.newFixedThreadPool(Musial.THREADS);
        for (String sampleId : sampleIdsAll) {
            if (!sampleIdsUpdate.contains(sampleId) && variantsDictionary.samples.containsKey(sampleId)) {
                // TODO: Remove sample information.
            } else {
                SampleEntry sampleEntry = cliarguments.samples.get(sampleId);
                SampleEntry.imputeVCFFileReader(sampleEntry);
                variantsDictionary.samples.put(sampleEntry.name, sampleEntry);
                for (String featureId : specifiedFeatures) {
                    executor.execute(
                            new SampleAnalyzerRunnable(sampleEntry, cliarguments.features.get(featureId), variantsDictionary)
                    );
                }
            }
        }
        executor.shutdown();
        //noinspection ResultOfMethodCallIgnored
        executor.awaitTermination(30, TimeUnit.MINUTES);
        Logging.logStatus("Collect feature information " + Logging.getDoneTag());

        // Run SnpEff annotation for variants.
        Logging.logStatus("Annotate novel variants with SnpEff " + Logging.getStartTag());
        File tmpDirectory = new File("./tmp/");
        try {
            IO.generateDirectory(tmpDirectory);
            File variantsFile = new File("./tmp/variants.vcf");
            IO.generateFile(variantsFile);
            IO.writeVcf(variantsFile, variantsDictionary);
            SnpEffAnnotator.runSnpEff(
                    tmpDirectory,
                    variantsFile,
                    cliarguments.referenceFASTA,
                    cliarguments.referenceGFF,
                    variantsDictionary.chromosome
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
        Logging.logStatus("Annotate novel variants with SnpEff " + Logging.getDoneTag());

        // Infer feature alleles from collected variant information.
        Logging.logStatus("Infer alleles " + Logging.getStartTag());
        for (FeatureEntry featureEntry : variantsDictionary.features.values()) {
            featureEntry.inferAlleleInformation(variantsDictionary);
        }
        Logging.logStatus("Infer alleles " + Logging.getDoneTag());

        // Infer coding feature variant/proteoform information; For all coding features.
        Logging.logStatus("Infer protein variants " + Logging.getStartTag());
        ConcurrentSkipListMap<String, String> variants;
        long noFeaturesWithAssignedProteins = variantsDictionary.features.keySet().stream().filter(
                feKey -> variantsDictionary.features.get(feKey).isCodingSequence
        ).count();
        if (noFeaturesWithAssignedProteins > 0) {
            for (String featureId : variantsDictionary.features.keySet()) {
                for (String sampleId : variantsDictionary.samples.keySet()) {
                    variants = Bio.inferProteoform(variantsDictionary, featureId, sampleId);
                    variantsDictionary.features.get(featureId).addProteoform(sampleId, variants, variantsDictionary);
                }
            }
        }
        Logging.logStatus("Infer protein variants " + Logging.getDoneTag());

        // Infer variant frequencies.
        Logging.logStatus("Infer variant statistics " + Logging.getStartTag());
        int totalNoVariants = 0;
        for (ConcurrentSkipListMap<String, NucleotideVariantEntry> perPositionNucleotideVariants : variantsDictionary.nucleotideVariants.values()) {
            totalNoVariants += perPositionNucleotideVariants.size();
        }
        for (FeatureEntry feature : variantsDictionary.features.values()) {
            totalNoVariants += feature.aminoacidVariants.size();
        }
        float noSamples = variantsDictionary.samples.size();
        float noOccurrences;
        // Compute frequencies of nucleotide variants.
        for (ConcurrentSkipListMap<String, NucleotideVariantEntry> perPositionNucleotideVariants : variantsDictionary.nucleotideVariants.values()) {
            for (NucleotideVariantEntry nucleotideVariant : perPositionNucleotideVariants.values()) {
                noOccurrences = 0;
                for (String occurrence : nucleotideVariant.occurrence) {
                    String featureId = occurrence.split(VariantsDictionary.FIELD_SEPARATOR_1)[0];
                    String alleleId = occurrence.split(VariantsDictionary.FIELD_SEPARATOR_1)[1];
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
        Logging.logStatus("Infer variant statistics " + Logging.getDoneTag());

        // Write updated database to file.
        Logging.logStatus("Dump updated file " + Logging.getStartTag());
        String outputFile = cliarguments.outputFile.getAbsolutePath();
        File vDictOutfile = new File(outputFile);
        if (!vDictOutfile.exists()) {
            IO.generateFile(vDictOutfile);
        }
        variantsDictionary.dump(vDictOutfile);
        Logging.logStatus("Dump updated file " + Logging.getDoneTag());
    }

    /*
     **
     * Extract genotype and/or proteoform sequences as alignment or unaligned in .fasta format from an existing variants dictionary.
     * <p>
     * Sequences are written in coding direction.
     *
     * @param cliarguments {@link CLIParametersInferSequences} instance yielding parameter specification for the MUSIAL infer sequences module.
     * @throws MusialIOException  If the generation of output directories or files fails.
     * @throws MusialBioException If sequence extraction procedures from the specified {@link VariantsDictionary} instance fail.
     *
    private static void runInferSequences(CLIParametersInferSequences cliarguments) throws MusialIOException, MusialBioException {
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

     **
     * @param cliarguments
     *
    private static void runStatistics(CLIParametersStatistics cliarguments) throws MusialIOException, IOException {
        // Run SnpEff on all variants wrt. the samples and features specified as input.
        File tempDir = new File("./temp/");
        try {
            Logging.logStatus("Annotating variants with SnpEff.");
            IO.generateDirectory(tempDir);
            File variantsFile = new File("./temp/variants.vcf");
            IO.generateFile(variantsFile);
            HashSet<String> excludedFeatures = new HashSet<>();
            if (cliarguments.features.size() == 0) {
                excludedFeatures.addAll(cliarguments.inputVDict.features.keySet());
                excludedFeatures.removeAll(cliarguments.features);
            }
            HashSet<String> excludedSamples = new HashSet<>();
            if (cliarguments.samples.size() == 0) {
                excludedSamples.addAll(cliarguments.inputVDict.samples.keySet());
                excludedSamples.removeAll(cliarguments.samples);
            }
            IO.writeVcf(variantsFile, cliarguments.inputVDict, excludedFeatures, excludedSamples);
            SnpEffAnnotator.runSnpEff(
                    tempDir,
                    variantsFile,
                    cliarguments.referenceFASTA,
                    cliarguments.referenceGFF,
                    cliarguments.inputVDict.chromosome
            );
            ArrayList<String> outputFileLineContent = new ArrayList<>();
            outputFileLineContent.add("POSITION\tREF_CONTENT\tALT_CONTENT\tSAMPLES\tINFO\n");
            List<String> annotations = IO.readLinesFromFile("./temp/annotated_variants.vcf").stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());
            String[] splitAnnotation;
            String position;
            String referenceContent;
            String alternativeContent;
            String[] annotationFields;
            ConcurrentSkipListMap<String, String> occurrence;
            ArrayList<String> samples = new ArrayList<>();
            boolean isRejected;
            ArrayList<Float> qualities = new ArrayList<>();
            ArrayList<Float> frequencies = new ArrayList<>();
            ArrayList<Float> coverages = new ArrayList<>();
            for (String annotation : annotations) {
                splitAnnotation = annotation.split("\t");
                position = splitAnnotation[1];
                referenceContent = splitAnnotation[3];
                alternativeContent = splitAnnotation[4];
                if (referenceContent.length() > alternativeContent.length()) {
                    alternativeContent += String.valueOf(Bio.DELETION_AA1).repeat(referenceContent.length() - alternativeContent.length());
                }
                if (splitAnnotation[7].equals(".")) {
                    continue;
                }
                annotationFields = splitAnnotation[7].split("=")[1].split("\\|");
                occurrence = cliarguments.inputVDict.variants.get(Integer.parseInt(position)).get(alternativeContent).occurrence;
                samples.clear();
                qualities.clear();
                frequencies.clear();
                coverages.clear();
                for (Map.Entry<String, String> occurenceEntry : occurrence.entrySet()) {
                    isRejected = occurenceEntry.getValue().split("\\|")[0].equals("true");
                    if (!isRejected) {
                        samples.add(occurenceEntry.getKey());
                        qualities.add(Float.parseFloat(occurenceEntry.getValue().split("\\|")[2]));
                        frequencies.add(Float.parseFloat(occurenceEntry.getValue().split("\\|")[3]));
                        coverages.add(Float.parseFloat(occurenceEntry.getValue().split("\\|")[4]));
                    }
                }
                if (samples.size() > 0) {
                    outputFileLineContent.add(
                            position + "\t"
                                    + referenceContent + "\t"
                                    + alternativeContent + "\t"
                                    + String.join(",", samples) + "\t"
                                    + "TYPE=" + annotationFields[1] + ","
                                    + "IMPACT=" + annotationFields[2] + ","
                                    + "MEAN_QUALITY=" + qualities.stream().mapToDouble(f -> f).sum() / qualities.size() + ","
                                    + "MEAN_FREQUENCY=" + frequencies.stream().mapToDouble(f -> f).sum() / frequencies.size() + ","
                                    + "MEAN_COVERAGE=" + coverages.stream().mapToDouble(f -> f).sum() / coverages.size() + "\n"
                    );
                }
            }
            File outputFile = new File(cliarguments.outputDirectory.getAbsolutePath() + "/" + "nucleotideVariants.tsv");
            IO.generateFile(outputFile);
            IO.writeFile(outputFile, outputFileLineContent);
        } finally {
            // IO.deleteDirectory(tempDir);
        }
    }
    */
}