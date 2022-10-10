package main;

import cli.CLIParameters;
import cli.CLIParametersInferSequences;
import cli.CLIParametersStatistics;
import cli.CLIParametersUpdateVDict;
import components.*;
import datastructure.*;
import exceptions.MusialException;
import me.tongfei.progressbar.ProgressBar;
import runnables.SampleAnalyzerRunnable;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

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
     * {@link String} representing the logo of the software.
     */
    public static String LOGO = "";
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
     * Boolean value to indicate if debugging mode is enabled.
     */
    public static boolean DEBUG = false;
    /**
     * Number of threads to use.
     */
    public static int THREADS = 1;
    /**
     * Whether to compress output files.
     */
    public static boolean COMPRESS = false;

    /**
     * Main method of MUSIAL; loads meta-information and invokes methods dependent on the user specified module.
     *
     * @param args {@link String[]} comprising the arguments parsed from the command line.
     */
    public static void main(String[] args) {
        try {
            if (args.length == 0) {
                loadMetadata();
                System.out.println("| no module specified.");
                printModuleInformation();
            }
            // 1. Retain user specified module and parse command line arguments.
            for (MusialModules mm : MusialModules.values()) {
                if (mm.name().equals(args[0])) {
                    MODULE = mm;
                    break;
                }
            }
            if (MODULE == null) {
                loadMetadata();
                System.out.println("| unknown module " + args[0] + ".");
                printModuleInformation();
            } else {
                loadMetadata();
            }
            CLIParameters arguments = parseCLIArguments(args);
            // 2. Invokes methods dependent on the user specified module.
            switch (MODULE) {
                case updateVDict -> runUpdateVDict((CLIParametersUpdateVDict) arguments);
                //case inferSequences -> runInferSequences((CLIParametersInferSequences) arguments);
                //case statistics -> runStatistics((CLIParametersStatistics) arguments);
            }
        } catch (Exception e) {
            e.printStackTrace();
            Logging.logError(e.getMessage());
        }
    }

    /**
     * Prints information about available modules to the user and exits.
     */
    public static void printModuleInformation() {
        System.out.println("| available MUSIAL modules:");
        System.out.println("\tupdateVDict : Generate a new or update an existing variants dictionary JSON file.");
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
        Musial.LOGO = properties.getProperty("logo");
        Musial.VERSION = properties.getProperty("version");
        Musial.CONTACT = properties.getProperty("contact");
        Musial.LICENSE = properties.getProperty("license");
        // Print information to stdout.
        System.out.println(Musial.LOGO + Musial.VERSION);
        System.out.println(Musial.LICENSE + ", Contact: " + Musial.CONTACT);
    }

    /**
     * Parses command line arguments by instantiating an object implementing the {@link CLIParameters} interface.
     * <p>
     * The command line arguments specified by the user and accepted by the `main` method are passed and used to
     * instantiate an object implementing the {@link CLIParameters} interface which is returned if all argument
     * validation steps were successful.
     * If any error arises during the parsing step the program exits.
     *
     * @param args Arguments parsed from the command line.
     * @return An instance of the {@link CLIParametersUpdateVDict} class.
     * @throws MusialException See documentation of {@link CLIParameters} implementing classes.
     */
    private static CLIParameters parseCLIArguments(String[] args)
            throws MusialException {
        CLIParameters cliParameters = null;
        switch (MODULE) {
            case updateVDict -> cliParameters = new CLIParametersUpdateVDict(args);
            case inferSequences -> cliParameters = new CLIParametersInferSequences(args);
            case statistics -> cliParameters = new CLIParametersStatistics(args);
        }
        assert cliParameters != null;
        return cliParameters;
    }

    /**
     * Generates a new or updates an existing variants dictionary JSON file based on the specifications parsed from a variants dictionary specification JSON file.
     *
     * @param cliarguments {@link CLIParametersUpdateVDict} instance yielding parameter specification for the MUSIAL update variants dictionary module.
     * @throws InterruptedException Thrown if any {@link Runnable} implementing class instance is interrupted by the user.
     * @throws IOException          Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialException      If any method fails wrt. biological context, i.e. parsing of unknown symbols; If any method fails wrt. internal logic, i.e. assignment of proteins to genomes; If any input or output file is missing or unable to being generated.
     */
    private static void runUpdateVDict(CLIParametersUpdateVDict cliarguments)
            throws InterruptedException, MusialException, IOException {
        try (
                ProgressBar progressBar_DumpData = new ProgressBar("Dump Updated Data", 0);
                ProgressBar progressBar_InferCodingFeatureInformation = new ProgressBar("Infer Coding Feature Information", 0);
                ProgressBar progressBar_InferFeatureAlleles = new ProgressBar("Infer Feature Alleles", 0);
                ProgressBar progressBar_UpdateSampleInformation = new ProgressBar("Update Sample/Variant Information", 0);
                ProgressBar progressBar_UpdateFeatureInformation = new ProgressBar("Update Feature Information", 0);
        ) {
            // Read-in existing variants dictionary or build new one.
            VariantsDictionary variantsDictionary = VariantsDictionaryFactory.build(cliarguments);

            // Collect feature information.
            Set<String> featureIdsUpdate = cliarguments.features.keySet();
            HashSet<String> featureIdsAll = new HashSet<>();
            featureIdsAll.addAll(featureIdsUpdate);
            featureIdsAll.addAll(variantsDictionary.features.keySet());
            progressBar_UpdateFeatureInformation.maxHint(featureIdsAll.size());
            FastaContainer referenceChromosomeFastaContainer = null;
            for (FastaContainer fastaContainer : IO.readFastaToSet(cliarguments.referenceFASTA)) {
                String chromosome = fastaContainer.getHeader().split(" ")[0].trim();
                if (variantsDictionary.chromosome.equals(chromosome)) {
                    referenceChromosomeFastaContainer = fastaContainer;
                    break;
                }
            }
            if (referenceChromosomeFastaContainer == null) {
                throw new MusialException("(Bio) Failed to match feature chromosome " + variantsDictionary.chromosome + " to specified reference .fasta");
            }
            for (String featureId : featureIdsAll) {
                if (!featureIdsUpdate.contains(featureId) && variantsDictionary.features.containsKey(featureId)) {
                    // TODO: Remove feature information.
                    progressBar_UpdateFeatureInformation.step();
                } else {
                    if (!variantsDictionary.features.containsKey(featureId)) {
                        cliarguments.features.get(featureId).imputeNucleotideSequence(referenceChromosomeFastaContainer);
                        if (cliarguments.features.get(featureId).pdbFile != null) {
                            cliarguments.features.get(featureId).imputeProteinInformation();
                        }
                        variantsDictionary.features.put(featureId, cliarguments.features.get(featureId));
                    }
                }
                progressBar_UpdateFeatureInformation.step();
            }
            progressBar_UpdateFeatureInformation.setExtraMessage(Logging.getDoneMessage());

            // Collect information about sample variants.
            Set<String> sampleIdsUpdate = cliarguments.samples.keySet();
            HashSet<String> sampleIdsAll = new HashSet<>();
            sampleIdsAll.addAll(sampleIdsUpdate);
            sampleIdsAll.addAll(variantsDictionary.samples.keySet());
            progressBar_UpdateSampleInformation.maxHint(sampleIdsAll.size());
            ExecutorService executor = Executors.newFixedThreadPool(Musial.THREADS);
            for (String sampleId : sampleIdsAll) {
                if (!sampleIdsUpdate.contains(sampleId) && variantsDictionary.samples.containsKey(sampleId)) {
                    // TODO: Remove sample information.
                    progressBar_UpdateSampleInformation.step();
                } else {
                    SampleEntry sampleEntry = cliarguments.samples.get(sampleId);
                    SampleEntry.imputeVCFFileReader(sampleEntry);
                    variantsDictionary.samples.put(sampleEntry.name, sampleEntry);
                    for (String featureId : featureIdsUpdate) {
                        executor.execute(
                                new SampleAnalyzerRunnable(sampleEntry, cliarguments.features.get(featureId), variantsDictionary, progressBar_UpdateSampleInformation)
                        );
                    }
                }
            }
            executor.shutdown();
            //noinspection ResultOfMethodCallIgnored
            executor.awaitTermination(30, TimeUnit.MINUTES);
            progressBar_UpdateSampleInformation.setExtraMessage(Logging.getDoneMessage());

            // Infer feature alleles from collected variant information.
            progressBar_InferFeatureAlleles.maxHint(variantsDictionary.features.keySet().size());
            for (FeatureEntry featureEntry : variantsDictionary.features.values()) {
                featureEntry.inferAlleleInformation(variantsDictionary);
                progressBar_InferFeatureAlleles.step();
            }
            progressBar_InferFeatureAlleles.setExtraMessage(Logging.getDoneMessage());

            // Infer coding feature variant/proteoform information; For all coding features.
            ConcurrentSkipListMap<String, String> variants;
            long noFeaturesWithAssignedProteins = variantsDictionary.features.keySet().stream().filter(
                    feKey -> variantsDictionary.features.get(feKey).isCodingSequence
            ).count();
            progressBar_InferCodingFeatureInformation.maxHint(noFeaturesWithAssignedProteins * variantsDictionary.samples.size());
            if (noFeaturesWithAssignedProteins > 0) {
                for (String featureId : variantsDictionary.features.keySet()) {
                    for (String sampleId : variantsDictionary.samples.keySet()) {
                        variants = Bio.inferProteoform(variantsDictionary, featureId, sampleId);
                        variantsDictionary.features.get(featureId).addProteoform(sampleId, variants, variantsDictionary);
                        progressBar_InferCodingFeatureInformation.step();
                    }
                }
            }
            progressBar_InferCodingFeatureInformation.setExtraMessage(Logging.getDoneMessage());

            // Write updated database to file.
            progressBar_DumpData.maxHint(1);
            String outputFile = cliarguments.outputFile.getAbsolutePath();
            File vDictOutfile = new File(outputFile);
            if (!vDictOutfile.exists()) {
                IO.generateFile(vDictOutfile);
            }
            variantsDictionary.dump(vDictOutfile);
            progressBar_DumpData.step();
            progressBar_DumpData.setExtraMessage(Logging.getDoneMessage());
        }
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