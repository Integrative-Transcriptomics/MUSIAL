package main;

import cli.CLIParameters;
import cli.CLIParametersInferSequences;
import cli.CLIParametersStatistics;
import cli.CLIParametersUpdateVDict;
import components.*;
import datastructure.*;
import exceptions.MusialBioException;
import exceptions.MusialCLException;
import exceptions.MusialIOException;
import exceptions.MusialIntegrityException;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
import runnables.SampleAnalyzerRunnable;

import java.io.*;
import java.util.*;
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
     * Factory to generate {@link ProgressBar} instances with predefined properties.
     */
    private static final ProgressBarBuilder progressBarBuilder =
            new ProgressBarBuilder().setStyle(ProgressBarStyle.ASCII).setMaxRenderedLength(90);
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
    public static boolean COMPRESS = true;

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
                case inferSequences -> runInferSequences((CLIParametersInferSequences) arguments);
                case statistics -> runStatistics((CLIParametersStatistics) arguments);
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
        System.out.println("\tinferSequences : Infer sequence information from an existing variants dictionary.");
        System.out.println("\tstatistics : Run SnpEff annotation and compute various statistics.");
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
     * Generates and returns {@link ProgressBar} instance.
     *
     * @return Instance of {@link ProgressBar} with the options specified by the progressBarBuilder property.
     */
    public static ProgressBar buildProgress() {
        return progressBarBuilder.build();
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
     * @throws MusialCLException        See documentation of {@link CLIParametersUpdateVDict} class.
     * @throws MusialIOException        See documentation of {@link CLIParametersUpdateVDict} class.
     * @throws MusialBioException       See documentation of {@link CLIParametersUpdateVDict} class.
     * @throws MusialIntegrityException See documentation of {@link CLIParametersUpdateVDict} class.
     */
    private static CLIParameters parseCLIArguments(String[] args)
            throws MusialCLException, MusialIOException, MusialBioException, MusialIntegrityException {
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
     * @throws InterruptedException     Thrown if any {@link Runnable} implementing class instance is interrupted by the user.
     * @throws MusialIOException        Thrown if any input or output file is missing or unable to being generated (caused by any MUSIAL method).
     * @throws IOException              Thrown if any input or output file is missing or unable to being generated (caused by any native Java method).
     * @throws MusialIntegrityException Thrown if any method fails wrt. internal logic, i.e. assignment of proteins to genomes.
     * @throws MusialBioException       Thrown if any method fails wrt. biological context, i.e. parsing of unknown symbols.
     */
    private static void runUpdateVDict(CLIParametersUpdateVDict cliarguments)
            throws InterruptedException, MusialIOException, MusialIntegrityException, MusialBioException, IOException {
        // Read-in existing variants dictionary or build new one.
        VariantsDictionary variantsDictionary = VariantsDictionaryFactory.build(cliarguments);
        // Collect feature entry information.
        Set<String> featureIdsToUpdate;
        /* TODO: Collect String Set of features to be removed; delete all information unique for removed features.
         */
        Logging.logStatus("Updating feature information.");
        try (ProgressBar pb = buildProgress()) {
            featureIdsToUpdate = cliarguments.features.keySet();
            pb.maxHint(featureIdsToUpdate.size());
            FastaContainer referenceChromosomeFastaContainer = null;
            for (FastaContainer fC : IO.readFastaToSet(cliarguments.referenceFASTA)) {
                String chromosome = fC.getHeader().split(" ")[0].trim();
                if (variantsDictionary.chromosome.equals(chromosome)) {
                    referenceChromosomeFastaContainer = fC;
                    break;
                }
            }
            if (referenceChromosomeFastaContainer == null) {
                throw new MusialIOException("Failed to match feature chromosome " + variantsDictionary.chromosome + " to specified reference .fasta");
            }
            for (String fId : featureIdsToUpdate) {
                if (!variantsDictionary.features.containsKey(fId)) {
                    FeatureEntry.imputeSequence(cliarguments.features.get(fId), referenceChromosomeFastaContainer);
                    if (cliarguments.features.get(fId).pdbFile != null) {
                        FeatureEntry.imputeProtein(cliarguments.features.get(fId));
                    }
                    variantsDictionary.features.put(fId, cliarguments.features.get(fId));
                }
                pb.step();
            }
            pb.setExtraMessage(Logging.getDoneMessage());
        }

        // Collect sample variants information.
        /* TODO: Collect String Set of samples to be removed; delete all information unique for removed samples.
         */
        Logging.logStatus("Updating variants from samples.");
        try (ProgressBar pb = buildProgress()) {
            pb.maxHint(cliarguments.samples.size());
            ExecutorService executor = Executors.newFixedThreadPool(Musial.THREADS);
            for (SampleEntry sampleEntry : cliarguments.samples.values()) {
                SampleEntry.imputeVCFFileReader(sampleEntry);
                variantsDictionary.samples.put(sampleEntry.name, sampleEntry);
                for (String fId : featureIdsToUpdate) {
                    executor.execute(
                            new SampleAnalyzerRunnable(sampleEntry, cliarguments.features.get(fId), variantsDictionary, pb)
                    );
                }
                pb.step();
            }
            executor.shutdown();
            //noinspection ResultOfMethodCallIgnored
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            pb.setExtraMessage(Logging.getDoneMessage());
        }

        // Infer proteoform information, if features have assigned .pdb files.
        ConcurrentSkipListMap<String, String> variants;
        long noFeaturesWithAssignedProteins = variantsDictionary.features.keySet().stream().filter(
                feKey -> variantsDictionary.features.get(feKey).pdbFile != null
        ).count();
        if (noFeaturesWithAssignedProteins > 0) {
            try (ProgressBar pb = buildProgress()) {
                Logging.logStatus("Infer proteoform information.");
                pb.maxHint(noFeaturesWithAssignedProteins * variantsDictionary.samples.size());
                for (String fId : variantsDictionary.features.keySet()) {
                    if (variantsDictionary.features.get(fId).pdbFile == null) {
                        continue;
                    }
                    for (String sId : variantsDictionary.samples.keySet()) {
                        pb.setExtraMessage("Processing sample " + sId + " / feature " + fId);
                        variants = Bio.inferProteoform(variantsDictionary, fId, sId);
                        variantsDictionary.features.get(fId).allocatedProtein
                                .addProteoform(variantsDictionary.features.get(fId), sId, variants);
                        pb.step();
                    }
                }
            }
        }

        // Write updated database to file.
        try (ProgressBar pb = buildProgress()) {
            Logging.logStatus("Write variants dictionary to file.");
            pb.maxHint(1);
            String outputFile = cliarguments.outputFile.getAbsolutePath();
            File vDictOutfile = new File(outputFile);
            if (!vDictOutfile.exists()) {
                IO.generateFile(vDictOutfile);
            }
            variantsDictionary.dump(vDictOutfile);
            pb.step();
            pb.setExtraMessage(Logging.getDoneMessage());
        }
    }

    /**
     * Extract genotype and/or proteoform sequences as alignment or unaligned in .fasta format from an existing variants dictionary.
     * <p>
     * Sequences are written in coding direction.
     *
     * @param cliarguments {@link CLIParametersInferSequences} instance yielding parameter specification for the MUSIAL infer sequences module.
     * @throws MusialIOException  If the generation of output directories or files fails.
     * @throws MusialBioException If sequence extraction procedures from the specified {@link VariantsDictionary} instance fail.
     */
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
                            sequence = cliarguments.inputVDict.getNucleotideSequence(featureEntry.name, sampleEntry.name);
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
                                variants = cliarguments.inputVDict.getSampleVariants(featureEntry.name, sampleEntry.name);
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
                            variants = cliarguments.inputVDict.getProteoformVariants(featureEntry.name, proteoform.name);
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

    /**
     * @param cliarguments
     */
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

}