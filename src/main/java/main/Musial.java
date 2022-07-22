package main;

import cli.CLIParameters;
import cli.CLIParametersInferSequences;
import cli.CLIParametersUpdateVDict;
import components.*;
import datastructure.FastaContainer;
import datastructure.FeatureEntry;
import datastructure.SampleEntry;
import datastructure.VariantsDictionary;
import exceptions.MusialBioException;
import exceptions.MusialCLException;
import exceptions.MusialIOException;
import exceptions.MusialIntegrityException;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
import runnables.SampleAnalyzerRunnable;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;
import java.util.Set;
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

        /* FIXME: SnpEff annotation is currently disabled.
        Logging.logStatus("Annotating novel variants with SnpEff.");
        File tempDir = new File("./temp/");
        try {
            IO.generateDirectory(tempDir);
            File novelVariantsOut = new File("./temp/novelVariants.vcf");
            IO.generateFile(novelVariantsOut);
            IO.writeVcf(novelVariantsOut, variantsDictionary.novelVariants, variantsDictionary.chromosome);
            SnpEffAnnotator.runSnpEff(tempDir, novelVariantsOut, cliarguments.referenceFASTA, cliarguments.referenceGFF,
                    variantsDictionary.chromosome);
            for (String line : IO.readLinesFromFile("./temp/SnpEff.vcf")) {
                if (line.startsWith("#")) {
                    continue;
                }
                String[] lineFields = line.split("\t");
                String position = lineFields[1];
                String refContent = lineFields[3];
                String altContent = lineFields[4];
                if (refContent.length() > altContent.length()) {
                    altContent = altContent + "-".repeat(refContent.length() - altContent.length());
                }
                variantsDictionary.addVariantAnnotation(
                        Integer.parseInt(position),
                        altContent,
                        SnpEffAnnotator.convertAnnotation(lineFields[7].split(";")[0].split("=")[1])
                );
            }
        } catch (Exception e) {
            throw e;
        } finally {
            IO.deleteDirectory(tempDir);
        }
         */

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

    private static void runInferSequences(CLIParametersInferSequences cliarguments) throws MusialIOException, MusialBioException {
        // TODO: Generalize methods; Currently code is copied for debugging.
        /*
        Write genotype sequences.
         */
        if (cliarguments.GS) {
            Logging.logStatus("Write genotype sequences.");
            try (ProgressBar pb = buildProgress()) {
                pb.maxHint(cliarguments.inputVDict.features.size());
                for (FeatureEntry featureEntry : cliarguments.inputVDict.features.values()) {
                    HashMap<String, ArrayList<String>> genotypeSequences = new HashMap<>();
                    // Access reference nucleotide sequence.
                    String featureWildTypeNucleotideSequence = featureEntry.nucleotideSequence;
                    genotypeSequences.put(featureWildTypeNucleotideSequence, new ArrayList<>());
                    genotypeSequences.get(featureWildTypeNucleotideSequence).add(VariantsDictionary.WILD_TYPE_SAMPLE_ID);
                    String genotypeSequence;
                    String sampleVSwab;
                    for (SampleEntry sampleEntry : cliarguments.inputVDict.samples.values()) {
                        sampleVSwab = sampleEntry.annotations.get(featureEntry.name + VariantsDictionary.ATTRIBUTE_VARIANT_SWAB_NAME);
                        if (sampleVSwab == null) {
                            genotypeSequences.get(featureWildTypeNucleotideSequence).add(sampleEntry.name);
                        } else {
                            genotypeSequence = cliarguments.inputVDict.getNucleotideSequence(featureEntry.name, sampleEntry.name);
                            if (genotypeSequences.containsKey(genotypeSequence)) {
                                genotypeSequences.get(genotypeSequence).add(sampleEntry.name);
                            } else {
                                genotypeSequences.put(genotypeSequence, new ArrayList<>());
                                genotypeSequences.get(genotypeSequence).add(sampleEntry.name);
                            }
                        }
                    }
                    int gtIndex = 1;
                    for (String s : genotypeSequences.keySet()) {
                        genotypeSequences.get(s).add(0, "GT" + gtIndex);
                        gtIndex += 1;
                    }
                    File gtsOutDir = new File(cliarguments.outputDirectory.getAbsolutePath() + "/GenotypeSequences/");
                    if (!Validation.isDirectory(gtsOutDir)) {
                        IO.generateDirectory(gtsOutDir);
                    }
                    IO.writeFasta(
                            new File(cliarguments.outputDirectory.getAbsolutePath() + "/GenotypeSequences/" + featureEntry.chromosome + "_" + featureEntry.name + "_" + cliarguments.samples.size() + "_genotypeSequences.fasta"),
                            genotypeSequences
                    );
                    pb.step();
                }
            }
        }
        /*
        Write proteoform sequences.
         */
        if (cliarguments.PS) {
            Logging.logStatus("Write proteoform sequences.");
            try (ProgressBar pb = buildProgress()) {
                pb.maxHint(cliarguments.inputVDict.features.size());
                for (FeatureEntry featureEntry : cliarguments.inputVDict.features.values()) {
                    HashMap<String, ArrayList<String>> proteoformSequences = new HashMap<>();
                    String proteoformSequence;
                    for (String pfId : cliarguments.inputVDict.features.get(featureEntry.name).allocatedProtein.proteoforms.keySet()) {
                        proteoformSequence = cliarguments.inputVDict.getProteoformSequence(featureEntry.name, pfId);
                        proteoformSequences.put(proteoformSequence, new ArrayList<>());
                        proteoformSequences.get(proteoformSequence).add(pfId);
                        proteoformSequences.get(proteoformSequence).add("VP=" + featureEntry.allocatedProtein.proteoforms.get(pfId).annotations.get("VP"));
                        proteoformSequences.get(proteoformSequence).addAll(featureEntry.allocatedProtein.proteoforms.get(pfId).samples);
                    }
                    File pfsOutDir = new File(cliarguments.outputDirectory.getAbsolutePath() + "/ProteoformSequences/");
                    if (!Validation.isDirectory(pfsOutDir)) {
                        IO.generateDirectory(pfsOutDir);
                    }
                    IO.writeFasta(
                            new File(cliarguments.outputDirectory.getAbsolutePath() + "/ProteoformSequences/" + featureEntry.chromosome + "_" + featureEntry.name + "_" + cliarguments.samples.size() + "_proteoformSequences.fasta"),
                            proteoformSequences
                    );
                    pb.step();
                }
            }
        }
        /*
        Write genotype sequences MSA.
        for (FeatureEntry featureEntry : variantsDictionary.features.values()) {
            Logging.logStatus("Write " + featureEntry.name + " MSA");
            ConcurrentSkipListMap<String, HashMap<String, Character>> proteoformSequences =
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
            String pfContent;
            String pfPosInfo;
            String posKey;
            int insPos;
            int insPosTotal;
            int pfVariantStart;
            List<Integer> pfInsPos;
            TreeSet<String> faEntryIds = new TreeSet<>();
            for (Map.Entry<String, String> chainEntry : featureEntry.allocatedProtein.chainSequences.entrySet()) {
                sequenceChars = chainEntry.getValue().toCharArray();
                for (int i = 0; i < sequenceChars.length; i++) {
                    if (!proteoformSequences.containsKey((i + 1) + "+0")) {
                        proteoformSequences.put((i + 1) + "+0", new HashMap<>());
                    }
                    proteoformSequences.get((i + 1) + "+0").put("Chain" + chainEntry.getKey(),
                            Character.isLowerCase(sequenceChars[i]) ? '.' : sequenceChars[i]);
                }
                faEntryIds.add("Chain" + chainEntry.getKey());
            }
            for (ProteoformEntry proteoform : featureEntry.allocatedProtein.proteoforms.values()) {
                if (proteoform.annotations.containsKey("PT") && proteoform.annotations.get("PT").equals("true")) {
                    // Skip PT proteoform.
                    continue;
                }
                if (proteoform.name.equals(AllocatedProteinEntry.WILD_TYPE_PROTEOFORM_ID)) {
                    sequenceChars =
                            Bio.translateNucSequence(featureEntry.nucleotideSequence, true, true, featureEntry.isSense).toCharArray();
                    for (int i = 0; i < sequenceChars.length; i++) {
                        if (!proteoformSequences.containsKey((i + 1) + "+0")) {
                            proteoformSequences.put((i + 1) + "+0", new HashMap<>());
                        }
                        proteoformSequences.get((i + 1) + "+0").put("WildType", sequenceChars[i]);
                    }
                    faEntryIds.add("WildType");
                } else {
                    String[] proteoformVariants = proteoform.annotations.get("vSwab").split("\\|");
                    for (String proteoformVariant : proteoformVariants) {
                        pfContent = proteoformVariant.split("@")[0];
                        pfContent = pfContent.replace(Bio.GAP, Bio.DELETION_AA1);
                        pfPosInfo = proteoformVariant.split("@")[1];
                        pfVariantStart = pfPosInfo.contains("+") ? Integer.parseInt(pfPosInfo.split("\\+")[0]) : Integer.parseInt(pfPosInfo);
                        pfInsPos = pfPosInfo.contains("+") ? Arrays.stream(pfPosInfo.split("\\+")[1].split(",")).map(Integer::valueOf).collect(Collectors.toList()) : new ArrayList<>();
                        sequenceChars = pfContent.toCharArray();
                        insPos = 0;
                        insPosTotal = 0;
                        for (int i = 0; i < sequenceChars.length; i++) {
                            if (pfInsPos.contains(i)) {
                                insPos += 1;
                                insPosTotal += 1;
                                posKey = (pfVariantStart + i - insPosTotal) + "+" + insPos;
                            } else {
                                insPos = 0;
                                posKey = (pfVariantStart + i - insPosTotal) + "+0";
                            }
                            if (!proteoformSequences.containsKey(posKey)) {
                                proteoformSequences.put(posKey, new HashMap<>());
                            }
                            proteoformSequences.get(posKey).put(proteoform.name, sequenceChars[i]);
                        }
                    }
                    faEntryIds.add(proteoform.name);
                }
            }
            HashMap<String, ArrayList<String>> faSequences = new HashMap<>();
            ArrayList<String> faHeadersList;
            StringBuilder faSequenceBuilder = new StringBuilder();
            for (String faEntryId : faEntryIds) {
                faHeadersList = new ArrayList<>();
                faSequenceBuilder.setLength(0);
                for (String pos : proteoformSequences.keySet()) {
                    if (proteoformSequences.get(pos).containsKey(faEntryId)) {
                        faSequenceBuilder.append(proteoformSequences.get(pos).get(faEntryId));
                    } else {
                        if (proteoformSequences.get(pos).containsKey("WildType")) {
                            faSequenceBuilder.append(proteoformSequences.get(pos).get("WildType"));
                        } else {
                            faSequenceBuilder.append(Bio.GAP);
                        }
                    }
                }
                faHeadersList.add(faEntryId);
                if (!faEntryId.startsWith("Chain")) {
                    if (faEntryId.equals("WildType")) {
                        faHeadersList.addAll(featureEntry.allocatedProtein.proteoforms.get(AllocatedProteinEntry.WILD_TYPE_PROTEOFORM_ID).samples);
                    } else {
                        faHeadersList.addAll(featureEntry.allocatedProtein.proteoforms.get(faEntryId).samples);
                    }
                }
                String sequence = faSequenceBuilder.toString();
                if (faSequences.containsKey(sequence)) {
                    faSequences.get(sequence).add(faEntryId);
                } else {
                    faSequences.put(sequence, faHeadersList);
                }
            }
            IO.writeFasta(
                    new File(cliarguments.outputFile.getParent() + "/" + featureEntry.chromosome + "_" + featureEntry.name
                            + ".fasta"), faSequences);
        }
         */
    }

}