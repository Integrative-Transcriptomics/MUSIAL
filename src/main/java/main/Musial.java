package main;

import cli.CLIParameters;
import cli.CLIParametersUpdateVDict;
import components.SnpEffAnnotator;
import components.VariantsDictionaryFactory;
import datastructure.FastaContainer;
import datastructure.FeatureEntry;
import datastructure.ProteoformEntry;
import datastructure.SampleEntry;
import datastructure.VariantsDictionary;
import exceptions.MusialBioException;
import exceptions.MusialCLException;
import exceptions.MusialIOException;
import exceptions.MusialIntegrityException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
import org.checkerframework.checker.units.qual.A;
import org.javatuples.Triplet;
import runnables.SampleAnalyzerRunnable;
import utility.Bio;
import utility.IO;
import utility.Logging;

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
     * Logo `String` of the software (MUSIAL ASII-art).
     */
    public static String LOGO = "";
    /**
     * Name `String` of the software (MUSIAL).
     */
    public static String NAME = "MUSIAL";
    /**
     * Version of the software.
     */
    public static String VERSION = "";
    /**
     * Contact of the software.
     */
    public static String CONTACT = "";
    /**
     * License of the software.
     */
    public static String LICENSE = "";
    /**
     * Which module of MUSIAL is run.
     */
    public static MusialModules MODULE = null;
    /**
     * Original out stream.
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
     * TODO
     *
     * @param args Arguments parsed from the command line.
     */
    public static void main(String[] args) {
        try {
            if (args.length == 0) {
                loadMetadata();
                System.out.println("\nNo module specified.");
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
                System.out.println("\nUnknown module " + args[0] + ".");
                printModuleInformation();
            } else {
                loadMetadata();
            }
            CLIParameters arguments = parseCLIArguments(args);
            // 2. Delegate pipeline to execute dependant on chosen module.
            //noinspection SwitchStatementWithTooFewBranches
            switch (MODULE) {
                case updateVDict -> runUpdateVDict((CLIParametersUpdateVDict) arguments);
            }
        } catch (Exception e) {
            if (DEBUG) {
                e.printStackTrace();
            }
            Logging.logError(e.getMessage());
        }
    }

    /**
     * Prints information about available modules to the user and exits.
     */
    private static void printModuleInformation() {
        System.out.println("\n" + "~".repeat(2) + " AVAILABLE MODULES " + "~".repeat(2));
        System.out.println("buildDB : TODO");
        System.exit(0);
    }

    /**
     * Loads metadata, such as the projects title and version from resources and prints the information to stdout.
     *
     * @throws IOException If any resource can not be loaded.
     */
    private static void loadMetadata() throws IOException {
        Properties properties = new Properties();
        // Load title from resources.
        InputStream in = Musial.class.getResourceAsStream("/info.properties");
        properties.load(in);
        Musial.LOGO = properties.getProperty("logo");
        Musial.VERSION = properties.getProperty("version");
        Musial.CONTACT = properties.getProperty("contact");
        Musial.LICENSE = properties.getProperty("license");
        // Print information into stdout.
        System.out.println(Musial.LOGO + Musial.VERSION);
        System.out.println(Musial.LICENSE + ", Contact: " + Musial.CONTACT + "\n");
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
     * @throws MusialCLException Catches exceptions from {@link CLIParametersUpdateVDict} class.
     * @throws MusialIOException Catches exceptions from {@link CLIParametersUpdateVDict} class.
     */
    private static CLIParameters parseCLIArguments(String[] args)
            throws MusialCLException, MusialIOException, FileNotFoundException, MusialBioException {
        CLIParameters cliParameters = null;
        //noinspection SwitchStatementWithTooFewBranches
        switch (MODULE) {
            case updateVDict -> cliParameters = new CLIParametersUpdateVDict(args);
        }
        assert cliParameters != null;
        return cliParameters;
    }

    private static void runUpdateVDict(CLIParametersUpdateVDict cliarguments)
            throws InterruptedException, MusialIOException, MusialIntegrityException, IOException, MusialBioException {
        // Read-in existing variants dictionary or build new one.
        VariantsDictionary variantsDictionary = VariantsDictionaryFactory.build(cliarguments);

        // Store/Update feature entry information.
        Set<String> featureIdsToUpdate;
        Set<String> featureIdsToRemove;
        try (ProgressBar pb = buildProgress()) {
            Logging.logStatus("Updating feature information.");
            featureIdsToUpdate = cliarguments.features.keySet();
            pb.maxHint(featureIdsToUpdate.size());
      /* TODO: Remove variants that only occur in this feature.
      for (String fId : featureIdsToRemove) {
        variantsDictionary.features.remove(fId);
        pb.step();
      }
       */
            FastaContainer referenceChromosomeFastaContainer = null;
            for (FastaContainer fC : IO.readFastaToSet(cliarguments.referenceFASTA)) {
                String chromosome = fC.getHeader().split(" ")[0].trim();
                if (variantsDictionary.chromosome.equals(chromosome)) {
                    referenceChromosomeFastaContainer = fC;
                    break;
                }
            }
            assert referenceChromosomeFastaContainer != null;
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

        // Store/Update sample variants information.
        // TODO: Remove samples.
        Logging.logStatus("Updating variants from samples.");
        try (ProgressBar pb = buildProgress()) {
            pb.maxHint(cliarguments.samples.size());
            for (SampleEntry sampleEntry : cliarguments.samples.values()) {
                SampleEntry.imputeVCFFileReader(sampleEntry);
                variantsDictionary.samples.put(sampleEntry.name, sampleEntry);
                pb.step();
            }
            pb.setExtraMessage(Logging.getDoneMessage());
        }
        try (ProgressBar pb = buildProgress()) {
            ExecutorService executor = Executors.newFixedThreadPool(Musial.THREADS);
            pb.maxHint((long) featureIdsToUpdate.size() * cliarguments.samples.size());
            for (String fId : featureIdsToUpdate) {
                for (SampleEntry sampleEntry : cliarguments.samples.values()) {
                    executor.execute(
                            new SampleAnalyzerRunnable(sampleEntry, cliarguments.features.get(fId), variantsDictionary, pb)
                    );
                }
            }
            executor.shutdown();
            //noinspection ResultOfMethodCallIgnored
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            pb.setExtraMessage(Logging.getDoneMessage());
        }
        File novelVariantsOut = new File("./temp/novelVariants.vcf");
        File tempDir = new File("./temp/");
        IO.generateDirectory(tempDir);
        try {
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
                try {
                    variantsDictionary.addVariantAnnotation(
                            Integer.parseInt(position),
                            altContent,
                            SnpEffAnnotator.convertAnnotation(lineFields[7].split(";")[0].split("=")[1])
                    );
                } catch (Exception e) {
                    continue;
                }
            }
        } finally {
            IO.deleteDirectory(tempDir);
        }

        // Infer proteoforms.
        ArrayList<Triplet<String, String, ArrayList<String>>> variableSegments;
        try (ProgressBar pb = buildProgress()) {
            Logging.logStatus("Infer protein variants.");
            pb.maxHint((long) variantsDictionary.features.size() * variantsDictionary.samples.size());
            for (String fId : variantsDictionary.features.keySet()) {
                if (variantsDictionary.features.get(fId).pdbFile == null) {
                    continue;
                }
                for (String sId : variantsDictionary.samples.keySet()) {
                    pb.setExtraMessage("Sample: " + sId + ", Feature: " + fId);
                    // FIXME: Bug if no protein is allocated, null pointer exception is thrown -> Validate that .pdb exists.
                    variableSegments = Bio.inferProteoform(variantsDictionary, fId, sId);
                    variantsDictionary.features.get(fId).allocatedProtein
                            .addProteoform(variantsDictionary.features.get(fId), sId, variableSegments);
                    pb.step();
                }
            }
        }

        // TODO: Migrate to separate module. Build MSA.
        /*
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
                if (proteoform.name.equals("WildTypex0.00")) {
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
                        faHeadersList.addAll(featureEntry.allocatedProtein.proteoforms.get("WildTypex0.00").samples);
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

        // Write built database to file.
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

}