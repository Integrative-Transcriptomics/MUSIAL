package main;

import cli.CLIParameters;
import cli.CLIParametersUpdateVDict;
import components.Bio;
import components.IO;
import components.Logging;
import components.VariantsDictionaryFactory;
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
                System.out.println(">no module specified.");
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
                System.out.println("unknown module " + args[0] + ".");
                printModuleInformation();
            } else {
                loadMetadata();
            }
            CLIParameters arguments = parseCLIArguments(args);
            // 2. Invokes methods dependent on the user specified module.
            //noinspection SwitchStatementWithTooFewBranches
            switch (MODULE) {
                case updateVDict -> runUpdateVDict((CLIParametersUpdateVDict) arguments);
            }
        } catch (Exception e) {
            e.printStackTrace();
            Logging.logError(e.getMessage());
        }
    }

    /**
     * Prints information about available modules to the user and exits.
     */
    private static void printModuleInformation() {
        System.out.println("available MUSIAL modules:");
        System.out.println("\tupdateVDict : Generate a new or update an existing variants dictionary JSON file.");
        System.out.println("specify -h for any module to obtain more information.");
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
        //noinspection SwitchStatementWithTooFewBranches
        switch (MODULE) {
            case updateVDict -> cliParameters = new CLIParametersUpdateVDict(args);
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

}