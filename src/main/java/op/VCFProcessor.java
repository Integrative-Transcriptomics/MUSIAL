package op;

import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Tuple;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFUtils;
import model.Sample;
import model.Storage;
import model.VariantCall;
import org.apache.commons.lang3.StringUtils;
import util.Bio;
import util.Constants;
import util.Logging;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;

/**
 * The {@code VCFProcessor} class is responsible for processing Variant Call Format (VCF) files.
 * <p>
 * This class handles the analysis of VCF files, including the extraction of variant data, imputation of contigs, and integration of the
 * processed data into the storage system. It provides methods to analyze VCF files, process variant contexts, and track statistics such as
 * the number of processed, ignored, and filtered variant calls.
 */
public class VCFProcessor {

    /**
     * The storage object for managing processed genomic data.
     */
    private final Storage storage;

    /**
     * A list of file paths to the VCF files to be processed.
     * <p>
     * This field stores the paths to the Variant Call Format (VCF) files that will be read and processed.
     */
    private final List<Path> paths;

    /**
     * A flag indicating whether to impute contigs from the VCF files.
     * <p>
     * If set to {@code true}, the program will infer and add contigs to the storage based on the contigs present in the VCF files. This is
     * useful when the contigs are not explicitly defined in the input data.
     */
    private final boolean imputeContigs;

    /**
     * Tracks the total number of processed variant calls.
     * <p>
     * This field is used to count the number of variant calls that have been processed during the execution of the program. It is
     * incremented each time a variant call is processed, regardless of its outcome.
     */
    private long processedCallsCount = 0;

    /**
     * Tracks the total number of ignored variant calls.
     * <p>
     * This field is used to count the number of variant calls that were ignored during processing. A call may be ignored for various
     * reasons, such as missing data or being classified as a reference call.
     */
    private long ignoredCallsCount = 0;

    /**
     * Tracks the total number of filtered variant calls.
     * <p>
     * This field is used to count the number of variant calls that were filtered out during processing. Filtering may occur due to low
     * coverage, low frequency, or other criteria defined in the program.
     */
    private long filteredCallsCount = 0;

    /**
     * Constructs a new {@link VCFProcessor} instance for processing VCF files.
     * <p>
     * This constructor initializes the processor with the specified list of VCF file paths, a storage object for managing genomic data, and
     * a flag indicating whether to impute contigs from the VCF files. The {@code imputeContigs} flag determines if contigs should be
     * inferred and added to the storage based on the VCF data.
     *
     * @param paths         A {@link List} of {@link Path} objects representing the file paths to the VCF files to be processed.
     * @param storage       The {@link Storage} object used to store and manage processed genomic data.
     * @param imputeContigs A {@code boolean} flag indicating whether to infer and add contigs from the VCF files to the storage.
     */
    public VCFProcessor(List<Path> paths, Storage storage, boolean imputeContigs) {
        this.paths = paths;
        this.storage = storage;
        this.imputeContigs = imputeContigs;
    }

    /**
     * Analyzes the provided VCF files and processes their variant data.
     * <p>
     * This method iterates through the list of VCF file paths, creates temporary indexed VCF files, and processes their variant contexts.
     * If the `imputeContigs` flag is set, it infers contigs from the VCF files and adds them to the storage. The method processes variant
     * contexts either for specific features in the storage or for all variants if no features are defined.
     *
     * @throws IOException If an I/O error occurs during file operations or VCF processing.
     */
    public void analyzeFiles() throws IOException {
        for (Path path : paths) {
            // Create temporary indexed VCF files for processing.
            File file = VCFUtils.createTemporaryIndexedVcfFromInput(path.toFile(), String.valueOf(path.hashCode()));
            File temporaryIndexedVcfFile = VCFUtils.createTemporaryIndexedVcfFromInput(file, String.valueOf(file.hashCode()));
            try (VCFFileReader vcfFileReader = new VCFFileReader(temporaryIndexedVcfFile)) {
                if (imputeContigs) {
                    // Extract unique intervals from the VCF file.
                    IntervalList intervalList = vcfFileReader.toIntervalList().uniqued();
                    Set<String> vcfContigs = new HashSet<>();
                    // Collect contig names from the intervals.
                    intervalList.getIntervals().forEach(interval -> vcfContigs.add(interval.getContig()));
                    // Add inferred contigs to the storage.
                    for (String contig : vcfContigs) {
                        storage.addContig(contig, Constants.EMPTY);
                    }
                }

                // Process variant contexts based on feature availability.
                if (!storage.getFeatures().isEmpty()) {
                    // Process variants for each feature in the storage.
                    storage.getFeatures().forEach(feature -> process(vcfFileReader.query(feature.contig, feature.start, feature.end),
                            path.toAbsolutePath()));
                } else {
                    // Process all variants if no features are defined.
                    process(vcfFileReader.iterator(), path.toAbsolutePath());
                }
            } finally {
                // Delete the temporary indexed VCF file.
                //noinspection ResultOfMethodCallIgnored
                temporaryIndexedVcfFile.delete();
            }
        }
    }

    /**
     * Processes a set of variant contexts from a VCF file and updates the storage with variant calls.
     * <p>
     * This method iterates through the provided {@link VariantContext} objects, processes each genotype, and extracts relevant information
     * such as reference and alternative alleles, allelic depths, and coverage data. The processed data is then added to the storage as
     * variant calls.
     *
     * @param variantContextIterator An {@link Iterator} of {@link VariantContext} objects representing the variants to process.
     * @param path                   The {@link Path} to the VCF file being processed, used for logging and error reporting.
     */
    private void process(Iterator<VariantContext> variantContextIterator, Path path) {
        while (variantContextIterator.hasNext()) {
            VariantContext variantContext = variantContextIterator.next();
            String contigIdentifier = variantContext.getContig();
            int position = variantContext.getStart();

            // Skip positions excluded by the build configuration.
            if (storage.parameters.isPositionMasked(contigIdentifier, variantContext.getStart())) continue;

            // Process each genotype in the current VariantContext.
            for (Genotype genotype : variantContext.getGenotypes()) {
                processedCallsCount++; // Increment the count of processed genotype records.

                if (genotype.isNoCall()) {
                    ignoredCallsCount++; // Skip and count no-call genotypes.
                    continue;
                }

                // Log a warning if AD and DP attributes are missing.
                if (!(genotype.hasDP() && (genotype.hasAD() || genotype.hasAnyAttribute("COV")))) {
                    Logging.logWarningOnce("MISSING_AD_DP_ATTRIBUTES",
                            String.format("Some variants may be ignored as AD/COV/DP4 and DP attributes are unavailable. Detected at" +
                                            " site %s %d for sample %s in file %s.",
                                    contigIdentifier, variantContext.getStart(), genotype.getSampleName(), path));
                }

                // Extract the sample identifier and ensure the sample exists in storage.
                String sampleIdentifier = genotype.getSampleName().split("\\$")[0];
                Sample sample = storage.addSample(sampleIdentifier);

                // Extract allelic depth (AD) information for the genotype.
                int[] ADs;
                if (genotype.hasAD()) {
                    ADs = genotype.getAD();
                } else if (genotype.hasAnyAttribute("COV")) {
                    ADs = (int[]) genotype.getAnyAttribute("COV");
                } else if (genotype.getAlleles().size() == 2
                        && variantContext.getNSamples() == 1
                        && variantContext.hasAttribute("DP4")) {
                    // Use DP4 attribute as a fallback if AD is missing.
                    List<Integer> DP4 = variantContext.getAttributeAsIntList("DP4", 0);
                    ADs = new int[]{
                            DP4.get(0) + DP4.get(1), // Reference allele coverage.
                            DP4.get(2) + DP4.get(3)  // Alternative allele coverage.
                    };
                } else {
                    ignoredCallsCount++; // Skip genotypes with insufficient information.
                    continue;
                }

                // Compute the total depth of coverage (ADSum) from allelic depths.
                int ADSum = Arrays.stream(ADs).sum();

                // Log warnings for discrepancies between AD sum and DP.
                if (genotype.hasDP()) {
                    if (ADSum > genotype.getDP()) {
                        Logging.logWarningOnce("AD_SUM_GREATER_THAN_DP",
                                String.format("Possible error in genotype data. Summed allelic depth (%d) is greater than total depth" +
                                                " (%d) at site %s %d for sample %s in file %s.",
                                        ADSum, genotype.getDP(), contigIdentifier, variantContext.getStart(), sampleIdentifier,
                                        path));
                    } else if (ADSum < 0.5 * genotype.getDP()) {
                        Logging.logWarningOnce("AD_SUM_LOWER_THAN_DP",
                                String.format("Allegedly low-quality data. Summed allelic depth (%d) is much lower than total depth " +
                                                "(%d) at site %s %d for sample %s in file %s.",
                                        ADSum, genotype.getDP(), contigIdentifier, variantContext.getStart(), sampleIdentifier,
                                        path));
                    }
                }

                // Process reference and alternative alleles for the genotype.
                String REF;
                String ALT;
                int AD;

                // Create a list to store alternatives for the current genotype.
                List<VariantCall.CallAlternative> alternatives = new ArrayList<>(ADs.length);

                for (int i = 0; i < ADs.length; i++) {
                    AD = ADs[i]; // Retrieve the allelic depth for the current allele.

                    if (i == 0) {
                        // For the reference allele, set REF to the first base and ALT to a dot.
                        REF = variantContext.getReference().getBaseString().substring(0, 1);
                        ALT = Constants.DOT;
                    } else {
                        // Skip alternative alleles with zero depth.
                        if (AD == 0) continue;

                        // For alternative alleles, retrieve the full reference and alternative sequences.
                        REF = variantContext.getReference().getBaseString();
                        ALT = variantContext.getAlleles().get(i).getBaseString();

                        // Handle upstream-deletion cases and ensure canonical formatting.
                        if (ALT.equals("*")) {
                            REF = REF.substring(0, 1);
                        } else if (Bio.isCanonicalVariant(REF, ALT)) {
                            REF = Bio.padGaps(REF, ALT.length());
                            ALT = Bio.padGaps(ALT, REF.length());
                        } else {
                            String commonSuffix = StringUtils.reverse(
                                    StringUtils.getCommonPrefix(StringUtils.reverse(REF), StringUtils.reverse(ALT))
                            );
                            REF = StringUtils.removeEnd(REF, commonSuffix);
                            ALT = StringUtils.removeEnd(ALT, commonSuffix);

                            if (Bio.isCanonicalVariant(REF, ALT)) {
                                REF = Bio.padGaps(REF, ALT.length());
                                ALT = Bio.padGaps(ALT, REF.length());
                            } else {
                                Tuple<String, String> alignment = Bio.globalNucleotideSequenceAlignment(
                                        REF, ALT, 5, 2, true, false, 0
                                );
                                REF = alignment.a;
                                ALT = alignment.b;
                            }
                        }
                    }

                    // Add the alternative allele to the list.
                    alternatives.add(new VariantCall.CallAlternative(REF, ALT, AD));
                }

                // Add the variant call to the storage and handle the result.
                VariantCall.Flag flag = sample.addVariantCall(contigIdentifier, position, alternatives, storage.parameters, path);

                switch (flag) {
                    case PASS -> {
                        // Variant was added successfully.
                    }
                    case REFERENCE_CALL -> ignoredCallsCount++;
                    case LOW_COVERAGE, LOW_FREQUENCY, MISSING_UPSTREAM_DELETION -> filteredCallsCount++;
                    default -> throw new IllegalStateException("Unexpected value: " + flag);
                }
            }
        }
    }

    /**
     * Retrieves the total number of processed variant calls.
     *
     * @return The total number of processed variant calls as a {@code long}.
     */
    public long getProcessedCallsCount() {
        return processedCallsCount;
    }

    /**
     * Retrieves the total number of ignored variant calls.
     * <p>
     * This method returns the count of variant calls that were ignored during processing. A variant call may be ignored for reasons such as
     * missing data, being classified as a reference call, or lacking sufficient information for analysis.
     *
     * @return The total number of ignored variant calls as a {@code long}.
     */
    public long getIgnoredCallsCount() {
        return ignoredCallsCount;
    }

    /**
     * Retrieves the total number of filtered variant calls.
     * <p>
     * This method returns the count of variant calls that were filtered out during processing. Filtering may occur due to criteria such as
     * low coverage, low frequency, or other conditions defined in the program.
     *
     * @return The total number of filtered variant calls as a {@code long}.
     */
    public long getFilteredCallsCount() {
        return filteredCallsCount;
    }

}
