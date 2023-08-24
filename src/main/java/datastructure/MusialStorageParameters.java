package datastructure;

/**
 * Stores parameters used for variant call filtering.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.1
 */
public class MusialStorageParameters {

    /**
     * The minimal coverage of a nucleotide variant to be accepted (extracted from the DP attribute of .vcf file format).
     */
    public final double minimalCoverage;
    /**
     * The minimal hom. frequency of a nucleotide variant to be accepted (extracted from the AF attribute of .vcf file format).
     */
    public final double minimalHomozygousFrequency;
    /**
     * The minimal het. frequency of a nucleotide variant to be accepted (extracted from the AF attribute of .vcf file format).
     */
    public final double minimalHeterozygousFrequency;
    /**
     * The maximal het. frequency of a nucleotide variant to be accepted (extracted from the AF attribute of .vcf file format).
     */
    public final double maximalHeterozygousFrequency;
    /**
     * The minimal quality of a nucleotide variant to be accepted (extracted from the Phred scaled GQ attribute of .vcf file format).
     */
    public final double minimalQuality;

    /**
     * Constructor of {@link MusialStorageParameters}.
     *
     * @param minimalCoverage              {@link Double}; Minimal read coverage to use for variant filtering.
     * @param minimalFrequency             {@link Double}; Minimal allele frequency to use for hom. variant filtering.
     * @param minimalHeterozygousFrequency {@link Double}; Minimal allele frequency to use for het. variant filtering.
     * @param maximalHeterozygousFrequency {@link Double}; Maximal allele frequency to use for het. variant filtering.
     * @param minimalQuality               {@link Double}; Minimal Phred scaled genotyping call quality to use for variant filtering.
     */
    public MusialStorageParameters(Double minimalCoverage, Double minimalFrequency, Double minimalHeterozygousFrequency, Double maximalHeterozygousFrequency, Double minimalQuality) {
        this.minimalCoverage = minimalCoverage;
        this.minimalHomozygousFrequency = minimalFrequency;
        this.minimalHeterozygousFrequency = minimalHeterozygousFrequency;
        this.maximalHeterozygousFrequency = maximalHeterozygousFrequency;
        this.minimalQuality = minimalQuality;
    }
}
