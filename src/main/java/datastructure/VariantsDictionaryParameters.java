package datastructure;

/**
 * Stores the main parameters used for variant filtering for a {@link VariantsDictionary}.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class VariantsDictionaryParameters {

    /**
     * The minimal coverage of a variant to be accepted (extracted from the DP attribute of .vcf file format).
     */
    public final double minCoverage;
    /**
     * The minimal hom. frequency of a variant to be accepted (extracted from the AF attribute of .vcf file format).
     */
    public final double minFrequency;
    /**
     * The minimal het. frequency of a variant to be accepted (extracted from the AF attribute of .vcf file format).
     */
    public final double minHetFrequency;
    /**
     * The maximal hem. frequency of a variant to be accepted (extracted from the AF attribute of .vcf file format).
     */
    public final double maxHetFrequency;
    /**
     * The minimal quality a variant to be accepted (extracted from the Phred scaled GQ attribute of .vcf file format).
     */
    public final double minQuality;

    /**
     * Constructor of {@link VariantsDictionaryParameters}.
     *
     * @param minCoverage  {@link Double}; Minimal read coverage to use for variant filtering.
     * @param minFrequency {@link Double}; Minimal allele frequency to use for hom. variant filtering.
     * @param minHet       {@link Double}; Minimal allele frequency to use for het. variant filtering.
     * @param maxHet       {@link Double}; Maximal allele frequency to use for het. variant filtering.
     * @param minQuality   {@link Double}; Minimal Phred scaled genotyping call quality to use for variant filtering.
     */
    public VariantsDictionaryParameters(Double minCoverage, Double minFrequency, Double minHet, Double maxHet, Double minQuality) {
        this.minCoverage = minCoverage;
        this.minFrequency = minFrequency;
        this.minHetFrequency = minHet;
        this.maxHetFrequency = maxHet;
        this.minQuality = minQuality;
    }
}
