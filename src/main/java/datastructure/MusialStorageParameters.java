package datastructure;

/**
 * Stores parameters used for variant call filtering.
 *
 * @author Simon Hackl
 * @version 2.3
 * @since 2.1
 */
public class MusialStorageParameters {

    /**
     * The minimal coverage of an allele to be accepted (extracted from the DP attribute of .vcf file format).
     */
    public final double minimalCoverage;
    /**
     * The minimal frequency of an allele to be accepted (extracted from the AF attribute of .vcf file format).
     */
    public final double minimalFrequency;


    /**
     * Constructor of {@link MusialStorageParameters}.
     *
     * @param minimalCoverage  {@link Double}; Minimal read coverage to use for allele filtering.
     * @param minimalFrequency {@link Double}; Minimal allele frequency to use for allele filtering.
     */
    public MusialStorageParameters(Double minimalCoverage, Double minimalFrequency) {
        this.minimalCoverage = minimalCoverage;
        this.minimalFrequency = minimalFrequency;
    }
}
