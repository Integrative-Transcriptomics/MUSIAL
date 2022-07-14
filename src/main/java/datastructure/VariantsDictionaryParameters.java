package datastructure;

public class VariantsDictionaryParameters {

    public final double minCoverage;

    public final double minFrequency;

    public final double minHetFrequency;

    public final double maxHetFrequency;

    public final double minQuality;

    public VariantsDictionaryParameters(Double minCoverage, Double minFrequency, Double minHet, Double maxHet, Double minQuality) {
        this.minCoverage = minCoverage;
        this.minFrequency = minFrequency;
        this.minHetFrequency = minHet;
        this.maxHetFrequency = maxHet;
        this.minQuality = minQuality;
    }
}
