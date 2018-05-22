/**
 * 
 */
package datastructure;

/**
 * @author Alexander Seitz
 *
 */
public class SNPStatistics {
	
	private String SampleName;
	private Long numSNPs;
	private Long numHetSNPs;
	private Long totalCalls;
	private Long refCalls;
	private Long noCalls;
	private Double foldCoverage;
	private Double percentCoverage;
	/**
	 * @param sampleName
	 * @param numSNPs
	 * @param numHetSNPs
	 * @param refLength
	 * @param refCalls
	 * @param noCalls
	 * @param foldCoverage
	 * @param percentCoverage
	 */
	public SNPStatistics(String sampleName, Long numSNPs, Long numHetSNPs, Long totalCalls, Long refCalls, Long noCalls, Double foldCoverage, Double percentCoverage) {
		SampleName = sampleName;
		this.numSNPs = numSNPs;
		this.numHetSNPs = numHetSNPs;
		this.totalCalls = totalCalls;
		this.refCalls = refCalls;
		this.noCalls = noCalls;
		this.foldCoverage = foldCoverage;
		this.percentCoverage = percentCoverage;
	}
	
	public String generateLine() {
		StringBuffer result = new StringBuffer();
		String separator = "\t";
		result.append(this.SampleName);
		result.append(separator);
		result.append(this.numSNPs);
		result.append(separator);
		result.append(this.numHetSNPs);
		result.append(separator);
		result.append(this.foldCoverage);
		result.append(separator);
		result.append(this.percentCoverage);
		result.append(separator);
		result.append(this.refCalls);
		result.append(separator);
		result.append(this.totalCalls);
		result.append(separator);
		result.append(this.noCalls);
		return result.toString();
	}

}
