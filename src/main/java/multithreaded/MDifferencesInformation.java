/**
 * 
 */
package multithreaded;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import io.myVCFFileReader;

/**
 * @author Alexander Seitz
 *
 */
public class MDifferencesInformation extends AMultithreadedVCF<MDifferencesInformation> {
	
	private List<Integer> positionsToReAnalyze;
	private Map<Integer, String> detailedInformation;
	String sampleName;
	
	public MDifferencesInformation(String sampleName, myVCFFileReader vcfReader, List<Integer> positionsToReAnalyze) {
		super(vcfReader);
		this.sampleName = sampleName;
		this.positionsToReAnalyze = positionsToReAnalyze;
		this.detailedInformation = new TreeMap<Integer, String>();
	}

	/* (non-Javadoc)
	 * @see java.util.concurrent.Callable#call()
	 */
	@Override
	public MDifferencesInformation call() {
		this.detailedInformation = this.vcfReader.reAnalyzePositions(this.positionsToReAnalyze);
		return this;
	}
	
	public Map<Integer, String> getDetailedInformation(){
		return this.detailedInformation;
	}
	
	public String getSampleName(){
		return this.sampleName;
	}

}
