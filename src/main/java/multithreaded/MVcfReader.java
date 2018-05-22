/**
 * 
 */
package multithreaded;

import io.myVCFFileReader;
import tools.Analyzer;

/**
 * @author Alexander Seitz
 *
 */
public class MVcfReader extends AMultithreadedVCF<MVcfReader> {
	
	Analyzer analyzer;

	/**
	 * 
	 */
	public MVcfReader(Analyzer analyzer, myVCFFileReader vcfReader) {
		super(vcfReader);
		this.analyzer = analyzer;
	}

	/* (non-Javadoc)
	 * @see java.util.concurrent.Callable#call()
	 */
	@Override
	public MVcfReader call() {
		long start = System.currentTimeMillis();
		try {
			this.vcfReader.parseVCFFile();
		}catch(Exception e) {
			e.printStackTrace();
		}
		long end = System.currentTimeMillis();
		long timeInS = (end - start) / 1000;
		analyzer.finished(this.vcfReader.getSampleName(), timeInS);
		return this;
	}

}
