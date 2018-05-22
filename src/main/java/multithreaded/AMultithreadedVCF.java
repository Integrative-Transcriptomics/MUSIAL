/**
 * 
 */
package multithreaded;

import io.myVCFFileReader;

/**
 * @author Alexander Seitz
 *
 */
public abstract class AMultithreadedVCF<T> implements IMultithreaded<T>{
	
	protected myVCFFileReader vcfReader;
	

	/**
	 * @param vcfReader2
	 */
	public AMultithreadedVCF(myVCFFileReader vcfReader) {
		this.vcfReader = vcfReader;
	}


	public myVCFFileReader getReader(){
		return this.vcfReader;
	}
	
}
