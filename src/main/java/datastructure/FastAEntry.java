/**
 * 
 */
package datastructure;

/**
 * @author seitz
 *
 */
public class FastAEntry {
	
	private String header = "";
	private String sequence = "";
	
	public FastAEntry(String header, String sequence){
		this.header = header;
		this.sequence = sequence;
	}
	
	public String getIthChar(Integer i){
		return ""+this.sequence.charAt(i);
	}
	
	public String getHeader() {
		return header;
	}
	public void setHeader(String header) {
		this.header = header;
	}
	public String getSequence() {
		return sequence;
	}
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	public int getSequenceLength(){
		return this.sequence.length();
	}
	
	public String toString(){
		StringBuffer result = new StringBuffer();
		result.append(">");
		result.append(this.header);
		result.append("\n");
		String[] seq = this.sequence.split("(?<=\\G.{80})");
		for(int i=0; i<seq.length; i++){
			result.append(seq[i]);
			result.append("\n");
		}
		return result.toString();
	}
}
