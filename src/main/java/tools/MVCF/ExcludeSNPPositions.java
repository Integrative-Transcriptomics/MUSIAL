/**
 * 
 */
package tools.MVCF;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * @author Alexander Seitz
 *
 */
public class ExcludeSNPPositions {
	
	private String mvcfFolder = "";
	private String alignmentFile = "";
	private String snpTableFile = "";
	private List<Integer> toExclude = new LinkedList<Integer>();

	private List<String> header;
	private HashMap<String,String> sequences;
	private List<Integer> snpsPos;
	private HashMap<Integer,Integer> snpPosNum;
	private HashMap<String, List<String>> snps;
	
	public ExcludeSNPPositions(String mvcfFolder, String toExcludeFile) {
		List<Integer> lis = new LinkedList<Integer>();
		try {
			@SuppressWarnings("resource")
			BufferedReader br = new BufferedReader(new FileReader(toExcludeFile));
			String line = "";
			while((line = br.readLine()) != null){
				Integer i = Integer.parseInt(line.trim());
				lis.add(i);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		new ExcludeSNPPositions(mvcfFolder, lis);
	}
	
	public ExcludeSNPPositions(String mvcfFolder, List<Integer> toExclude) {
		File folder = new File(mvcfFolder);
		if(!folder.exists()){
			System.err.println("MVCF-Folder does not exist:\n"+mvcfFolder);
		}
		this.mvcfFolder = folder.getAbsolutePath();
		this.toExclude = toExclude;
		this.header = new LinkedList<String>();
		this.sequences = new HashMap<String,String>();
		this.snpsPos = new ArrayList<Integer>();
		this.snpPosNum = new HashMap<Integer, Integer>();
		this.snps = new HashMap<String, List<String>>();
		analyzeMVCFFolder();
		run();
	}

	/**
	 * 
	 */
	private void run() {
		readAlignment();
		readTable();
		excludePositions();
		writeNewAlignment();
	}

	private void writeNewAlignment() {
		String outFile = this.alignmentFile.replace(".fasta", "");
		outFile = outFile + "_filtered.fasta";
		try {
			PrintWriter outMerged = new PrintWriter(new BufferedWriter(new FileWriter(outFile, false)));
			for(String head: this.header){
				String seq = this.sequences.get(head);
				String[] formattedSeq = seq.split("(?<=\\G.{80})");
				outMerged.println(">"+head);
				for(int i=0; i<formattedSeq.length; i++){
					outMerged.println(formattedSeq[i]);
				}
				outMerged.flush();
			}
			outMerged.flush();
			outMerged.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void excludePositions() {
		Set<Integer> removed = new HashSet<Integer>();
		List<Integer> excludePositions = new LinkedList<Integer>();
		List<Integer> notFound = new LinkedList<Integer>();
		for(Integer ex: this.toExclude){
			if(this.snpPosNum.containsKey(ex) && !excludePositions.contains(ex)){
				excludePositions.add(this.snpPosNum.get(ex));
			}else{
				notFound.add(ex);
			}
		}
		for(String head: this.header){
			String seq = this.sequences.get(head);
			StringBuilder newSeq = new StringBuilder();
			for(int j=0; j<seq.length(); j++){
				char c = seq.charAt(j);
				if(!excludePositions.contains(j+1)){
					newSeq.append(c);
				}else{
					removed.add(j+1);
				}
			}
			this.sequences.put(head, newSeq.toString());
		}
	}

	private void readTable() {
		try {
			@SuppressWarnings("resource")
			BufferedReader br = new BufferedReader(new FileReader(this.snpTableFile));
			String line = "";
			Integer lineNum = 0;
			String[] names = new String[0];
			while((line = br.readLine()) != null){
				String[] splitted = line.split("\t");
				if(lineNum>0){
					Integer SNP = Integer.parseInt(splitted[0]);
					this.snpPosNum.put(SNP, lineNum);
					this.snpsPos.add(SNP);
					for(int i=0; i<names.length; i++){
						String sample = names[i];
						String snp = splitted[i+1];
						if(this.snps.containsKey(sample)){
							List<String> lis = this.snps.get(sample);
							lis.add(snp);
							this.snps.put(sample, lis);
						}else{
							List<String> lis = new ArrayList<String>();
							lis.add(snp);
							this.snps.put(sample, lis);
						}
					}
				}else{
					names = Arrays.copyOfRange(splitted, 1, splitted.length);
				}
				lineNum++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void readAlignment() {
		try {
			@SuppressWarnings("resource")
			BufferedReader br = new BufferedReader(new FileReader(this.alignmentFile));
			String line = "";
			String currHeader = "";
			StringBuffer currSeq = new StringBuffer();
			while((line = br.readLine()) != null){
				if(line.startsWith(";")){
					continue;
				}else if(line.startsWith(">")){
					if("".equals(currHeader)){
						currHeader = line.substring(1);
					}else{
						this.header.add(currHeader);
						this.sequences.put(currHeader, currSeq.toString());
						currHeader = line.substring(1);
						currSeq = new StringBuffer();
					}
				}else{
					currSeq.append(line);
				}
			}
			this.header.add(currHeader);
			this.sequences.put(currHeader, currSeq.toString());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * 
	 * 
	 */
	private void analyzeMVCFFolder() {
		this.alignmentFile = this.mvcfFolder+"/snpAlignment.fasta";
		this.snpTableFile = this.mvcfFolder+"/snpTable.tsv";
		if(!new File(this.alignmentFile).exists()){
			System.err.println("snpAlignment.fasta does not exist in folder:\n"+this.mvcfFolder);
		}
		if(!new File(this.snpTableFile).exists()){
			System.err.println("snpTable.tsv does not exist in folder:\n"+this.mvcfFolder);
		}
	}
	
	public static void main(String[] args) {
		List<Integer> exclude = new LinkedList<Integer>();
		exclude.add(73);
		new ExcludeSNPPositions("/share/projects/Lepra/results/Verena/GC/EAGER_Results_ALN_01_DD/EAGER_Results/MVCF/VCFMulti_20160926-141225", exclude);
	}
}
