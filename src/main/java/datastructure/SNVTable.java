/**
 * 
 */
package datastructure;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import io.myVCFFileReader;

/**
 * @author Alexander Seitz
 *
 */
public class SNVTable {

	private String inputFile;

	private Map<Integer,Integer> snpPosNum = new TreeMap<Integer, Integer>();
	private Map<String,List<String>> snps = new TreeMap<String, List<String>>();
	//	Map<Integer, Set<String>> deletions = new HashMap<Integer, Set<String>>();
	private List<Integer> snpPos = new ArrayList<Integer>();
	private List<String> sampleNames = new LinkedList<String>();
	private FastAEntry reference;

	/**
	 * 
	 */
	public SNVTable(String inputFile) {
		this.inputFile = inputFile;
		readSNPTable();
	}

	//	public SNPTable(String[] vcfFiles, String file){
	//		this.inputFile = file;
	//		for(String sample: vcfFiles){
	//			File vcf = new File(sample);
	//			String sampleName = vcf.getParent();
	//			
	//		}
	//	}

	/**
	 * @param table
	 * @param fastAEntry
	 */
	public SNVTable(Map<String, Map<Integer, String>> table, FastAEntry fastAEntry) {
		this.reference = fastAEntry;
		for(String sampleName : table.keySet()){
			this.sampleNames.add(sampleName);
			Map<Integer, String> snpMap = table.get(sampleName);
			List<Integer> snpPositions = new ArrayList<Integer>(snpMap.keySet());
			Collections.sort(snpPositions);
			if(this.snpPos.isEmpty()){
				this.snpPos = snpPositions;
				for(int i=1; i<=snpPositions.size(); i++){
					this.snpPosNum.put(snpPositions.get(i-1), i);
				}
			}
			List<String> snps = new ArrayList<String>();
			for(Integer snpP: snpPositions){
				snps.add(snpMap.get(snpP));
			}
			this.snps.put(sampleName, snps);
		}
		List<String> snps = new ArrayList<String>();
		for(Integer snpP: this.snpPos){
			snps.add(fastAEntry.getIthChar(snpP-1));
		}
		this.snps.put("Ref", snps);
		Set<Integer> toRemove = new TreeSet<Integer>();
		
		// remove entries that are no SNPs
		for(Integer pos: this.snpPos){
			Set<String> currSnpCalls = new HashSet<String>();
			String refSNP = getReferenceSnp(pos);
			for(String sample: this.sampleNames){
				String currSNP = getSnp(pos,sample);
				if(!"N".equals(currSNP) && !".".equals(currSNP) && !refSNP.equals(currSNP)) {
					currSnpCalls.add(currSNP);
				}
			}
			if(currSnpCalls.size()==0) {
//				System.out.println("Remove SNP: "+pos);
				toRemove.add(pos);
			}
		}
		for(Integer removePos: toRemove) {
			this.snpPos.remove(removePos);
			Integer removPosNum = this.snpPosNum.get(removePos);
			this.snpPosNum.remove(removePos);
			for(Integer snpPos: this.snpPosNum.keySet()) {
				Integer currSnpPosNum = this.snpPosNum.get(snpPos);
				if(currSnpPosNum > removPosNum) {
					this.snpPosNum.put(snpPos, currSnpPosNum-1);
				}
			}
			for(String sample: this.sampleNames) {
				List<String> currSnps = this.snps.get(sample);
				currSnps.remove((int)removPosNum-1);
				this.snps.put(sample, currSnps);
			}
			List<String> currSnps = this.snps.get("Ref");
			currSnps.remove((int)removPosNum-1);
			this.snps.put("Ref", currSnps);
		}
	}

	private void readSNPTable() {
		try {
			@SuppressWarnings("resource")
			BufferedReader br = new BufferedReader(new FileReader(this.inputFile));
			String line = "";
			Integer lineNum = 0;
			String[] names = new String[0];
			while((line = br.readLine()) != null){
				String[] splitted = line.split("\t");
				if(lineNum>0 && line.length()>0){
					Integer SNP = Integer.parseInt(splitted[0]);
					this.snpPosNum.put(SNP, lineNum);
					this.snpPos.add(SNP);
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
					for(String n: names){
						if(!"Ref".equals(n)){
							this.sampleNames.add(n);
						}
					}
				}
				lineNum++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * @return
	 */
	public List<Integer> getSNVs() {
		return this.snpPos;

	}

	/**
	 * @param snpPos
	 * @return
	 */
	public String getSnp(Integer snpPos, String sample) {
		if(this.snpPosNum.containsKey(snpPos)){
			Integer snpNum = this.snpPosNum.get(snpPos);
			if(this.snps.get(sample).get(snpNum-1) == null) {
//				System.out.println(sample+"\t"+snpPos+"\twas null");
				return "N";
			}
			return this.snps.get(sample).get(snpNum-1);
		}else{
			return "N";
		}
	}

	public String getReferenceSnp(Integer snpPos){
		return getSnp(snpPos, "Ref");
	}
	
	public String getReferencePos(Integer snpPos){
		if(snpPos==144242) {
			System.out.println("snpPos: "+this.reference==null);
		}
		if(this.reference != null) {
			return this.reference.getIthChar(snpPos-1);
		}
		return getReferenceSnp(snpPos);
	}

	/**
	 * @return
	 */
	public List<String> getSampleNames(){
		return this.sampleNames;
	}

	public String getRefCall(Integer pos){
		return this.getSnp(pos, "Ref");
	}

	/**
	 * @param name
	 * @param snps2
	 */
	public void add(String name, Map<Integer, String> newSnps) {
		if(!this.sampleNames.contains(name)){
			this.sampleNames.add(name);
			List<String> sampleSNPs = new ArrayList<String>();
			for(Integer pos: this.snpPos){
				if(newSnps.containsKey(pos)){
					sampleSNPs.add(newSnps.get(pos));
				}else{
					sampleSNPs.add("N");
				}
			}
			this.snps.put(name, sampleSNPs);
		}else{
			System.err.println("could not add sample to snpTable");
			System.err.println("sample name already exists: "+name);
		}
	}

	public void add(String name, myVCFFileReader newSnps) {
		if(!this.sampleNames.contains(name)){
			this.sampleNames.add(name);
			List<String> sampleSNPs = new ArrayList<String>();
			for(Integer pos: this.snpPos){
				sampleSNPs.add(newSnps.get(pos));
			}
			this.snps.put(name, sampleSNPs);
		}else{
			System.err.println("could not add sample to snpTable");
			System.err.println("sample name already exists: "+name);
		}
	}

	/**
	 * @param snpEffOutputFile
	 * @return
	 */
	private Map<Integer, Set<String>> parseSnpEffOutput(String snpEffOutputFile) {
		Map<Integer, Set<String>> result = new HashMap<Integer, Set<String>>();
		try {
			@SuppressWarnings("resource")
			BufferedReader br = new BufferedReader(new FileReader(snpEffOutputFile));
			String line = "";
			while((line = br.readLine()) != null) {
				if(line.startsWith("# Chromo")) {
					String[] splitted = line.split("\t");
					StringBuffer headerLine = new StringBuffer();
					for(int i=3; i<splitted.length; i++) {
						headerLine.append("\t");
						headerLine.append(splitted[i]);
					}
					Set<String> header = new TreeSet<String>();
					header.add(headerLine.toString());
					result.put(-1, header);
				}else if(line.startsWith("#")) {
					continue;
				}else {
					String[] splitted = line.split("\t");
					StringBuffer currLine = new StringBuffer();
					for(int i=3; i<splitted.length; i++) {
						currLine.append("\t");
						currLine.append(splitted[i]);
					}
					Integer pos = Integer.parseInt(splitted[1]);
					Set<String> set = new TreeSet<String>();
					if(result.containsKey(pos)) {
						set = result.get(pos);
					}
					set.add(currLine.toString());
					result.put(pos, set);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return result;
	}

	/**
	 * @param snpEffOutputFile
	 * @return
	 */
	public String generateSNPTableWithSNPEff(String snpEffOutputFile) {
		StringBuffer result = new StringBuffer();
		Map<Integer, Set<String>> snpEffOutput = parseSnpEffOutput(snpEffOutputFile);
		// write header
		result.append("Position");
		result.append("\t");
		result.append("Ref");
		for(String sample: this.sampleNames){
			result.append("\t");
			result.append(sample);
		}
		result.append(snpEffOutput.get(-1).iterator().next());
		result.append("\n");
		for(Integer pos: this.snpPos) {
			if(snpEffOutput.containsKey(pos)) {
				for(String snpEff: snpEffOutput.get(pos)) {
					result.append(pos);
					result.append("\t");
					result.append(getReferenceSnp(pos));
					for(String sample: this.sampleNames){
						result.append("\t");
						result.append(getSnp(pos,sample));
					}
					result.append(snpEff);
					result.append("\n");
				}
			}else {
				result.append(pos);
				result.append("\t");
				result.append(getReferenceSnp(pos));
				for(String sample: this.sampleNames){
					result.append("\t");
					result.append(getSnp(pos,sample));
				}
				result.append("\n");
			}
		}
		return result.toString();
	}

	/**
	 * @param snpEffOutputFile
	 * @return
	 */
	private Map<Integer, List<String>> parseSNPEffOutput(String snpEffOutputFile) {
		Map<Integer, List<String>> result = new TreeMap<Integer, List<String>>();
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(snpEffOutputFile));
			String line = "";
			while((line = br.readLine()) != null){
				if(line.startsWith("#")) {
					if(line.contains("ID=ANN")) {
						String[] splitted = line.split(":");
						if(splitted.length==2) {
							String header = splitted[1].trim();
							header = header.replaceAll("'", "");
							header = header.replaceAll(">", "");
							header = header.replaceAll(" ", "");
							String[] splittedHeader = header.split("\\|");
							StringBuffer headerLine = new StringBuffer();
							for(String column:splittedHeader) {
								headerLine.append(column);
								headerLine.append("\t");
							}
							List<String> head = new LinkedList<String>();
							head.add(headerLine.toString().trim());
							result.put(-1, head);
						}
					}else {
						continue;
					}
				}else {
					String[] splitted = line.split("\t");
					Integer pos = Integer.parseInt(splitted[1]);
					for(String s: splitted) {
						if(s.startsWith("ANN=")) {
							s = s.substring(4);
							String[] splittedEntries = s.split(",");
							for(String entry: splittedEntries) {
								String[] splittedData = entry.split("\\|");
								StringBuffer dataLine = new StringBuffer();
								for(String data: splittedData) {
									dataLine.append(data);
									dataLine.append("\t");
								}
								List<String> dat = new LinkedList<String>();
								if(result.containsKey(pos)) {
									dat = result.get(pos);
								}
								dat.add(dataLine.toString().trim());
								result.put(pos, dat);
							}
							break;
						}
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return result;
	}

	public String toStringWithSNPEffInfos(String snpEffOutputFile) {
		StringBuffer result = new StringBuffer();
		Map<Integer, List<String>> snpEffOutput = parseSNPEffOutput(snpEffOutputFile);
		result.append("Position");
		result.append("\t");
		result.append("Ref");
		for(String sample: this.sampleNames){
			result.append("\t");
			result.append(sample);
		}
		result.append("\t");
		result.append(snpEffOutput.get(-1).get(0));
		result.append("\n");
		for(Integer pos: this.snpPos){
			if(snpEffOutput.containsKey(pos)) {
				List<String> snpEffDat = snpEffOutput.get(pos);
				for(String dat: snpEffDat) {
					result.append(pos);
					result.append("\t");
					result.append(getReferenceSnp(pos));
					for(String sample: this.sampleNames){
						result.append("\t");
						result.append(getSnp(pos,sample));
					}
					result.append("\t");
					result.append(dat);
					result.append("\n");
				}
			}else {
				result.append(pos);
				result.append("\t");
				result.append(getReferenceSnp(pos));
				for(String sample: this.sampleNames){
					result.append("\t");
					result.append(getSnp(pos,sample));
				}
				result.append("\n");
			}
		}
		return result.toString();
	}
	
	public String toStringHeterozygousOnly() {
		StringBuffer result = new StringBuffer();
		result.append("Position");
		result.append("\t");
		result.append("Ref");
		for(String sample: this.sampleNames){
			result.append("\t");
			result.append(sample);
		}
		result.append("\n");
		Set<String> possibleSNPs = new HashSet<String>();
		possibleSNPs.add("A");
		possibleSNPs.add("C");
		possibleSNPs.add("G");
		possibleSNPs.add("T");
		possibleSNPs.add("N");
		possibleSNPs.add(".");
		for(Integer pos: this.snpPos){
			StringBuffer currLine = new StringBuffer();
			Boolean hasHetCall = false;
			currLine.append(pos);
			currLine.append("\t");
			String refSNP = getReferenceSnp(pos);
			currLine.append(refSNP);
			for(String sample: this.sampleNames){
				currLine.append("\t");
				String currSNP = getSnp(pos,sample);
				currLine.append(currSNP);
				if(!possibleSNPs.contains(currSNP)) {
					hasHetCall = true;
				}
			}
			currLine.append("\n");
			if(hasHetCall) {
				result.append(currLine.toString());
			}
		}
		return result.toString();
	}

	public String toString(){
		StringBuffer result = new StringBuffer();
		result.append("Position");
		result.append("\t");
		result.append("Ref");
		for(String sample: this.sampleNames){
			result.append("\t");
			result.append(sample);
		}
		result.append("\n");
		for(Integer pos: this.snpPos){
			result.append(pos);
			result.append("\t");
			String refSNP = getReferenceSnp(pos);
			result.append(refSNP);
			for(String sample: this.sampleNames){
				result.append("\t");
				String currSNP = getSnp(pos,sample);
				result.append(currSNP);
			}
			result.append("\n");
		}
		return result.toString();
	}

	public String toSnpTableForSnpEffString(String referenceName) {

	String refName = referenceName.replace(" ", "_");
		StringBuffer result = new StringBuffer();
		for(Integer pos: this.snpPos) {
			String refSNP = getReferenceSnp(pos); 
			Set<String> snps = new TreeSet<String>();
			for(String sample: this.getSampleNames()) {
				String snp = getSnp(pos,sample);
				if(!"N".equals(snp) && !".".equals(snp)) {
					if(".".equals(snp)){
						snp = getSnp(pos, "Ref");
					}
					snps.add(snp);
				}
			}
			for(String snp: snps) {
				result.append(refName);
				result.append("\t");
				result.append(pos);
				result.append("\t");
				result.append(refSNP);
				result.append("\t");
				result.append(snp);
				result.append("\n");
			}
		}
		return result.toString();
	}

	public String toSnpVcfForSnpEffString(String referenceName) {
		StringBuffer result = new StringBuffer();
		Set<String> alreadyDone = new HashSet<String>();
		for(Integer pos: this.snpPos) {
			String refSNV = getReferenceSnp(pos); 
			Set<String> snps = new TreeSet<String>();
			for(String sample: this.getSampleNames()) {
				StringBuffer intermediateResult = new StringBuffer(); 
				if(alreadyDone.contains(sample+pos)) {
					continue;
				}
				String snp = getSnp(pos,sample);
				if(!"N".equals(snp) && !".".equals(snp)) {
					if(".".equals(snp)){
						snp = getSnp(pos, "Ref");
					}
					if("-".equals(snp)) {
						snp = getReferencePos(pos-1);
						refSNV = getReferencePos(pos-1);
						Integer tmpPos = pos;
						while(this.snpPos.contains(tmpPos) && "-".equals(getSnp(tmpPos, sample))) {
							refSNV += getReferencePos(tmpPos);
							tmpPos++;
							alreadyDone.add(sample+tmpPos);
						}
						pos--;
					}
//					snps.add(snp);
					intermediateResult.append(referenceName);
					intermediateResult.append("\t");
					intermediateResult.append(pos);
					intermediateResult.append("\t");
					intermediateResult.append(".");
					intermediateResult.append("\t");
					intermediateResult.append(refSNV);
					intermediateResult.append("\t");
					intermediateResult.append(snp);
					intermediateResult.append("\n");
					snps.add(intermediateResult.toString());
				}
			}
			for(String snp: snps) {
				result.append(snp);
//				result.append(referenceName);
//				result.append("\t");
//				result.append(pos);
//				result.append("\t");
//				result.append(".");
//				result.append("\t");
//				result.append(refSNV);
//				result.append("\t");
//				result.append(snp);
//				result.append("\n");
			}
		}
		return result.toString();
	}

//	public String toAlignment(){
//		return toAlignment(new HashSet<Integer>(), new TreeMap<Integer, Map<String, Integer>>());
//	}
	
	public String toAlignment(Map<Integer, Map<String, Integer>> dels) {
		return toAlignment(new HashSet<Integer>(), dels);
	}
	
//	public String toAlignment(Set<Integer> toExclude) {
//		return toAlignment(toExclude, new TreeMap<Integer, Map<String, Integer>>());
//	}

	/**
	 * @param toExclude
	 * @return
	 */
	public String toAlignment(Set<Integer> toExclude, Map<Integer, Map<String, Integer>> dels) {
		StringBuffer result = new StringBuffer();
		for(String sample:this.sampleNames){
			result.append(">");
			result.append(sample.replace(" ", "_"));
			result.append("\n");
			int snpNum=0;
			for(Integer pos: this.snpPos){
				if(toExclude.contains(pos)){
					continue;
				}
				String snp = getSnp(pos, sample);
				if(".".equals(snp)){
					result.append(getSnp(pos, "Ref"));
				}else{
					result.append(snp);
				}
				snpNum++;
				// handle the insertions and deletions
				int max = 0;
				if(dels.containsKey(pos)){
					int min=0;
					for(String sampleDel: dels.get(pos).keySet()) {
						Integer numDeletions = dels.get(pos).get(sampleDel);
						if(sampleDel.equals(sample)) {
							min = numDeletions;
						}
						if(numDeletions > max) {
							max = numDeletions;
						}
					}
					while(min<max) {
						result.append("-");
						min++;
						snpNum++;
					}
				}
				if(snpNum%80==0){
					result.append("\n");
				}				
			}
			result.append("\n");
		}
		return result.toString();
	}

	public List<String> compareTo(SNVTable table){
		List<String> result = new LinkedList<String>();
		List<String> otherSampleNames = table.getSampleNames();
		Set<String> bothSampleNames = new TreeSet<String>();
		result.add("Difference\tPosition\tSample\tCurent\tOther");
		bothSampleNames.addAll(this.sampleNames);
		bothSampleNames.addAll(otherSampleNames);
//		for(String sample: this.sampleNames){
//			if(otherSampleNames.contains(sample)){
//				bothSampleNames.add(sample);
//			}
//		}
		List<Integer> othersnps = table.getSNVs();
		TreeSet<Integer> allPositions = new TreeSet<Integer>();
		allPositions.addAll(this.snpPos);
		allPositions.addAll(othersnps);
		for(Integer pos: allPositions){
			if(this.snpPos.contains(pos) && othersnps.contains(pos)){
				for(String sample: bothSampleNames){
					String thisSNP = "";
					String otherSNP = "";
					if(this.sampleNames.contains(sample)) {
						thisSNP = getSnp(pos, sample);
					}
					if(otherSampleNames.contains(sample)) {
						otherSNP = table.getSnp(pos, sample);
					}
					Set<String> snps = new HashSet<String>();
					snps.add(thisSNP);
					snps.add(otherSNP);
					Boolean meaningfullDifference = false;
					for(String snp: snps) {
						if(!"".equals(snp) && !"N".equals(snp)) { //TODO && !".".equals(snp)
							meaningfullDifference = true;
						}
					}
					if(!thisSNP.equals(otherSNP) && meaningfullDifference){
						result.add("=\t"+pos+"\t"+sample+"\t"+thisSNP+"\t"+otherSNP);
					}
				}
			}else if(this.snpPos.contains(pos)){
				for(String sample: this.sampleNames){
					String snp = getSnp(pos, sample);
					if(!"N".equals(snp) && !".".equals(snp)){
						result.add("+\t"+pos+"\t"+sample+"\t"+snp+"\t");
					}
				}
			}else{
				for(String sample: table.getSampleNames()){
					String snp = table.getSnp(pos, sample);
					if(!"N".equals(snp) && !".".equals(snp)){
						result.add("-\t"+pos+"\t"+sample+"\t\t"+snp);
					}
				}
			}
		}
		return result;
	}

}
