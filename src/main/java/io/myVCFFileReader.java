/**
 *
 */
package io;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import datastructure.FastAEntry;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import utilities.Pair;
import utilities.Statistics;

/**
 * @author Alexander Seitz
 *
 */
public class myVCFFileReader {

	public String filename;
	private String sampleName;
	private Integer genomeSize;
	private Map<Integer, String> genotypes;
	private Map<Integer, String> maybeGenotypes;
	private Map<Integer, Double> coverage;
	private Set<Integer> lowCoverageRegions;
	private String contigName = "";

//	private Set<Integer> insertions;
	private Set<Pair<Integer,Integer>> deletions;

	private Double minCovGoodCall = 5d;
	private Double minCovAdditionalCall = 3d;
	private double minFreq = 0.9;
	private Double genQual = 30d;
	private Set<Integer> varPos;
	private Set<Integer> heterozygousPos;
	private double mean;
	private double stdev;
	private int currDeletions = 0;

	private Boolean callHeterozygous = false;
	private Double minHet = 0.25;
	private Double maxHet = 0.75;
//	private Long totalCalls = 0l;
//	private Long refCalls = 0l;
//	private Long noCalls = 0l;


	public myVCFFileReader(String file, String sampleName, Integer genomeSize, Double minCovGood, Double minCovAdd, Double minFreq, Double genQual, Boolean callHeterozygous, Double minHet, Double maxHet){
		this.filename = file;
		this.sampleName = sampleName;
		this.genomeSize = genomeSize;
		this.minCovGoodCall = minCovGood;
		this.minCovAdditionalCall = minCovAdd;
		this.minFreq = minFreq;
		this.genQual = genQual;
		this.genotypes = new HashMap<Integer, String>();
		this.maybeGenotypes = new HashMap<Integer, String>();
		this.coverage = new HashMap<Integer, Double>();
		this.lowCoverageRegions = new HashSet<Integer>();
		this.varPos = new TreeSet<Integer>();
		this.heterozygousPos = new TreeSet<Integer>();
//		this.insertions = new TreeSet<Integer>();
		this.deletions = new TreeSet<Pair<Integer,Integer>>();
		this.mean = 0d;
		this.stdev = 0d;
		this.callHeterozygous = callHeterozygous;
		this.minHet = minHet;
		this.maxHet = maxHet;
		//		run();
	}


	public myVCFFileReader(String referenceAlignmentName, String refSequence, Integer genomeSize, Double minCovGood) {
		this.sampleName = referenceAlignmentName;
		this.genomeSize = genomeSize;
		this.genotypes = new HashMap<Integer, String>();
		this.maybeGenotypes = new HashMap<Integer, String>();
		this.coverage = new HashMap<Integer, Double>();
		this.lowCoverageRegions = new HashSet<Integer>();
		this.varPos = new TreeSet<Integer>();
		this.heterozygousPos = new TreeSet<Integer>();
//		this.insertions = new TreeSet<Integer>();
		this.deletions = new TreeSet<Pair<Integer,Integer>>();
		for(int i=1; i<=genomeSize; i++) {
			this.genotypes.put(i, ".");
			this.coverage.put(i, minCovGood);
		}
	}


	/**
	 *
	 */
	public void parseVCFFile() {
		@SuppressWarnings("resource")
		VCFFileReader vcfFileReader = new VCFFileReader(new File(this.filename), false);
		Iterator<VariantContext> iter = vcfFileReader.iterator();
		Integer lastPos=0;
		String previousGenotype = "N";
		Boolean haplotypeCaller = false;
		if(vcfFileReader.getFileHeader().toString().contains("HaplotypeCaller")) {
			haplotypeCaller = true;
		}
		while(iter.hasNext()){
			VariantContext variantContext = iter.next();
			if(this.contigName.length() == 0) {
				this.contigName = variantContext.getContig();
			}
//			this.totalCalls++;
			int pos = variantContext.getStart();
			if(!(pos-lastPos==1)){
				while(this.currDeletions>0 && !(pos-lastPos==1)) {
					this.genotypes.put(lastPos+1,"-");
					lastPos++;
					this.currDeletions--;
					this.varPos.add(lastPos);
				}
				addGenos(lastPos,pos,previousGenotype);
			}
			if(this.currDeletions>0) {
				System.out.println("now: "+this.sampleName+" "+pos);
				this.genotypes.put(pos, "-");
				this.currDeletions--;
				lastPos = pos;
				this.varPos.add(lastPos);
				continue;
			}
			String currentGenotype = "N";
			Double def = 0.5;
			Double coverage = variantContext.getAttributeAsDouble("DP", def);
			if(coverage.equals(def)) {
				coverage = getRealCoverage(variantContext);
			}
			this.coverage.put(pos, coverage);
			if(coverage < this.minCovAdditionalCall){
				currentGenotype="N";
			}else if(coverage >= this.minCovGoodCall){
				currentGenotype = calculateGenotypeCall(variantContext, pos, this.minCovGoodCall);
			}else{ // maybe calls
				currentGenotype="N";
				String maybeGenotype = calculateGenotypeCall(variantContext, pos, this.minCovAdditionalCall);
				if(".".equals(maybeGenotype)){
					currentGenotype=".";
				}else if(!"N".equals(maybeGenotype)){
					this.maybeGenotypes.put(pos, maybeGenotype);
				}
			}
			// if we used the haplotyper caller and not the genotype caller
			// if a position is missing, use the last one and not an "N"
			if(haplotypeCaller) {
				previousGenotype = currentGenotype;
			}
			this.genotypes.put(pos, currentGenotype);
			lastPos=pos;
		}
		addGenos(lastPos, this.genomeSize+1, previousGenotype);
		Statistics stats = new Statistics(this.coverage);
		this.mean = stats.getMean();
		this.stdev = stats.getStdDev();
		// search for low coverage regions
		for(Integer pos: this.genotypes.keySet()){
			if(this.coverage.containsKey(pos) && (this.coverage.get(pos) < this.minCovAdditionalCall || this.coverage.get(pos) <= this.mean - this.stdev)){
				this.lowCoverageRegions.add(pos);
			}
		}

		// search for high coverage regions and remove them
//		if(stdev>mean){
//			System.out.println(new File(new File(this.filename).getParent()).getName()+"\t"+mean+"\t"+stdev);
//			for(Integer pos: this.genotypes.keySet()){
//				if(this.coverage.get(pos) >= ( this.mean+this.stdev)){
//					if(!".".equals(this.genotypes.get(pos))){
//						this.genotypes.put(pos, "N");
//					}
//				}
//			}
//		}
	}


	/**
	 * @param variantContext
	 * @return the coverage of the variantContext
	 */
	private Double getRealCoverage(VariantContext variantContext) {
		Double result = 0.0;
		for(Genotype gc: variantContext.getGenotypes()) {
			if(gc.hasDP()) {
				result += (double) gc.getDP();
			}
		}
		return result;
	}


	/**
	 * @param variantContext
	 * @param pos
	 * @return
	 */
	private String calculateGenotypeCall(VariantContext variantContext, Integer pos, Double minCov) {
		String currentGenotype = "N";
		// we only work with one genotype
		if(variantContext.getGenotypes().size() >1){
			currentGenotype="N";
		}else{
			Genotype g = variantContext.getGenotype(0);
			Double phredQual = variantContext.getCommonInfo().getPhredScaledQual();
			if(phredQual == null || phredQual < 0) {
				phredQual=100.0; //TODO
//				phredQual = 0.0;
//				int[] qualities = g.getPL();
//				if(qualities != null) {
//					for(int qual: qualities) {
//						phredQual = Math.max(phredQual, qual);
//					}
//				}
			}
			if(phredQual < this.genQual){
				currentGenotype="N";
			}else if(g.isHomRef()){ // homozygot reference call
				currentGenotype=".";
				// homozygot variant call
			}else if(g.isHomVar()){
				int[] ad = g.getAD();
				int sum = 0;
				int max=0;
				int maxid=0;
				int currid=0;
				for(int v: ad){
					if(v>max){
						max=v;
						maxid=currid;
					}
					currid++;
					sum += v;
				}
				if(sum<g.getDP()){
					if(sum<minCov){
						currentGenotype="N";
					}else if((double)max/(double)g.getDP()>=this.minFreq){
						if(maxid==0){
							currentGenotype=".";
						}else if(maxid==1){
							String[] splitted = g.getGenotypeString().split("/");
							currentGenotype=splitted[0];
							this.varPos.add(pos);
						}else{
							//TODO check everything
							System.out.println("Problem here");
							System.out.println(pos);
							System.out.println(g);
							currentGenotype="N";
						}
					}else{
						currentGenotype="N";
					}
				}else if(ad[0]>1 && ((double) ad[1]/(double) g.getDP() < this.minFreq)){
					currentGenotype="N";
				}else if(ad[0]==0 && ad[1] < g.getDP()){
					currentGenotype="N";
				}else{
					if(minCov < this.minCovGoodCall){
						currentGenotype="N";
					}else{
						String[] splitted = g.getGenotypeString().split("/");
						currentGenotype=splitted[0];
						for(Allele allele: variantContext.getAlleles()) {
							if(allele.isReference()) {
								if(allele.getBaseString().length()>currentGenotype.length()) {
									this.currDeletions = allele.getBaseString().length()-currentGenotype.length();
									System.out.println("longer snp here?\t"+pos+"\t"+(allele.getBaseString().length()-currentGenotype.length()));
								}else if(allele.getBaseString().length()<currentGenotype.length()) {
//										Integer refLength = allele.getBaseString().length();
										Integer lengthDiff = (currentGenotype.length()-allele.getBaseString().length());
										Pair<Integer, Integer> p = new Pair<Integer, Integer>(pos, lengthDiff);
										this.deletions.add(p);
								}
							}
						}
						this.varPos.add(pos);
					}
				}
				// heterozygot call
			}else if(g.isHet()){
				int[] ad = g.getAD();
				String[] splitted = g.getGenotypeString().split("/");
				if(ad.length==2){
					int ref = ad[0];
					int var = ad[1];
					if(ref==1){
						currentGenotype=splitted[1];
						this.varPos.add(pos);
					}else if(var==1){
						currentGenotype=".";
					}else if((double)ref/(double)g.getDP()>=this.minFreq){
						currentGenotype=".";
					}else if((double)var/(double)g.getDP()>=this.minFreq){
						currentGenotype=splitted[1];
						this.varPos.add(pos);
					}else if(this.callHeterozygous && (g.getDP() > 2*this.minCovGoodCall) && (((double)ref/(double)g.getDP()<=this.maxHet && (double)ref/(double)g.getDP()>=this.minHet)
						|| ((double)var/(double)g.getDP()<=this.maxHet && (double)var/(double)g.getDP()>=this.minHet))){
							String genotype1 = splitted[0];
							String genotype2 = splitted[1];
							currentGenotype=getCorrectUPACCode(genotype1, genotype2);
							this.varPos.add(pos);
							this.heterozygousPos.add(pos);
					}else{
						currentGenotype="N";
					}
				}else{ // heterozygous, more than one alternative allele
					currentGenotype="N";
					for(int i=0; i<ad.length; i++){
						double freq = (double)ad[i]/(double)g.getDP();
						if(freq >= this.minFreq){
							if(i==0){
								currentGenotype = ".";
							}else{
								if(splitted.length>=i-1){
									currentGenotype=splitted[i-1];
								}
							}
							break;
						}
					}
				}
			}else{ // neither homozygot var/ref or heterozygot
				currentGenotype="N";
			}
		}
//		if(".".equals(currentGenotype)) {
//			this.refCalls++;
//		}else if("N".equals(currentGenotype)) {
//			this.noCalls++;
//		}
		return currentGenotype;
	}


	/**
	 * @param genotype1
	 * @param genotype2
	 * @return
	 */
	private String getCorrectUPACCode(String genotype1, String genotype2) {
		String genotype = "N";
		String combined = genotype1+genotype2;
		if("AG".equals(combined) || "GA".equals(combined)) {
			genotype="R";
		}else if("CT".equals(combined) || "TC".equals(combined)) {
			genotype="Y";
		}else if("GC".equals(combined) || "CG".equals(combined)) {
			genotype="S";
		}else if("AT".equals(combined) || "TA".equals(combined)) {
			genotype="W";
		}else if("GT".equals(combined) || "TG".equals(combined)) {
			genotype="K";
		}else if("AC".equals(combined) || "CA".equals(combined)) {
			genotype="M";
		}
		return genotype;
	}


	/**
	 * @param lastPos
	 * @param pos
	 */
	private void addGenos(Integer lastPos, int pos, String genotype) {
		for(int i=lastPos+1; i<pos; i++){
			this.coverage.put(i, 0d);
			this.genotypes.put(i, genotype);
//			this.totalCalls++;
		}
	}

	public Set<Integer> getVariablePositions(){
		return this.varPos;
	}


	/**
	 * @param position
	 * @return
	 */
	public String get(Integer position) {
		if(this.genotypes.containsKey(position)) {
			return this.genotypes.get(position);
		}else {
			return "N";
		}
	}

	/**
	 *
	 * @return
	 */
	public Map<Integer, String> getMaybeGenotypes(){
		return this.maybeGenotypes;
	}


	/**
	 * @return the mean
	 */
	public double getMean() {
		return mean;
	}


	/**
	 * @return the stdev
	 */
	public double getStdev() {
		return stdev;
	}


	/**
	 * @return the genomeSize
	 */
	public Integer getGenomeSize() {
		return genomeSize;
	}

	public String getGeneSequence(FastAEntry reference, Map<Integer, Map<String, Integer>> dels, Integer start, Integer end) {
		StringBuffer genomeSequence = new StringBuffer();
		for(int i=1; i<=this.genomeSize; i++){
			if(i<start || i > end) {
				continue;
			}
			int max = 0;
			if(dels.containsKey(i)){
				int min=0;
				for(String sample: dels.get(i).keySet()) {
					Integer numDeletions = dels.get(i).get(sample);
					if(sample.equals(this.sampleName)) {
						min = numDeletions;
					}
					if(numDeletions > max) {
						max = numDeletions;
					}
				}
				while(min<max) {
					genomeSequence.append("-");
					min++;
				}
			}
			String genotype = this.genotypes.get(i);
			if(".".equals(genotype)){
				genotype = ""+reference.getSequence().charAt(i-1);
			}
			genomeSequence.append(genotype);
		}
		return genomeSequence.toString();
	}


	/**
	 * @param fastAEntry
	 * @return the genomeSequence
	 */
	public String getGenomeSequence(FastAEntry reference, Map<Integer, Map<String, Integer>> dels) {
		StringBuffer genomeSequence = new StringBuffer();
		for(int i=1; i<=this.genomeSize; i++){
			int max = 0;
			if(dels.containsKey(i)){
				int min=0;
				for(String sample: dels.get(i).keySet()) {
					Integer numDeletions = dels.get(i).get(sample);
					if(sample.equals(this.sampleName)) {
						min = numDeletions;
					}
					if(numDeletions > max) {
						max = numDeletions;
					}
				}
				while(min<max) {
					genomeSequence.append("-");
					min++;
				}
			}
			String genotype = this.genotypes.get(i);
			if(".".equals(genotype)){
				genotype = ""+reference.getSequence().charAt(i-1);
			}
			genomeSequence.append(genotype);
		}
		return genomeSequence.toString();
	}


	/**
	 * @return the lowCoverageRegions
	 */
	public Set<Integer> getLowCoverageRegions() {
		return lowCoverageRegions;
	}


	/**
	 * @param pos the genotype position of the SNP
	 * @param string the SNP
	 * @return true if this sample also has the respective SNP (also if only in lower coverage positions)
	 */
	public boolean canResolveCall(Integer pos, String otherCall) {
		// this sample has the respective SNP called
		if(this.genotypes.containsKey(pos)){
			String snp = this.genotypes.get(pos);
			if(snp.equals(otherCall)){
				return true;
			}
		}
		// this sample has the respective SNP maybe called
		if(this.maybeGenotypes.containsKey(pos)){
			String snp = this.maybeGenotypes.get(pos);
			if(snp.equals(otherCall)){
				return true;
			}
		}
		return false;
	}


	/**
	 * @param pos
	 * @param string
	 */
	public void resolveCall(Integer pos, String snp) {
		this.genotypes.put(pos, snp);
	}


	/**
	 * @param posToReAnalyze
	 * @return
	 */
	public Map<Integer, String> reAnalyzePositions(List<Integer> posToReAnalyze) {
		Map<Integer, String> vcfInfo = new TreeMap<Integer,String>();
		@SuppressWarnings("resource")
		VCFFileReader vcfFileReader = new VCFFileReader(new File(this.filename), false);
		Iterator<VariantContext> iter = vcfFileReader.iterator();
		while(iter.hasNext()){
			VariantContext variantContext = iter.next();
			int pos = variantContext.getStart();
			if(posToReAnalyze.contains(pos)){
				vcfInfo.put(pos, variantContext.toString());
			}
		}
		return vcfInfo;
	}

	public Set<Pair<Integer,Integer>> getDeletions(){
		return this.deletions;
	}

	public String getSampleName() {
		return this.sampleName;
	}


	/**
	 * @return
	 */
	public Long getNumSNPs() {
//		return (long) this.varPos.size();
		long result = 0l;
		for(Integer pos: this.genotypes.keySet()) {
			String snp = this.genotypes.get(pos);
			if(!"N".equals(snp) && !".".equals(snp)) {
				result++;
			}
		}
		return result;
	}


	/**
	 * @return
	 */
	public Long getNumHetSNPs() {
//		return (long) this.heterozygousPos.size();
		long result = 0l;
		for(Integer pos: this.genotypes.keySet()) {
			String snp = this.genotypes.get(pos);
			if(!"N".equals(snp) && !".".equals(snp) && !"A".equals(snp) && !"T".equals(snp) && !"G".equals(snp) && !"C".equals(snp)) {
				result++;
			}
		}
		return result;
	}

	public Set<Integer> getHeterozygousPos(){
		return this.heterozygousPos;
	}


	/**
	 * @return
	 */
	public Long getTotalCalls() {
		return (long) this.genotypes.size();
	}


	/**
	 * @return
	 */
	public Long getRefCalls() {
//		return this.refCalls ;
		long result = 0l;
		for(Integer pos: this.genotypes.keySet()) {
			String snp = this.genotypes.get(pos);
			if(".".equals(snp)) {
				result++;
			}
		}
		return result;
	}


	/**
	 * @return
	 */
	public Long getNoCalls() {
//		return this.noCalls;
		long result = 0l;
		for(Integer pos: this.genotypes.keySet()) {
			String snp = this.genotypes.get(pos);
			if("N".equals(snp)) {
				result++;
			}
		}
		return result;
	}


	/**
	 * @return
	 */
	public Double getFoldCoverage() {
		Double coverage = 0d;
		for(Integer pos: this.coverage.keySet()) {
			coverage += this.coverage.get(pos);
		}
		return coverage/this.coverage.size();
	}


	/**
	 * @return
	 */
	public Double getPercentCoverage() {
		Double coverage = 0d;
		for(Integer pos: this.coverage.keySet()) {
			if(this.coverage.get(pos) >= this.minCovGoodCall) {
				coverage++;
			}
		}
		return coverage/this.coverage.size();
	}

	public Double getCoverage(Integer pos) {
		return this.coverage.get(pos);
	}


	/**
	 * @return
	 */
	public List<Pair<Integer, Integer>> getUncovered() {
		List<Pair<Integer,Integer>> uncovered = new LinkedList<Pair<Integer,Integer>>();
		Integer start = Integer.MAX_VALUE;
		Integer end = Integer.MAX_VALUE;
		for(Integer pos: this.coverage.keySet()) {
			Double coverage = this.coverage.get(pos);
			if(coverage == 0) {
				if(start == Integer.MAX_VALUE) {
					start = pos;
					end = pos;
				}else if(end == pos-1) {
					end = pos;
				}
			}else if(!(start == Integer.MAX_VALUE)) {
				uncovered.add(new Pair<Integer,Integer>(start,end));
				start = Integer.MAX_VALUE;
				end = Integer.MAX_VALUE;
			}
		}
		if(!(start == Integer.MAX_VALUE)) {
			uncovered.add(new Pair<Integer,Integer>(start,end));
			start = Integer.MAX_VALUE;
			end = Integer.MAX_VALUE;
		}
		return uncovered;
	}

	public String getContigName() {
		return this.contigName;
	}

}
