package datastructure;
import java.util.Set;

//import datastructure.protein.CodonsAminoacids;

public class Gene implements Comparable<Gene>
{
	String id;
	String name;
	String panId;
	String blockId;
	int start;
	int end;
	char strand;

	String genomeID;
	String description;
	String type;
	String source;
	String origin; //which strain
	String idAsParent;
	String note;
	String consensusDescription;
	String consensusName;

	boolean mappabilityIssue = false;

	String tigr;

	Gene nextGene;
	Gene previousGene;

	//	List<TSS> utrTsss;

	public String corrected;
	/**
	 * Normal Constructor
	 * @param source
	 * @param origin
	 * @param id
	 * @param type
	 * @param start
	 * @param end
	 * @param strand
	 * @param description
	 */
	public Gene(String source, String origin, String id,String name, String type, int start, int end, char strand, String genomeID, String description) 
	{
		super();
		this.source = source;
		this.origin = origin;
		this.id = id;
		this.name= name;
		this.panId = "";
		this.type = type;
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.genomeID = genomeID;
		this.description=description;

		//		utrTsss =new LinkedList<TSS>();
	}
	/**
	 *  Constructor with superID
	 * @param source
	 * @param origin
	 * @param id
	 * @param superId
	 * @param type
	 * @param start
	 * @param end
	 * @param strand
	 * @param description
	 */
	public Gene(String source, String origin, String id, String name, String panId, String type, int start, int end, char strand, String genomeID,String description) 
	{
		super();
		this.source = source;
		this.origin = origin;
		this.id = id;
		this.name= name;
		this.panId = panId;
		this.type = type;
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.genomeID = genomeID;
		this.description=description;

		//		utrTsss =new LinkedList<TSS>();
	}

	/**
	 *  Constructor with superID and previousGene, nextGene
	 * @param source
	 * @param origin
	 * @param id
	 * @param superId
	 * @param type
	 * @param start
	 * @param end
	 * @param strand
	 * @param description
	 */
	public Gene(String source, String origin, String id, String name, String panId, String type, int start, int end, char strand,String genomeID, String description, Gene previousGene, Gene nextGene) 
	{
		super();
		this.source = source;
		this.origin = origin;
		this.id = id;
		this.name= name;
		this.panId = panId;
		this.type = type;
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.description=description;
		this.genomeID = genomeID;
		this.nextGene = nextGene;
		this.previousGene = previousGene;
		//		utrTsss =new LinkedList<TSS>();
	}


	/**
	 * Constructor with idAsParent
	 * @param source
	 * @param origin
	 * @param id
	 * @param type
	 * @param start
	 * @param end
	 * @param strand
	 * @param description
	 * @param idAsParent
	 */
	public Gene(String source, String origin,String id,  String name,  String type, int start, int end, char strand, String genomeID,String description, String idAsParent) 
	{
		//this(source, origin, id, name , type, start, end, strand, genomeID, description);
		this.source = source;
		this.origin = origin;
		this.id = id;
		this.name= name;
		this.type = type;
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.description=description;
		this.genomeID = genomeID;
		this.panId = "";
		this.idAsParent = idAsParent;

	}

	/**
	 * Constructor with consensusDescription and note
	 * @param source
	 * @param origin
	 * @param id
	 * @param type
	 * @param start
	 * @param end
	 * @param strand
	 * @param description
	 * @param consensusDescription
	 * @param note
	 */
	//	public Gene(String source, String origin, String id, String name, String type, int start, int end, char strand, String description , String note) 
	//	{
	//		super();
	//		this.source = source;
	//		this.origin = origin;
	//		this.id = id;
	//		this.name= name;
	//		this.panId = "";
	//		this.type = type;
	//		this.start = start;
	//		this.end = end;
	//		this.strand = strand;
	//		this.description=description;
	//		this.note = note;
	////		utrTsss =new LinkedList<TSS>();
	//	}

	/**
	 * Constructor with consensusDescription, note, nextGene and previousGene
	 * @param source
	 * @param origin
	 * @param id
	 * @param type
	 * @param start
	 * @param end
	 * @param strand
	 * @param description
	 * @param consensusDescription
	 * @param note 
	 * @param nextGene
	 * @param previousGene
	 */
	public Gene(String source, String origin, String id, String name, String type, int start, int end, char strand, String genomeID,String description, String note, Gene nextGene, Gene previousGene) 
	{
		super();
		this.source = source;
		this.origin = origin;
		this.id = id;
		this.name= name;
		this.panId = "";
		this.type = type;
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.description=description;
		this.note = note;
		this.nextGene = nextGene;
		this.previousGene = previousGene;
		//		utrTsss =new LinkedList<TSS>();
	}

	/**
	 * Normal Constructor with nextGene and previous Gene
	 * @param source
	 * @param origin
	 * @param id
	 * @param type
	 * @param start
	 * @param end
	 * @param strand
	 * @param description
	 * @param nextGene
	 * @param previousGene
	 */
	public Gene(String source, String origin, String id,String name, String type, int start, int end, char strand, String genomeID,String description, Gene nextGene, Gene previousGene) 
	{
		super();
		this.source = source;
		this.origin = origin;
		this.id = id;
		this.name= name;
		this.panId = "";
		this.type = type;
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.description=description;
		this.nextGene = nextGene;
		this.previousGene = previousGene;
		//		utrTsss =new LinkedList<TSS>();
	}



	public String getGenomeID() {

		return this.genomeID;
	}


	public String getBlockId() {
		return blockId;
	}
	public void setBlockId(String blockId) {
		this.blockId = blockId;
	}
	public String getName(){
		return this.name;
	}	

	public Gene getNextGene() {
		return nextGene;
	}
	public void setNextGene(Gene nextGene) {
		this.nextGene = nextGene;
	}
	public Gene getPreviousGene() {
		return previousGene;
	}
	public void setPreviousGene(Gene previousGene) {
		this.previousGene = previousGene;
	}
	public String getId() {
		return id;
	}

	public String getConsensusDescription(){
		return this.consensusDescription;
	}
	public String getConsensusName(){
		return this.consensusName;
	}

	public String getIdAsParent() {
		return idAsParent;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public char getStrand() {
		return strand;
	}

	public int getLength()
	{
		return end-start+1;
	}


	public void setConsensusDescription(String description){
		this.consensusDescription = description;
	}

	public void setConsensusName(String name){
		this.consensusName = name;
	}



	public String getDescription()
	{
		return(description);
	}


	public void setName(String name){
		this.name = name;
	}



	public String getPanId(){
		return this.panId;
	}

	public void setPanId(String id){
		this.panId = id;
	}


	public void markMappabilityIssue(){

		this.mappabilityIssue = true;

	}

	//	public void addUTRtss(TSS tss)
	//	{
	//		utrTsss.add(tss);
	//	}
	//
	//	public void classifyUTRtsss()
	//	{
	//		if(this.utrTsss.size()==0)
	//			return;
	//
	//		TSS best = this.utrTsss.get(0);
	//
	//		for(TSS tss : utrTsss)
	//			if(tss.getHeight()>best.getHeight())
	//				best=tss;
	//
	//		//		for(TSS tss : utrTsss)
	//		//			if(tss.distanceTo(this)<best.distanceTo(this))
	//		//				best=tss;
	//
	//		best.setPrimary(this);
	//	}

	public String getType() {
		return type;
	}

	public String getSource() {
		return source;
	}

	public void setSource(String source) {
		this.source = source;
	}

	public String toGFFString()
	{
		//XXX type is overwritten to gene (RNA tag gets lost)
		String res = source+"\t"+origin+"\t"+"gene"+"\t"+start+"\t"+end+"\t"+"."+"\t"+strand+"\t"+"."+"\t"+"locus_tag="+id;

		if(description.length()!=0)
			res = res+";product="+description;
		if(note!=null)
			res = res+";Note="+note;
		if(name!=null && name.length()!=0)
			res = res+ ";gene="+name;
		if(this.panId.length()!=0)
			res = res+";superId="+panId;
		return res;
	}

	public String getOrigin() {
		return origin;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public String getNote() {
		return note;
	}

	public void setNote(String note) {
		this.note = note;
	}

	public void setID(String id){
		this.id = id;
	}

	public void setStrand(char strand){
		this.strand = strand;
	}

	public void setType(String type){
		this.type = type;
	}

	public void setOrigin(String origin){
		this.origin =origin;
	}

	public boolean getMappabilityIssue(){
		return this.mappabilityIssue;
	}

	public void setTigr(String tigr){
		this.tigr =tigr;
	}

	public String getTigr(){
		if(this.tigr != null){
			return this.tigr;
		} else {
			return "";
		}
	}

	public int[] getOverlapStartEnd(Gene gene)
	{
		int oStart = Math.max(this.getStart(),gene.getStart());
		int oEnd   = Math.min(this.getEnd()  , gene.getEnd())  ;

		//no overlap
		if(oEnd<oStart)
		{
			return null; //getNumOverlappingBases and maybe other functions depend on this!
		}

		int[] res = {oStart,oEnd};

		return res;
	}

	public int getNumOverlappingBases(Gene gene)
	{
		int[] ov = getOverlapStartEnd(gene);

		if(ov == null)
			return 0;

		return ov[1]-ov[0]+1;
	}

	public String get_sRNA_asRNA_labelIfContainedInSet(Set<String> sRNAnameList, Set<String> asRNAnameList)
	{
		if(sRNAnameList.contains(id) && asRNAnameList.contains(id))
			return "sRNA/asRNA";
		if(sRNAnameList.contains(id))
			return "sRNA";
		if(asRNAnameList.contains(id))
			return "asRNA";
		else
			return "";
	}

	/**
	 * The equals method has been changed and only compares the id!	
	 */
	//	public boolean equals(Object obj) {
	//		if (obj == this) {
	//			return true;
	//		}
	//		SuperGenomifiedGene sg = new SuperGenomifiedGene(source, origin, id, name, type, start, end, strand, description,  start, end, strand );
	//		if (obj == null || (obj.getClass() != this.getClass()&& obj.getClass() != sg.getClass())) {
	//			return false;
	//		}
	//		if(obj.getClass() == sg.getClass()){
	//			SuperGenomifiedGene compare = (SuperGenomifiedGene) obj;
	//			return this.id.equals(compare.id);
	//		}
	//		Gene compare = (Gene) obj;
	//		return this.id.equals(compare.id);
	//	}
	public boolean equals(Object obj) {
		if (obj == this) {
			return true;
		}
		Gene compare = (Gene) obj;
		return this.id.equals(compare.id);
	}

	/**
	 * If this starts before o than -1 is returned. If this and o starts at the same position the end is used to compare. if they are also equal 0 is returned.
	 */
	public int compareTo(Gene o) {
		if (this.getStart()== o.getStart()) {
			if (this.getEnd()< o.getEnd()) {
				return -1;
			}else if (this.getEnd()> o.getEnd()) {
				return 1;
			}else{
								
				return this.id.compareTo(o.getId());
			}
		}else if (this.getStart()< o.getStart()) {
			return -1;
		}else{//(this.getStart()> o.getStart()) 
			return 1;
		}
	}


	public String toString(){

		String s = "";

		s = id + "\tstrand\t"+ origin + "\tName\t" + name + "\tstart\t" + start + "\tend\t" + end + "\tstrand\t" + strand + "\n";

		return s;
	}


	public String toStringForSnapshot() {

		String s = "";

		s += "origin:" + origin + "|";
		s += "id:" + id + "|";
		s += "name:" + name + "|";
		s += "type:" + type + "|";
		s += "start:" + start + "|";
		s += "end:" + end + "|";
		s += "strand:" + strand + "|";
		s += "genomeID:" + genomeID + "|";
		s += "description:" + description + "|";


		return s;

	}


//	public static String translateToProteinSeq(String seq){
//		String proteinSequence = "";
//		CodonsAminoacids mapping = new CodonsAminoacids();
//		
//		
//		if(seq.length() % 3 !=0){
//
//			System.out.println("The Sequence can not be translated (sequencelength modulo 3 != 0)");
//			return seq;
//
//		}
//
//
//		for(int i =0; i< seq.length();i=i+3){
//
//			CodonsAminoacids.Codon codon = mapping.new Codon(seq.charAt(i), seq.charAt(i+1),seq.charAt(i+2));			
//			proteinSequence = proteinSequence+ mapping.translate(codon).oneLetter;
//			
//		}
//
//		return proteinSequence;
//	}


}