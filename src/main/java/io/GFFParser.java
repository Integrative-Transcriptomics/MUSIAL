/**
 * 
 */
package io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import datastructure.Gene;

/**
 * @author Alexander Seitz
 *
 */
public class GFFParser {
	private List<Gene> geneList;
	private String inputGffFile;
	private String genomeID = "0";
	
	public GFFParser(String filename) throws Exception{
		this.inputGffFile = filename;
		parseGFF();
	}
	
	public GFFParser(String filename, String genomeID) throws Exception{
		this.inputGffFile = filename;
		this.genomeID = genomeID;
		parseGFF();
	}

	private List<Gene> parseGFF() throws Exception {
		this.geneList = new LinkedList<Gene>();
		try {
			Map<String,String> locusTag2DescMap = new HashMap<String, String>();
			Map<String,String> parentId2DescMap = new HashMap<String, String>();
			Map<String,String> position2DescMap = new HashMap<String, String>();
			Set<String> locusTag2ExonSet = new HashSet<String>();
			Set<String> parentId2ExonSet = new HashSet<String>();
			Set<String> position2ExonSet = new HashSet<String>();
			BufferedReader r = new BufferedReader(new FileReader(this.inputGffFile));
//			String[] tmp_name = filename.split("/");
			//  String origin = tmp_name[tmp_name.length];

			String line;
			String id;
			String name;
			String type;
			int start;
			int end;
			char strand;
			String desc;
			String source;
			String idAsParent;


			Map<String,String> attributes = new HashMap<String,String>();
			String[] cells;
			String[] attcells;
			for(line=r.readLine();line!=null;line=r.readLine())
			{
				line=line.trim();
				if(line.length()==0)
					continue;
				if(line.charAt(0)=='#')
					continue;
				cells = line.split("[\\t]");
				//if(!(cells[2].equalsIgnoreCase("CDS") || cells[2].equalsIgnoreCase("exon")))
				//continue;
				//source = 0
				source = cells[0];
				//type = 2
				type = cells[2];
				//extract tRNA
				//              if(type.equalsIgnoreCase("trna")||type.equalsIgnoreCase("rRNA")){
				//                  rna = true;
				//              }
				// start = 3
				start = Integer.parseInt(cells[3]);
				// end = 4
				end = Integer.parseInt(cells[4]);
				// strand = 6
				strand = cells[6].charAt(0);
				attributes.clear();
				if(cells.length>=9)
					for(String s : cells[8].split(";"))
					{
						attcells = s.split("=");
						if(attcells.length==2)
							attributes.put(attcells[0].trim(), attcells[1].trim());
						else if(attcells.length==1)
							attributes.put(attcells[0].trim(), "");                   
					}
				if(attributes.containsKey("locus_tag"))
					id=attributes.get("locus_tag");
				else if(attributes.containsKey("ID"))
					id=attributes.get("ID");
				else
					id="locus"+start+"-"+end;
				if(attributes.containsKey("ID"))
					idAsParent = attributes.get("ID");
				else
					idAsParent = null;
				if(attributes.containsKey("product"))
					desc=attributes.get("product");
				else
					desc="";
				if(attributes.containsKey("Name")){
					name = attributes.get("Name");
				}else{
					name = "";
				}
				//locus tag to desc mapping
				if(attributes.containsKey("locus_tag"))
				{
					if(attributes.containsKey("product"))
						locusTag2DescMap.put(attributes.get("locus_tag"), attributes.get("product"));
					if(attributes.containsKey("pseudo"))
						locusTag2DescMap.put(attributes.get("locus_tag"), "pseudo");
					if(type.equalsIgnoreCase("exon")){
						locusTag2ExonSet.add(attributes.get("locus_tag"));
					}
				}
				//parent id to desc mapping
				if(attributes.containsKey("Parent"))
				{
					if(attributes.containsKey("product"))
						parentId2DescMap.put(attributes.get("Parent"), attributes.get("product"));
					if(attributes.containsKey("pseudo"))
						parentId2DescMap.put(attributes.get("Parent"), "pseudo");
					if(type.equalsIgnoreCase("exon")){
						parentId2ExonSet.add(attributes.get("Parent"));
					}
				}
				//position to desc map
				if(attributes.containsKey("product"))
					position2DescMap.put(start+"_"+end+"_"+strand, attributes.get("product"));
				if(attributes.containsKey("pseudo"))
					position2DescMap.put(start+"_"+end+"_"+strand, "pseudo");
				if(type.equalsIgnoreCase("exon")){
					position2ExonSet.add(start+"_"+end+"_"+strand);
				}
				//create gene
				if(type.equalsIgnoreCase("gene") || type.equalsIgnoreCase("mRNA")){
					// set type to RNA to avoid lost of information
					this.geneList.add(new Gene(source,source,id,name,type,start,end,strand,this.genomeID,desc,idAsParent));
				}
			}
			r.close();
			//set descriptions
			for(Gene g: this.geneList){
				if(locusTag2ExonSet.contains(g.getId())){
					g.setType("exon");
				}else if(parentId2ExonSet.contains(g.getIdAsParent())){
					g.setType("exon");
				}else if(position2ExonSet.contains(g.getStart()+"_"+g.getEnd()+"_"+g.getStrand())){
					g.setType("exon");
				}
				if(locusTag2DescMap.containsKey(g.getId()))
					g.setDescription(locusTag2DescMap.get(g.getId()));
				else if(parentId2DescMap.containsKey(g.getIdAsParent()))
					g.setDescription(parentId2DescMap.get(g.getIdAsParent()));
				else if(position2DescMap.containsKey(g.getStart()+"_"+g.getEnd()+"_"+g.getStrand()))
					g.setDescription(position2DescMap.get(g.getStart()+"_"+g.getEnd()+"_"+g.getStrand()));
			}
		}catch(Exception e) {
			throw new Exception("A problem occured while parsing the GFF annotation:\n"+this.inputGffFile+"\n"+e.toString(),e);
		}
		if(this.geneList.size()==0){
			System.out.println(this.inputGffFile+" does not contain any genes!");
		}
		return this.geneList;
	}
	
	public Integer getNumberOfGenes(){
		return this.geneList.size();
	}
	
	public List<Gene> getGeneList(){
		return this.geneList;
	}
	
	public Set<Gene> getGenes(Set<String> geneNames){
		Set<Gene> result = new HashSet<Gene>();
		for(Gene gene: this.geneList){
			if(geneNames.contains(gene.getName())){
				result.add(gene);
			}
		}
		return result;
	}

}
