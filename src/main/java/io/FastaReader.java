/**
 * 
 */
package io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import datastructure.FastAEntry;

/**
 * @author Alexander Seitz
 *
 */
public class FastaReader {
	
	private String infile;
	private List<String> headers;
	private Map<String, FastAEntry> sequences;
	
	public FastaReader(String inFile){
		this.infile = inFile;
		this.sequences = new HashMap<String, FastAEntry>();
		this.headers = new LinkedList<String>();
		try {
			@SuppressWarnings("resource")
			BufferedReader br = new BufferedReader(new FileReader(this.infile));
			String currLine = "";
			String currHeader = "";
			StringBuffer currSequence = new StringBuffer();
			while((currLine=br.readLine()) != null){
				if(currLine.startsWith(">")){
					if(currHeader.length()>0){
						this.sequences.put(currHeader, new FastAEntry(currHeader, currSequence.toString()));
						currSequence = new StringBuffer();
					}
					currHeader = currLine.substring(1);
					this.headers.add(currHeader);
				}else if(!currLine.startsWith(";")){
					currSequence.append(currLine.trim());
				}
			}
			this.sequences.put(currHeader, new FastAEntry(currHeader, currSequence.toString()));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public List<String> getHeaders(){
		return this.headers;
	}
	
	public FastAEntry getEntry(String header){
		if(this.sequences.containsKey(header)){
			return this.sequences.get(header);
		}else{
			return null;
		}
	}
	
	public FastAEntry getFirstEntry(){
		if(this.headers.size()>0){
			return this.sequences.get(this.headers.get(0));
		}else{
			return null;
		}
	}

}
