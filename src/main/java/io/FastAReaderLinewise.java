package io;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import datastructure.FastAEntry;

public class FastAReaderLinewise {
	private BufferedReader br;
	private String currLine = "";
	
	public FastAReaderLinewise(String file){
		try {
			this.br = new BufferedReader(new FileReader(file));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	private void readOneLine(){
		currLine = null;
		try {
			currLine = this.br.readLine();
			if(currLine != null){
				currLine = currLine.trim();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * 
	 * @return the next fastA entry from the file
	 */
	public FastAEntry getOneFastAEntry(){
		if(currLine == null){
			return null;
		}
		while(currLine != null && !currLine.startsWith(">")){
			readOneLine();
		}
		if(currLine == null){
			return null;
		}
		String header = currLine.substring(1).trim();
		StringBuffer sequence = new StringBuffer();
		readOneLine();
		while(currLine != null && !currLine.startsWith(">")){
			if("".equals(currLine)){
				readOneLine();
				continue;
			}
			if(currLine.startsWith(";")){
				readOneLine();
				continue;//TODO
			}
			sequence.append(currLine);
			readOneLine();
		}
		return new FastAEntry(header, sequence.toString().trim());
	}
	
	public List<String> getAllHeaders(){
		LinkedList<String> result = new LinkedList<String>();
		FastAEntry fae = getOneFastAEntry();
		while(fae != null){
			result.add(fae.getHeader());
			fae = getOneFastAEntry();
		}
		return result;
	}
}
