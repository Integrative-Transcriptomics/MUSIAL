/**
 *
 */
package tools;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import utilities.Inputs;
import utilities.Utilities;

/**
 * @author Alexander Seitz
 *
 */
public abstract class ATool {

	protected String eagerFolder;
	protected String outputFolder;
	protected List<File> eagerSamples;
	protected Map<String, String> vcfFilePaths = new TreeMap<String, String>();

	protected Options helpOptions = new Options();
	protected Options options = new Options();
	protected List<Inputs> requiredInputs = new LinkedList<Inputs>();
	protected List<Inputs> optionalInputs = new LinkedList<Inputs>();

	public ATool(String[] args, String className){
		helpOptions.addOption("h", "help", false, "show this help page");
		options.addOption("h", "help", false, "show this help page");
		options.addOption("o", true, "the output folder, use input folder as default");
		options.addOption(Option.builder("i")
				.longOpt("input")
				.argName("INPUT")
				.desc("The EAGER folder used for as input")
				.hasArg()
				.hasArgs()
				.build());
//		options.addOption(OptionBuilder.withLongOpt("input")
//				.withArgName("INPUT")
//				.withDescription("The EAGER folder used for as input")
//				.isRequired()
//				.hasArg()
//				.create("i"));
		requiredInputs.add(Inputs.input);

		HelpFormatter helpformatter = new HelpFormatter();

		CommandLineParser parser = new DefaultParser();

		try {
			CommandLine cmd = parser.parse(helpOptions, args);
			if (cmd.hasOption('h')) {
				helpformatter.printHelp(className, options);
				System.exit(1);
			}
		} catch (ParseException e1) {
		}
	}

	protected abstract void run();


	/**
	 * 
	 *
	 */
	protected void getVCFFiles() {
		this.vcfFilePaths = new TreeMap<String, String>();
		for(File sample: this.eagerSamples){
			String vcfFile = findVCFFile(sample);
			if(vcfFile.length() >0) {
				vcfFilePaths.put(sample.getName(), vcfFile);
			}else {
				System.out.println("no vcf file for sample "+sample);
			}
		}
	}

	/**
	 * 
	 * find all folders that are EAGER output folders
	 *
	 */
	protected void getEagerSamples() {
		this.eagerSamples = Utilities.getEagerSamples(new File(this.eagerFolder));
	}

	/**
	 * @param eagerFolder
	 * @return
	 * 
	 */
	protected String findGenomeFile(File eagerFolder) {
		String[] eagerContent = eagerFolder.list();
		String genomeFolder = "";
		for(String file: eagerContent){
			if(file.startsWith("12-")){
				genomeFolder=eagerFolder.getAbsolutePath()+"/"+file;
				break;
			}
		}
		if(genomeFolder.length()==0){
//			throw new Exception("Could not find VCF2Genome folder for sample:\n"+eagerFolder.getAbsolutePath());
			return "";
		}
		String  genomeFile = "";
		String[] genomeFolderContent = new File(genomeFolder).list();
		for(String file:genomeFolderContent){
			if(file.endsWith(".fasta") && !file.endsWith("nr1234.fasta") && !file.endsWith("refMod.fasta")){
				genomeFile=genomeFolder+"/"+file;
			}
		}
		if(genomeFile.length()==0){
//			throw new Exception("Could not find genome fasta file for sample:\n"+eagerFolder.getAbsolutePath());
			return "";
		}
		return genomeFile;
	}

	/**
	 * @param eagerFolder
	 * @return
	 * 
	 */
	protected String findVCFFile(File eagerFolder) {
		String[] eagerContent = eagerFolder.list();
		String genomeFolder = "";
		for(String file: eagerContent){
			if(file.startsWith("10-")){
				genomeFolder=eagerFolder.getAbsolutePath()+"/"+file;
				break;
			}
		}
		if(genomeFolder.length()==0){
			return "";//TODO
//			throw new Exception("Could not find VCF2Genome folder for sample:\n"+eagerFolder.getAbsolutePath());
		}
		String  vcfFile = "";
		String[] gatkFolderContent = new File(genomeFolder).list();
		for(String file:gatkFolderContent){
			if(file.endsWith(".vcf") || file.endsWith(".vcf.gz")){
				vcfFile=genomeFolder+"/"+file;
			}
		}
		if(vcfFile.length()==0){
			return "";//TODO
//			throw new Exception("Could not find vcf file for sample:\n"+eagerFolder.getAbsolutePath()+"\n"+vcfFile);
		}
		return vcfFile;
	}

	public List<Inputs> getRequiredInputs(){
		return this.requiredInputs;
	}

	public List<Inputs> getOptionalInputs(){
		return this.optionalInputs;
	}

}
