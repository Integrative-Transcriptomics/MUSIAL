/**
 * 
 */
package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import datastructure.FastAEntry;
import datastructure.Gene;
import datastructure.SNPStatistics;
import datastructure.SNPTable;
import io.FastaReader;
import io.GFFParser;
import io.myVCFFileReader;
import main.Musial;
import multithreaded.MVcfReader;
import utilities.Pair;
import utilities.Utilities;

/**
 * @author Alexander Seitz
 *
 */
public class Analyzer extends ATool {
	
	private static final String CLASS_NAME = Musial.CLASS_NAME + " [options]";

	private Set<String> vcfFiles = new TreeSet<String>();
	private String outputFolder;
	private Map<String, String> sampleNames;
	private Map<String, String> outgroupSamples = new TreeMap<String, String>();
	private Map<String, myVCFFileReader> vcfFileReader;
	private Map<Integer, List<String>> lowCovRegions;
	private Map<Integer, Map<String, Integer>> deletions;
	private Map<String, SNPStatistics> statistics;
	private Integer threads = 1;

	private String gffFile;
	private String referenceFile;
	private Double minFreq = 0.9;
	private Double minCovGood = 5d;
	private Double minCovAdd = 5d;
	private Double genQual = 30d;

	private Double minHet = 0.45;
	private Double maxHet = 0.55;
	private Boolean callHeterozygous = false;
	
	private Boolean writeUncovered = false;
	private Boolean writeSharedAlleleFreqs = false;
	private Boolean addReference = false;
	private Boolean hasOutgroupName = false;
	private Boolean calculatePhylos = false;
	private Boolean runSnpEff = false;

//	private String snpEffDir = "";
	private String snpEffName = "mySample";
	private String snpEffOutputFile = "";
//	private File snpEffJar;
	String snpEffInput = "";
	
	private String referenceAlignmentName = "";
	private String outGroupName = "";

	private Set<Integer> toExclude = new TreeSet<Integer>();
	
	private Set<String> geneNames = new TreeSet<String>();
	private Set<Gene> genes = new HashSet<Gene>();

	private SNPTable snpTable;

	private FastaReader fr;

	private boolean writeVCFFileListFile = false;
	
	private Integer finished = 0;
	
//	private List<String> differences = new LinkedList<String>();

	/**
	 * @param args
	 * @param className
	 */
	public Analyzer(String[] args) {
		super(args, CLASS_NAME);
		Options helpOptions = new Options();
//		Options options = new Options();

		helpOptions.addOption("h", "help", false, "show this help page");

		options.addOption("h", "help", false, "show this help page");
		options.addOption("t", "threads", true, "number of threads to use ["+this.threads+"]");
		options.addOption("c", "coverage", true, "min coverage to make call ["+this.minCovGood+"]");
		options.addOption("ca", "covAdd", true, "minCoverage to make additional call ["+this.minCovAdd+"]");
		//		options.addOption("c", "compare", true, "compare to other SNPTable");
		options.addOption("f", "frequency", true, "minimum frequency to call a SNP ["+this.minFreq+"]");
		options.addOption("g", "gff", true, "the gff file");
		options.addOption("u", "uncovered", false, "write uncovered positions for each sample");
		options.addOption("s", "sharedAlleleFreq", false, "write shared allele frequencies for each position");
		options.addOption("ar", "addReference", true, "add reference under given name to alignment");
		options.addOption("ogf", "outgroupFile", true, "vcfFile for outgroup");
		options.addOption("ogn", "outgroupName", true, "name for outgroup (contained in input)");
		options.addOption("p", "phylogenetics", false, "calculate phylogenetic trees");
		options.addOption("sf", "snpEff", false, "run SNPEff");
		options.addOption(Option.builder("gn")
				.longOpt("geneNames")
				.argName("GENENames")
				.desc("the gene names for gene alignments")
				.hasArg()
				.hasArgs()
				.build());
		options.addOption(Option.builder("y")
				.longOpt("heterozygous")
				.argName("HETEROZYGOUS")
				.desc("call heterozygous positions ["+this.minHet+","+this.maxHet+"]")
				.numberOfArgs(2)
				.optionalArg(true)
				.build());
		options.addOption(Option.builder("il")
				.longOpt("inputList")
				.argName("INPUTLIST")
				.desc("the input vcf files als list")
				.hasArg()
				.hasArgs()
				.build());
		options.addOption(Option.builder("if")
				.longOpt("inputFile")
				.argName("INPUTFILE")
				.desc("the input vcf files from file")
				.hasArg()
				.hasArgs()
				.build());
		options.addOption(Option.builder("o")
				.longOpt("output")
				.argName("OUTPUT")
				.desc("the output folder")
				.required()
				.hasArg()
				.build());
		options.addOption(Option.builder("r")
				.longOpt("reference")
				.argName("REFERENCE")
				.desc("the reference fasta file")
				.required()
				.hasArg()
				.build());
		options.addOption(Option.builder("e")
				.longOpt("exclude")
				.argName("EXCLUDE")
				.desc("exclude given positions in Alignment")
				.hasArg()
				.hasArgs()
				.build());
		

		HelpFormatter helpformatter = new HelpFormatter();

		CommandLineParser parser = new DefaultParser();

		try {
			CommandLine cmd = parser.parse(helpOptions, args);
			if (cmd.hasOption('h')) {
				printHelp(options, helpformatter);
				System.exit(1);
			}
		} catch (ParseException e) {
			// do nothing
		}

		try{
			CommandLine cmd = parser.parse(options, args);
			this.referenceFile = cmd.getOptionValue("r");
			this.outputFolder = cmd.getOptionValue("o");
			File out = new File(this.outputFolder);
			if(!out.exists()){
				out.mkdirs();
			}
			if(!out.isDirectory()){
				System.err.println("output no directory");
				System.exit(1);
			}
			this.vcfFiles = new TreeSet<String>();
			if(cmd.hasOption("il")) {
				for(String f: cmd.getOptionValues("il")) {
					this.vcfFiles.add(f);
				}
			}
			if(cmd.hasOption("if")) {
				for(String f: cmd.getOptionValues("if")) {
					parseInputVcfFile(f);
				}
			}
			if(cmd.hasOption("i")) {
				for(String eagerFolder: cmd.getOptionValues("i")) {
					this.eagerFolder = eagerFolder;
					try {
						getEagerSamples();
						getVCFFiles();
						prepareOutput();
						this.writeVCFFileListFile  = true;
					} catch (Exception e) {
						e.printStackTrace();
					}
					
				}
			}
			if(this.vcfFiles.size() == 0) {
				System.err.println("no input files given");
				printHelp(options, helpformatter);
				System.exit(1);
			}
			if(cmd.hasOption("p")) {
				this.calculatePhylos = true;
			}
			if(cmd.hasOption("ogf")) {
				String outGroupFile = cmd.getOptionValue("ogf");
				if(new File(outGroupFile).exists()) {
					this.outgroupSamples.put(new File(outGroupFile).getParentFile().getName(), new File(outGroupFile).getAbsolutePath());
				}
				
			}
			if(cmd.hasOption("ogn")) {
				this.outGroupName = cmd.getOptionValue("ogn");
				this.hasOutgroupName = true;
			}
			if(cmd.hasOption("ar")) {
				this.addReference = true;
				this.referenceAlignmentName = cmd.getOptionValue("ar");
			}
			if(cmd.hasOption("s")) {
				this.writeSharedAlleleFreqs = true;
			}
			if(cmd.hasOption("u")) {
				this.writeUncovered = true;
			}
			if(cmd.hasOption("g")) {
				this.gffFile = cmd.getOptionValue("g");
			}
			if(cmd.hasOption("gn")) {
				if(!cmd.hasOption("g")) {
					System.err.println("Gene alignments only with option -g possible");
					System.err.println("skipping gene alignments");
				}else {
					for(String g: cmd.getOptionValues("gn")) {
						this.geneNames.add(g.trim());
					}
				}
			}
			if(cmd.hasOption("t")){
				if(Utilities.isInteger(cmd.getOptionValue("t"))){
					this.threads = Integer.parseInt(cmd.getOptionValue("t"));
				}else{
					System.err.println("given number of threads no Integer:\t"+cmd.getOptionValue("t"));
					System.err.println("running on "+this.threads+" threads");
				}
			}
			if(cmd.hasOption("e")){
				this.toExclude = new TreeSet<Integer>();
				for(String exclude: cmd.getOptionValues("e")){
					if(new File(exclude).exists()){
						List<String> excludeList = Utilities.readFileLinewise(exclude);
						for(String excl: excludeList){
							if(Utilities.isInteger(excl)){
								this.toExclude.add(Integer.parseInt(excl));
							}else{
								System.err.println("Not an Integer: "+excl);
								System.err.println("in file: "+exclude);
							}
						}
					}else if(Utilities.isInteger(exclude)){
						this.toExclude.add(Integer.parseInt(exclude));
					}else{
						System.err.println("could not use exclude information: "+exclude);
					}
				}
			}
			if(cmd.hasOption("c")) {
				if(Utilities.isInteger(cmd.getOptionValue("c"))) {
					this.minCovGood = (double) Integer.parseInt(cmd.getOptionValue("c"));
					this.minCovAdd = this.minCovGood;
				}else{
					System.err.println("Given value for coverage not Integer: "+cmd.getOptionValue("c"));
				}
			}
			if(cmd.hasOption("ca")) {
				if(Utilities.isInteger(cmd.getOptionValue("c"))) {
					this.minCovAdd = (double) Integer.parseInt(cmd.getOptionValue("ca"));
				}else{
					System.err.println("Given value for additional coverage not Integer: "+cmd.getOptionValue("ca"));
				}
			}
			//			if(cmd.hasOption("c")) {
			//				this.compareTo = cmd.getOptionValue("c");
			//			}else {
			//				this.compareTo = "";
			//			}
			if(cmd.hasOption("f")) {
				if(Utilities.isDouble(cmd.getOptionValue("f"))) {
					this.minFreq = Double.parseDouble(cmd.getOptionValue("f"));
					if(this.minFreq <= 0.0 || this.minFreq > 1) {
						System.err.println("Wrong parameter for minimum frequency: "+this.minFreq);
						System.err.println("please choose a value between (0.0,1.0]");
						System.exit(1);
					}
				}else {
					System.err.println("Value for minimum frequency no Double: "+this.minFreq);
					System.err.println("please choose a value between (0.0,1.0]");
					System.exit(1);
				}
			}
			if(cmd.hasOption("y")) {
				this.callHeterozygous = true;
				String[] percentages = cmd.getOptionValues("y");
				if(percentages != null) {
					if(percentages.length != 0 || percentages.length != 2) {
						System.err.println("Wrong number of parameters for option: -y");
						System.err.println("given: "+Arrays.asList(percentages));
						System.err.println("Call: -y <minPercentage> <maxPercentage>");
						System.err.println("Using default values: "+this.minHet+","+this.maxHet);
					}else if(percentages.length == 2) {
						if(Utilities.isDouble(percentages[0]) && Utilities.isDouble(percentages[1])) {
							Double min = Double.parseDouble(percentages[0]);
							Double max = Double.parseDouble(percentages[1]);
							if(min > 1 || max > 1 || min < 0 || max < 0) {
								System.err.println("Given percentages not within 0 and 1 for option: -y");
								System.err.println("Given: "+Arrays.asList(percentages));
								System.err.println("Using default values: "+this.minHet+","+this.maxHet);
							}
							if(min > max) {
								System.err.println("minimum percentage larger than maximum for option: -y");
								System.err.println("Given: "+Arrays.asList(percentages));
								System.err.println("Switching them around");
								double tmp = min;
								min = max;
								max = tmp;
								System.err.println("Using: "+min+","+max);
							}
							this.minHet = min;
							this.maxHet = max;
						}else {
							System.err.println("One of the given percentages no double for option: -y");
							System.err.println("Given: "+Arrays.asList(percentages));
							System.err.println("Using default values: "+this.minHet+","+this.maxHet);
						}
					}
				}
			}
			if(cmd.hasOption("sf")) {
				this.runSnpEff = true;
				if(!cmd.hasOption("g")) {
					System.err.println("For the SNPEff analysis you need to input a gff (-g)");
					System.exit(1);
				}
			}
		} catch (ParseException e) {
			printHelp(options, helpformatter);
			System.err.println(e.getMessage());
			System.exit(1);
		}
		run();
	}
	
	private void prepareOutput() {
		if(this.eagerFolder.equals(this.outputFolder)){
			this.outputFolder += "/MVCF";
			Utilities.createOutFolder(this.outputFolder);
		}
		String vcfFolder = this.outputFolder+"/vcfs";
		addVcfFilesToFolder(vcfFolder);
	}
	
	private void addVcfFilesToFolder(String vcfFolder) {
		for(String sampleName: this.vcfFilePaths.keySet()){
			String originalVCF = this.vcfFilePaths.get(sampleName);
			String vcfFile = vcfFolder+"/"+sampleName;
			Utilities.createOutFolder(vcfFile);
			vcfFile += "/input.vcf";
			if(originalVCF.endsWith(".gz"))
				vcfFile += ".gz";
			Utilities.generateSymlink(originalVCF, vcfFile);
			this.vcfFiles.add(vcfFile);
		}
	}

	/**
	 * 
	 */
	private void parseInputVcfFile(String f) {
		try {
			@SuppressWarnings("resource")
			BufferedReader br = new BufferedReader(new FileReader(f));
			String line = "";
			while((line = br.readLine()) != null) {
				if(line.length()>0) {
					this.vcfFiles.add(line.trim());
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param options
	 * @param helpformatter
	 */
	private void printHelp(Options options, HelpFormatter helpformatter) {
		helpformatter.printHelp(Analyzer.CLASS_NAME, options);
		System.err.println("At least one input option has to be given (-i, -il, -if)");
	}

	/* (non-Javadoc)
	 * @see tools.ATool#run()
	 */
	@Override
	protected void run() {

		// create the output directory
		File out = new File(this.outputFolder);
		String outdir = out.getAbsolutePath();
		this.outputFolder = outdir;
		Utilities.createOutFolder(outdir);
		// if desired, write the input file list
		if(this.writeVCFFileListFile) {
			StringBuffer vcfFileListFile = new StringBuffer();
			for(String vcfFile: this.vcfFiles) {
				vcfFileListFile.append(vcfFile);
				vcfFileListFile.append("\n");
			}
			Utilities.writeToFile(vcfFileListFile.toString(), new File(this.outputFolder+"/vcfFilesUsed.txt"));
		}

		// set variables

		//				this.differences = new LinkedList<String>();
		this.fr = new FastaReader(this.referenceFile);
//		boolean firstRun = true;

		for(String header: this.fr.getHeaders()) {
			this.vcfFileReader = new HashMap<String, myVCFFileReader>();
			this.sampleNames = new TreeMap<String,String>();
			
			System.out.println(header);

			// parse the vcf files

			System.out.println("parsing vcf files");
			long start = System.currentTimeMillis();
			parseVCFFiles(fr.getEntry(header).getSequenceLength(), header);
//			parseVCFFiles(fr.getFirstEntry().getSequenceLength());
			long end = System.currentTimeMillis();
			System.out.println("finished parsing vcf files");
			System.out.println("\ttime:\t"+ ((end-start)/1000) + "s");

			// calculate the differences TODO maybe add
			//		if(this.compareTo.length()>0 && new File(this.compareTo).exists()) {
			//			System.out.println("calculating differences");
			//			start = System.currentTimeMillis();
			//			this.differences = this.snpTable.compareTo(new SNPTable("/share/projects/Lepra/results/All_newMappingSameParameters/multivcf/2017_02_08_MVCF_results_GC_comb_newSamples/snpTable.tsv"));
			//			this.differences = addAdditionalInformationToDifferences();
			//			end = System.currentTimeMillis();
			//			System.out.println("finished calculating differences");
			//			System.out.println("\ttime:\t"+ ((end-start)/1000) + "s");			
			//		}

			// if set, get the desired genes

			if(this.geneNames.size()>0) {
				try {
					genes = new GFFParser(this.gffFile).getGenes(this.geneNames);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}

			// write the output

			System.out.println("generating output");
			start = System.currentTimeMillis();
			writeOutput(header);
			end = System.currentTimeMillis();
			System.out.println("finished writing output");
			System.out.println("\ttime:\t"+ ((end-start)/1000) + "s");

			//TODO

			if(this.runSnpEff){
				System.out.println("running SNPEff");
				start = System.currentTimeMillis();
				// download snpEff
				System.out.println("\tDownloading SNPEff");
				String[] downloadSnpEff = {"wget", "http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip"};
				Utilities.runCommand(downloadSnpEff, this.outputFolder+"/download.error", this.outputFolder+"/download.output", this.outputFolder);
				String snpEffZip = this.outputFolder+"/snpEff_latest_core.zip";
				// if download was unsuccessful, copy from internal jar
				if(!new File(snpEffZip).exists()) {
					Utilities.ExportResource("snpEff_latest_core.zip", snpEffZip);
				}
				// extract the SNPEff zip file
				System.out.println("\tExtracting SNPEff");
				String[] unpackSnpEff = {"unzip", snpEffZip};
				Utilities.runCommand(unpackSnpEff, "", "", this.outputFolder);
				String snpEffDir = this.outputFolder+"/snpEff";
				System.out.println("\tCreating SNPEff database");
				createSNPEffDatabase(snpEffDir);
				runSnpEff(snpEffDir);
				// create new SNPTable with SNPEff infos
				Utilities.writeToFile(this.snpTable.toStringWithSNPEffInfos(this.snpEffOutputFile), new File(this.outputFolder+"/snpTableWithSnpEffInfos.tsv"));
				end = System.currentTimeMillis();
				System.out.println("\ttime:\t"+ ((end-start)/1000) + "s");
			}
//			firstRun = false;
		}
	}

	/**
	 * @param snpEffDir
	 */
	private void createSNPEffDatabase(String snpEffDir) {
		// create data dir
		Utilities.createOutFolder(snpEffDir+"/data");
		// create directory for current sample
		Utilities.createOutFolder(snpEffDir+"/data/"+this.snpEffName);
		// copy the gff file
		Utilities.copyFileUsingIO(new File(this.gffFile), new File(snpEffDir+"/data/"+this.snpEffName+"/genes.gff"));
		// copy the reference file
		Utilities.copyFileUsingIO(new File(this.referenceFile), new File(snpEffDir+"/data/"+this.snpEffName+"/sequences.fa"));
		// add line to config
		try {
			List<String> lines = Files.readAllLines(new File(snpEffDir+"/snpEff.config").toPath(), StandardCharsets.UTF_8);
			lines.add(126, this.snpEffName+".genome : "+this.snpEffName);
			Files.write(new File(snpEffDir+"/snpEff.config").toPath(), lines, StandardCharsets.UTF_8);
		} catch (IOException e) {
			e.printStackTrace();
		}
		// create database
		String[] createSNPEffDatabase = {"java", "-jar", "snpEff.jar", "build", "-gff3", "-v", this.snpEffName};
		Utilities.runCommand(createSNPEffDatabase, "", "", snpEffDir);
		
	}

	private void runSnpEff(String snpEffDir) {
		this.snpEffOutputFile = this.outputFolder+"/snpeff_output.out";
		String[] runSNPEff = {"java", "-jar", "snpEff.jar", this.snpEffName, this.snpEffInput};
		Utilities.runCommand(runSNPEff, this.outputFolder+"/snpEff.error", this.snpEffOutputFile, snpEffDir);
	}
	
	private void writeOutput(String header) {
		// create the output folder
		File out = new File(this.outputFolder);
		String outdir = out.getAbsolutePath();
		this.outputFolder = outdir;
		// write out the files with the used vcf files
		Utilities.createOutFolder(outdir);
		if(this.writeVCFFileListFile) {
			StringBuffer vcfFileListFile = new StringBuffer();
			for(String sampleName: this.vcfFilePaths.keySet()) {
				vcfFileListFile.append(this.vcfFilePaths.get(sampleName));
				vcfFileListFile.append("\n");
			}
			Utilities.writeToFile(vcfFileListFile.toString(), new File(this.outputFolder+"/vcfFilesUsed.txt"));
		}
		// write the snpTable
		Utilities.writeToFile(this.snpTable.toString(), new File(outdir+"/"+header+"_snpTable.tsv"));
		
		// write the snpTable with only heterozygous calls
		if(this.callHeterozygous) {
			Utilities.writeToFile(this.snpTable.toStringHeterozygousOnly(), new File(outdir+"/"+header+"_snpTable_only_heterozybous.tsv"));
		}
		
		// write the vcf file for SNPEff
		File snpEffVcfFile = new File(outdir+"/"+header+"_snpVcfForSnpEff.vcf");
		this.snpEffInput = snpEffVcfFile.getAbsolutePath();
		Utilities.writeToFile(this.snpTable.toSnpVcfForSnpEffString(this.fr.getEntry(header).getHeader()), snpEffVcfFile);
		
		// write the SNPalignment File
		File snpAlignmentFile = new File(outdir+"/"+header+"_snpAlignment.fasta");
		Utilities.writeToFile(this.snpTable.toAlignment(this.deletions), snpAlignmentFile);
		
		// write the genome alignment File
		File genomeAlignmentFile = new File(outdir+"/"+header+"_genomeAlignment.fasta");
		writeGenomeAlignments(genomeAlignmentFile, header);
		
		// write gene alignments
		List<File> geneAlignmentFiles = new LinkedList<File>();
		for(Gene gene: this.genes) {
			if(gene.getSource().equals(header)) {
				File currGeneAlignmentFile = new File(outdir+"/"+header+"_gene_"+gene.getName()+"_Alignment.fasta");
				writeGeneAlignment(currGeneAlignmentFile, gene.getStart(), gene.getEnd(), header);
				geneAlignmentFiles.add(currGeneAlignmentFile);
			}
		}
		writeLowCoveragePositions(new File(outdir+"/"+header+"_lowCoveragePositions.tsv"));
		File excludedSNPAlignmentFile = null;
		if(this.toExclude.size()>0){
			excludedSNPAlignmentFile = new File(outdir+"/"+header+"_snpAlignment_excluded.fasta");
			Utilities.writeToFile(this.snpTable.toAlignment(this.toExclude, this.deletions), excludedSNPAlignmentFile);
		}
		if(this.writeUncovered) {
			for(String sample: vcfFileReader.keySet()) {
				writeUncovered(new File(outdir+"/"+header+"_uncovered_"+sample+".bed"), sample);
			}
		}
		if(writeSharedAlleleFreqs) {
			writeAlleleFrequencies(new File(outdir+"/"+header+"_alleleFrequencies.tsv"), header);
		}
//		Utilities.writeToFile(gatherDifferences(), new File(outdir+"/differences.txt"));
		Utilities.writeToFile(gatherStatistics(), new File(outdir+"/"+header+"_snpStatistics.tsv"));
		Map<Integer, String> heterozygousCalls = getHeterozygousCalls();
		if(heterozygousCalls.size()>0) {
			StringBuffer hetCallsOutput = new StringBuffer();
			for(Integer pos: heterozygousCalls.keySet()) {
				hetCallsOutput.append(pos);
				hetCallsOutput.append("\t");
				hetCallsOutput.append(heterozygousCalls.get(pos));
				hetCallsOutput.append("\n");
			}
			Utilities.writeToFile(hetCallsOutput.toString(), new File(outdir+"/"+header+"_heterozygousCalls.tsv"));
		}
		if(this.calculatePhylos) {
			System.out.println("\tcalculating phylogenies");
			long startPhylo = System.currentTimeMillis();
			calculatePhylogeny(snpAlignmentFile);
//			calculatePhylogeny(genomeAlignmentFile);
//			for(File f: geneAlignmentFiles) {
//				calculatePhylogeny(f);
//			}
//			if(excludedSNPAlignmentFile != null) {
//				calculatePhylogeny(excludedSNPAlignmentFile);
//			}
			long endPhylo = System.currentTimeMillis();
			System.out.println("\tfinished calculating phylogenies");
			System.out.println("\t\ttime:\t"+ ((endPhylo-startPhylo)/1000) + "s");
		}
	}
	
	/**
	 * @param snpTableFile
	 */
	private void calculatePhylogeny(File alignmentFile) {
		String[] runPhylogeny = new String[] {"raxml-ng-mpi", "--all", "--msa", alignmentFile.getAbsolutePath(), "--model", "GTR+G", "--tree", "pars{10}", "--bs-trees", "200"};
		Process process;
		try {
			process = new ProcessBuilder(runPhylogeny).redirectError(new File(alignmentFile.getAbsolutePath()+".error")).redirectOutput(new File(alignmentFile.getAbsolutePath()+".output")).start();
			process.waitFor();
		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
		}
		
	}

	/**
	 * @param file
	 */
	private void writeAlleleFrequencies(File file, String sequenceHeader) {
		StringBuffer result = new StringBuffer();
		result.append("Position\tFrequency\n");
		for(int i=1; i<=this.fr.getEntry(sequenceHeader).getSequence().length(); i++) {
			result.append(i);
			result.append("\t");
			Map<String, Integer> freq = new HashMap<String, Integer>();
			for(String sample: this.vcfFileReader.keySet()) {
				String base = this.vcfFileReader.get(sample).get(i);
				if(freq.containsKey(base)) {
					freq.put(base, freq.get(base)+1);
				}else {
					freq.put(base, 1);
				}
			}
			if(freq.size() == 1) {
				result.append(1d);
			}else {
				double max = 0;
				double all = 0;
				for(int value: freq.values()) {
					all += value;
					if(value > max) {
						max = value;
					}
				}
				result.append(max/all);
			}
			result.append("\n");
		}
		Utilities.writeToFile(result.toString(), file);
		
	}

	private void writeUncovered(File file, String sampleName) {
		StringBuffer result = new StringBuffer();
		List<Pair<Integer,Integer>> uncovered = this.vcfFileReader.get(sampleName).getUncovered();
		String refName = this.vcfFileReader.get(sampleName).getContigName();
		for(Pair<Integer, Integer> p: uncovered) {
			result.append(refName+"\t"+p.getFirst()+"\t"+p.getSecond()+"\n");
		}
		Utilities.writeToFile(result.toString(), file);
	}
	
	private void writeGeneAlignment(File file, Integer start, Integer end, String sequenceHeader) {
		StringBuffer alignment = new StringBuffer();
		for(String sample:this.vcfFileReader.keySet()){
			FastAEntry fe = new FastAEntry(sample, this.vcfFileReader.get(sample).getGeneSequence(this.fr.getEntry(sequenceHeader), this.deletions, start, end));
			alignment.append(fe.toString());
		}
		Utilities.writeToFile(alignment.toString(), file);
	}
	
	private void writeGenomeAlignments(File file, String sequenceHeader) {
		StringBuffer alignment = new StringBuffer();
		for(String sample:this.vcfFileReader.keySet()){
			FastAEntry fe = new FastAEntry(sample, this.vcfFileReader.get(sample).getGenomeSequence(this.fr.getEntry(sequenceHeader), this.deletions));
			alignment.append(fe.toString());
		}
		Utilities.writeToFile(alignment.toString(), file);
	}
	
	private void writeLowCoveragePositions(File file) {
		List<StringBuffer> regions = new LinkedList<StringBuffer>();
		int regionID = 0;
		Integer lastPos = 1;
		StringBuffer currRegion = new StringBuffer();
		for(Integer pos: this.lowCovRegions.keySet()){
			if(this.lowCovRegions.get(pos).size() >3){
				if(pos != lastPos+1){
					regionID++;
					regions.add(currRegion);
					currRegion = new StringBuffer();
				}
				lastPos = pos;
				currRegion.append(pos);
				currRegion.append("\t");
				currRegion.append(regionID);
				currRegion.append("\t");
				currRegion.append(this.lowCovRegions.get(pos).size());
				currRegion.append("\t");
				currRegion.append(this.sampleNames.size());
				currRegion.append("\t");
				for(String sample: this.lowCovRegions.get(pos)){
					currRegion.append(sample);
					currRegion.append(":");
					currRegion.append(this.vcfFileReader.get(sample).getCoverage(pos));
					currRegion.append(",");
				}
				currRegion.append("\n");
			}
		}
		StringBuffer result = new StringBuffer();
		result.append("position\tregionID\tnumLowCovSamples\tnumSamples\tSampleNames\n");
		for(StringBuffer sb: regions){
			if(sb.toString().split("\n").length>=10){
				result.append(sb.toString());
			}
		}
		Utilities.writeToFile(result.toString(), file);
	}
	
	private String gatherStatistics() {
		StringBuffer result = new StringBuffer();
		result.append("SNP Statistics for "+this.vcfFileReader.size()+" samples\n");
		result.append("Coverage Threshold: "+this.minCovGood+"\n");
		result.append("Minimum SNP allele frequency: "+this.minFreq+"\n"
				);
		result.append("sample\tSNP Calls (all)\tSNP Calls (het)\tCoverage (fold)\tCoverage(%)\tReference Calls\ttotal Calls\tno Calls\n");
		for(String sample: this.statistics.keySet()) {
			result.append(this.statistics.get(sample).generateLine());
			result.append("\n");
		}
		return result.toString();
	}
	
	private Map<Integer, String> getHeterozygousCalls() {
		Map<Integer, String> heterozygousCalls = new TreeMap<Integer, String>();
		for(String sample: this.vcfFileReader.keySet()) {
			Set<Integer> currHetPos = this.vcfFileReader.get(sample).getHeterozygousPos();
			for(Integer pos: currHetPos) {
				String hetCall = this.vcfFileReader.get(sample).get(pos);
				String samples = sample+"("+hetCall+")";
				if(heterozygousCalls.containsKey(pos)) {
					samples = heterozygousCalls.get(pos)+"\t"+sample+"("+hetCall+")";
				}
				heterozygousCalls.put(pos, samples);
			}
		}
		return heterozygousCalls;
	}
	
	private void parseVCFFiles(Integer genomeSize, String sequenceHeader) {
		// check for duplicate sample names
//		if(firstRun) {
			System.out.println("There are "+this.vcfFiles.size()+" input samples");
			for(String vcf: this.vcfFiles){
				vcf = new File(vcf).getAbsolutePath();
				String sampleName = new File(vcf).getParentFile().getName();
				if(this.hasOutgroupName && sampleName.equals(this.outGroupName)) {
					this.outgroupSamples.put(sampleName, vcf);
					continue;
				}
				if(this.sampleNames.containsKey(sampleName)){
					System.err.println("Error: duplicate sample name: "+sampleName);
					System.exit(1);
				}
				this.sampleNames.put(sampleName, vcf);
			}
//		}
		Set<Integer> variablePositions = new TreeSet<Integer>();
		this.deletions = new TreeMap<Integer, Map<String, Integer>>();
		//		Integer numSamples = this.sampleNames.size();
		Integer currSample = 0;
		// parse each vcf file
		List<MVcfReader> todo = new LinkedList<MVcfReader>();
		// add all to list to run in parallel
		for(String sampleName: this.sampleNames.keySet()){
			if(!new File(this.sampleNames.get(sampleName)).exists()) {
				System.err.println("File for sample "+sampleName+" does not exist:");
				System.err.println(this.sampleNames.get(sampleName));
				System.err.println("coninuing without");
				continue;
			}
			myVCFFileReader reader = new myVCFFileReader(this.sampleNames.get(sampleName), sampleName, genomeSize, minCovGood, minCovAdd, minFreq, genQual, callHeterozygous, minHet,maxHet, sequenceHeader);
			this.vcfFileReader.put(sampleName, reader);
			todo.add(new MVcfReader(this, reader));
		}
		ExecutorService es = Executors.newFixedThreadPool(this.threads);
		try {
			@SuppressWarnings("unused")
			List<Future<MVcfReader>> answers = es.invokeAll(todo);
			es.shutdown();
		} catch (Exception e) {
			e.printStackTrace();
		}
		// add reference
		if(this.addReference) {
			if(this.sampleNames.containsKey(this.referenceAlignmentName)) {
				System.err.println("Name given for reference already exists.");
				System.err.println("Reference will not be added");
			}else {
				myVCFFileReader reader = new myVCFFileReader(this.referenceAlignmentName, this.fr.getEntry(sequenceHeader).getSequence(), genomeSize, minCovGood);
				this.sampleNames.put(referenceAlignmentName, "");
				this.vcfFileReader.put(referenceAlignmentName, reader);
			}
		}
		for(MVcfReader mReader: todo){
			currSample++;
			variablePositions.addAll(mReader.getReader().getVariablePositions());
			for(Pair<Integer, Integer> deletePos: mReader.getReader().getDeletions()) {
				Map<String, Integer> samplesWithCurrDeletion = new HashMap<String, Integer>();
				if(deletions.containsKey(deletePos.getFirst())) {
					samplesWithCurrDeletion = deletions.get(deletePos.getFirst());
				}
				samplesWithCurrDeletion.put(mReader.getReader().getSampleName(), deletePos.getSecond());
				deletions.put(deletePos.getFirst(), samplesWithCurrDeletion);
			}
		}
		// resolve maybe positions
		System.out.println("resolving possible maybe calls");
		resolveMaybeCalls();

		System.out.println("generating SNP table");

		// generating SNPTable and statistics
		this.statistics = new TreeMap<String, SNPStatistics>();

		Map<String, Map<Integer,String>> table = new TreeMap<String, Map<Integer,String>>();
		for(String sampleName: this.sampleNames.keySet()){
			myVCFFileReader reader = this.vcfFileReader.get(sampleName);
			if(reader == null) {
				System.out.println("Error: "+sampleName);
				System.exit(1);
			}
			Map<Integer, String> snps = new TreeMap<Integer, String>();
			for(Integer pos: variablePositions){
				String snp = reader.get(pos);
				snps.put(pos, snp);
			}
			table.put(sampleName, snps);
			SNPStatistics snpStatistics = new SNPStatistics(sampleName, reader.getNumSNPs(), reader.getNumHetSNPs(), reader.getTotalCalls(), reader.getRefCalls(), reader.getNoCalls(), reader.getFoldCoverage(), reader.getPercentCoverage());
			this.statistics.put(sampleName, snpStatistics);
		}
		this.snpTable = new SNPTable(table, this.fr.getEntry(sequenceHeader));
		for(String outgroupName: this.outgroupSamples.keySet()) {
			System.out.println("adding outgroup: "+outgroupName);
			String outgroupFile = this.outgroupSamples.get(outgroupName);
			myVCFFileReader reader  = new myVCFFileReader(outgroupFile, outgroupName, genomeSize, minCovGood, minCovAdd, minFreq, genQual, callHeterozygous, minHet,maxHet, sequenceHeader);
			reader.parseVCFFile();
			this.snpTable.add(outgroupName, reader);
		}

		System.out.println("generating low coverage regions");

		this.lowCovRegions = new TreeMap<Integer, List<String>>();
		for(String sampleName: this.sampleNames.keySet()){
			myVCFFileReader reader = this.vcfFileReader.get(sampleName);
			for(Integer pos: reader.getLowCoverageRegions()){
				List<String> lowCovSample = new LinkedList<String>();
				if(this.lowCovRegions.containsKey(pos)){
					lowCovSample = this.lowCovRegions.get(pos);
				}
				lowCovSample.add(sampleName);
				this.lowCovRegions.put(pos, lowCovSample);
			}
		}
		//		
		//		Map<String, Map<Integer,String>> genometable = new TreeMap<String, Map<Integer,String>>();
		//		for(String sampleName: this.sampleNames.keySet()){
		//			myVCFFileReader reader = this.vcfFileReader.get(sampleName);
		//			Map<Integer, String> snps = new TreeMap<Integer, String>();
		//			for(Integer pos=0; pos<reader.getGenomeSize(); pos++){
		//				String snp = reader.get(pos);
		//				snps.put(pos, snp);
		//			}
		//			genometable.put(sampleName, snps);
		//		}
		//		this.genomeTable = new SNPTable(genometable, this.fr.getFirstEntry());
		//		
		//		System.out.println("generated genome Table");
	}
	
	private void resolveMaybeCalls() {
		for(String sample: this.vcfFileReader.keySet()){
			Map<Integer, String> maybeCalls = this.vcfFileReader.get(sample).getMaybeGenotypes();
			for(Integer pos: maybeCalls.keySet()){
				for(String otherSample: this.vcfFileReader.keySet()){
					if(sample.equals(otherSample)){
						continue;
					}
					if(this.vcfFileReader.get(otherSample).canResolveCall(pos, maybeCalls.get(pos))){
						this.vcfFileReader.get(sample).resolveCall(pos, maybeCalls.get(pos));
						System.out.println("resolved maybe call:\t" + sample + "\t" + pos + "\t" + maybeCalls.get(pos));
					}
				}
			}
		}
	}

	/**
	 * @param sampleName
	 * @param timeInS 
	 * @return 
	 */
	public synchronized void finished(String sampleName, long timeInS) {
		this.finished++;
		System.out.println("\tfinished "+this.finished+"/"+this.vcfFiles.size()+"\t"+sampleName+"\t"+timeInS);
		
	}

}
