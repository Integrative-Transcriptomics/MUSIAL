# MUSIAL - MUlti Sample varIant AnaLysis

calculate SNP, gene, and whole genome alignments, together with other relevant statistics based on vcf files.

## Dependencies:

- jdk7+
- gradle

## generating the jar file
A precombiled jar can be found in the folder "Releases"
This program can be built with gradle (https://gradle.org). for that just type

`gradle build`

The jar-files are then contained in the build/libs folder. Additionally the jar file in the "Releases" folder will be replaced with the newly built jar.

## Parameters
The parameters tools are available:
- -ar,--addReference <arg>           add reference under given name to alignment
- -c,--coverage <arg>                min coverage to make call [5.0]
- -ca,--covAdd <arg>                 minCoverage to make additional call [5.0]
- -e,--exclude <EXCLUDE>             exclude given positions in Alignment
- -f,--frequency <arg>               minimum frequency to call a SNP [0.9]
- -g,--gff <arg>                     the gff file
- -gn,--geneNames <GENENames>        the gene names for gene alignments
- -h,--help                          show this help page
- -i,--input <INPUT>                 The EAGER folder used for as input
- -if,--inputFile <INPUTFILE>        the input vcf files from file
- -il,--inputList <INPUTLIST>        the input vcf files als list
- -o,--output <OUTPUT>               the output folder
- -ogf,--outgroupFile <arg>          vcfFile for outgroup
- -ogn,--outgroupName <arg>          name for outgroup (contained in input)
- -p,--phylogenetics                 calculate phylogenetic trees
- -r,--reference <REFERENCE>         the reference fasta file
- -s,--sharedAlleleFreq              write shared allele frequencies for each position
- -sf,--snpEff                       run SNPEff
- -t,--threads <arg>                 number of threads to use [1]
- -u,--uncovered                     write uncovered positions for each sample
- -y,--heterozygous <HETEROZYGOUS>   call heterozygous positions
                                    [0.45,0.55]
The parameters -o and -r, together with at least one input option has to be given (-i, -il, -if)

With the option -il, you can give the different vcf files separated with a space.  
With the option -if, the file has to contain one path to each vcf file per line.  
Both commands will use the name of the folder containing the vcf file as sample name. So make sure that no duplications exist.  
For the option -i, a new directory called "vcfs" will be created with the sample names of the EAGER samples and symbolic links will be created to the corresponding vcf files.
