# MUSIAL - MUlti Sample varIant AnaLysis
Summarize SNV and indel information on single gene or genome level together with other relevant statistics based on _.vcf_ files.

## Dependencies:
- JDK 15+
- Gradle 7+

## Building:
A precombiled **.jar** can be found in "Releases". This program can be built with gradle (https://gradle.org).
For that just type `gradle clean build` in the projects root directory.
The **jar** file is then contained in `/releases` folder.

## Parameters
The parameters tools are available:
* `-a,--referenceAnnotation <arg>`
Path to a .gff file; The reference genome feature annotation.

* `-gf,--geneFeatures <arg>`
Text file specifying gene features to analyze exclusively instead of
the whole genome. Each line has to contain at least a string to match
with the NAME attribute of the specified .gff file and may contain
the following optional comma-separated fields in the given order:
<NAME> (a string representing the name to use for the feature)
<PDBPATH> (a string representing the path to a .pdb file specifying the protein product of the feature).

* `-h,--help`
Display help information.

* `-hf,--heterozygousFrequencies <arg>`
Two comma separated numbers in the interval [0.0,1.0] used for
heterozygous calls, i.e. a het. call will be accepted if its read
frequency lies in the specified interval. [0.45, 0.55]

* `-mc,--minCoverage <arg>`
Minimum coverage to call a SNV [5.0]

* `-mf,--minFrequency <arg>`
Minimum (read) frequency to call a SNV [0.9]

* `-mq,--minQual <arg>`
Minimum quality to call a SNV [30.0]
                                       
* `-nt,--numThreads <arg>`
Number of threads to use [1]
 
* `-o,--output <arg>`
Path to the directory at which results files shall be generated.
                                       
* `-r,--reference <arg>`
Path to a .fasta file; The reference genome sequence.
                                       
* `-s,--sampleSpecification <arg>`
Text file specifying the sample input files and sample meta-information.
Each line has to contain at least the path to the samples .vcf file and 
may contain the following optional comma-separated fields in the given order:
<NAME> (a string representing the name to use for the sample).
                                       
* `-sDir,--sampleDirectory <arg>`
Directory from which sample input files are collected.
The specified directory and all its sub-directories are scanned for
.vcf files which are collected as sample input files.
                                       
* `-sf,--snpEff`
Run SNPEff. [false]
 
* `-v,--verbose`
Verbose mode on, i.e. print stacktrace on errors. [false]

At least a reference genome (`-r`) and sample input (`-s` or `-sDir`) have to be specified.

## Examples:
An exemplary output directory can be found at `/examples/treponema_pallidum/BamA`.
The example output directory can also be used to test IVE - MUSIALs Integrative Visualization Extension.
Documentation on IVE will follow soon.