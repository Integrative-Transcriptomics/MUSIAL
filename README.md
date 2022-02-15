# **MUSIAL - MUlti Sample varIant AnaLysis**

`Contact:   simon.hackl@uni-tuebingen.de`

---
## **Description**:
Summarize SNV and indel information on single gene or genome level together with other relevant statistics based on **.vcf** files. Optional, assessed variant information can be allocated to reference protein structures (for single gene analysis). A companion tool to visualize allocated variants can be obtained [here](https://github.com/Integrative-Transcriptomics/MUSIAL-IVE).

---
## **Dependencies and Building**:
- JDK 15+
- Gradle 7+

A precombiled **.jar** can be found in "Releases". This program can be built with gradle (https://gradle.org).
For that just type `gradle clean build` in the projects root directory.
The **.jar** file is then contained in the `/releases` directory.

---
## **Usage**:
    usage: MUSIAL [OPTIONS]
     -a,--referenceAnnotation <arg>        Path to a .gff file; The reference
                                           genome feature annotation.
     -gf,--geneFeatures <arg>              Text file specifying gene features
                                           to analyze exclusively instead of
                                           the whole genome. Each line has to
                                           contain at least a string to match
                                           with the NAME attribute of the
                                           specified .gff file and may contain
                                           the following optional
                                           comma-separated fields in the given
                                           order:
                                           <NAME> (a string representing the
                                           name to use for the feature)
                                           <PDBPATH> (a string representing
                                           the path to a .pdb file specifying
                                           the protein product of the
                                           feature).
     -h,--help                             Display help information.
     -hf,--heterozygousFrequencies <arg>   Two comma separated numbers in the
                                           interval [0.0,1.0] used for
                                           heterozygous calls, i.e. a het.
                                           call will be accepted if its read
                                           frequency lies in the specified
                                           interval. [0.45, 0.55]
     -mc,--minCoverage <arg>               Minimum coverage to call a SNV
                                           [5.0]
     -mf,--minFrequency <arg>              Minimum (read) frequency to call a
                                           SNV [0.9]
     -mq,--minQual <arg>                   Minimum quality to call a SNV
                                           [30.0]
     -nt,--numThreads <arg>                Number of threads to use [1]
     -o,--output <arg>                     Path to the directory at which
                                           results files shall be generated.
     -r,--reference <arg>                  Path to a .fasta file; The
                                           reference genome sequence.
     -s,--sampleSpecification <arg>        Text file specifying the sample
                                           input files and sample
                                           meta-information. Each line has to
                                           contain at least the path to the
                                           samples .vcf file and may contain
                                           the following optional
                                           comma-separated fields in the given
                                           order:
                                           <NAME> (a string representing the
                                           name to use for the sample).
     -sDir,--sampleDirectory <arg>         Directory from which sample input
                                           files are collected. The specified
                                           directory and all its
                                           sub-directories are scanned for
                                           .vcf files which are collected as
                                           sample input files.
     -sf,--snpEff                          Run SNPEff. [false]
     -v,--verbose                          Verbose mode on, i.e. print
                                           stacktrace on errors. [false]

At least a reference genome (`-r`) and sample input (`-s` or `-sDir`) have to be specified.

---
## **Example Data**:
Example data to test the tool is available in the `./examples` directory. To run `MUSIAL` with the example data just execute the command stored in `/examples/treponema_pallidum_pallidum_bamA_Sun/command.txt` from the projects root directory.

---
## **Generated Output**:
- **runInfo.txt** stores information about the input parameters and files.
- **sequences.fasta** contains the full sequences of the reference and each sample with the variants – the most frequent allele, if multiples were detected called for the respective sample being included; note that the file does not represent a multiple sequence alignment.
- **variantContentTable.tsv** contains an overview of the per sample per variable position nucleotide content by extending the nucleotide alphabet by
the symbols ′.′ for a reference call or no information, ′−′ for a deletion with respect to the reference and ′∼′ for an insertion in any other sample.
- **variantAnnotationTable.tsv** yields the same structure as **variantContentTable.tsv** but instead of the nucleotide content it contains annotations of the called variant. Namely, the allele frequency, the call quality and depth of coverage as well as the SnpEff annotation, if present.
- **positionStatistics.tsv** and **sampleStatistics.tsv** yield total variant counts with respect to one variable position or sample, respectively.
- **snpEffSummary.tsv** contains, if enabled, each SnpEff annotated variant with a list of the respective samples in which the variant was detected.
- If enabled, the variantStructureAllocation directory provides the following content:
    - a copy of the PDB file used as input.
    - **residueDistanceMap.tsv** yielding the distance of each residue pair of the specified protein structure in Angstrom.
    - **variantStructureAllocation.json** yielding an allocation of the detected variants to the specified protein structure.

---
## **Remarks**:
- Although `MUSIAL` is capable of analyzing heterozygous variant calls, only the most frequent allele per sample per position is used for output file generation and optional structure assignment.
- **.vcf** files (see the official [VCF specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)) may contain ambiguous variant calls like the one in the example below:

      #CHROM POS ID REF ALT QUAL FILTER INFO

      A 4 . GCG G,GCGCG . PASS DP=100
    
  In such scenarios it is not clear how the reference content (GCG) is correctly aligned with the alternative content (GCGCG) and currently `MUSIAL` skips such content during the analysis.
