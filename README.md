<p align="center">
  <img width="256px" height="256px" src="https://github.com/Integrative-Transcriptomics/MUSIAL/blob/77ce6bca5c468d4a29f28a0b6c4c0f262026d1d4/media/logo_lucid_512.png">
</p>

# MUSIAL

**MUSIAL** (MUlti Sample varIant AnaLysis) is a Java command-line tool to analyze large sets of VCF files with prokaryotic single nucleotide variants (SNVs) and insertions/deletions (indels). It provides an interface for generating comprehensive statistics and alignments, as well as assessing variability at genome, gene and protein levels.

- **Integrates SnpEff and other Sequence Ontology compliant annotations** to help interpret variants.
- **Projection to genomic features (genes) facilitates allele- and proteoform-specific information** that supports the characterization of individual samples.
- **VCF based sequence reconstruction** at nucleotide and protein sequence level and tabular reports on sample, feature and variant statistics.

### 📖 Usage

An executable `jar` file (`Java 21`) is available from the [Releases](https://github.com/Integrative-Transcriptomics/MUSIAL/releases) section.
MUSIAL operates on a modular, task-based architecture that is primarily initiated by the `build` task, which creates a JSON file (_storage_) as its primary output; this is then used as input for all other tasks.

Details on the use of the software and tutorials can be found in the repository [Wiki](https://github.com/Integrative-Transcriptomics/MUSIAL/wiki). The general CLI usage is `java -jar MUSIAL-v2.4.2.jar <task>`, whereby the following tasks are available:

<details>
<summary><code>build</code> - Build a local database file (storage) in JSON format from variant calls; the mandatory input for other tasks.</summary>

```
Command line arguments of task build

 -C,--configuration <arg>   Path to a JSON file specifying the build task parameter configuration for MUSIAL. Visit the documentation for details.
```
</details>

<details>
<summary><code>expand</code> - Expand an existing storage file from variant call files and/or meta data.</summary>

```
Command line arguments of task expand

 -d,--dry-run          Only report on novel entries without writing the updated storage.
 -I,--storage <arg>    Path to a .json(.gz) file generated with the build task of MUSIAL.
 -m,--vcfMeta <arg>    Path to a .tsv or .csv file specifying sample annotations.
 -o,--output <arg>     Path to write the output file (default: overwrite input file).
 -V,--vcfFiles <arg>   List of file or directory paths. All files must be in VCF format.
```
</details>

<details>
<summary><code>view</code> - View the content (features, samples or variants; and their attributes) of a MUSIAL storage file.</summary>

```
Command line arguments of task view

 -C,--content <arg>   The content to view. One of FEATURES, SAMPLES, VARIANTS (case-insensitive).
 -I,--storage <arg>   Path to a .json(.gz) file generated with the build task of MUSIAL.
 -o,--output <arg>    Path to write the output file. If not provided, a default file will be created based on the input file (default). If `print` or
                      `stdout` is specified, the output will be printed to the console.
 -q,--query <arg>     One or multiple identifiers or genomic ranges (contig:start-end) to query.
```
</details>

<details>
<summary><code>profile</code> - Profile samples with respect to variants, alleles, or proteoforms.</summary>

```
Command line arguments of task profile

 -C,--content <arg>   The content to view. One of VARIANTS, ALLELES, PROTEOFORMS (case-insensitive).
 -I,--storage <arg>   Path to a .json(.gz) file generated with the build task of MUSIAL.
 -o,--output <arg>    Path to write the output file. If not provided, a default file will be created based on the input file (default). If `print` or
                      `stdout` is specified, the output will be printed to the console.
 -q,--query <arg>     One or multiple identifiers or genomic ranges (contig:start-end) to consider.
 -x,--reduced         Represent entries in a reduced format, i.e., sequence types as numbers with 0 as the reference or synonymous sequence and
                      variants without detailed call information.
```
</details>

<details>
<summary><code>sequence</code> - Generate and write sequence data.</summary>

```
Command line arguments of task sequence

 -a,--align             Whether to align sequences (optional, default: false).
 -c,--content <arg>     Whether to generate NUCLEOTIDE or AMINOACID sequences (optional, case-insensitive, default: NUCLEOTIDE).
 -f,--split <arg>       Whether to split output files by FEATURE, SAMPLE, BOTH, or NONE (optional, case-insensitive, default: FEATURE).
 -I,--storage <arg>     Path to a .json(.gz) file generated with the build task of MUSIAL.
 -l,--locations <arg>   One or multiple feature identifiers or genomic ranges (contig:start-end) to generate sequence data of. If none are provided,
                        all features or full contig ranges will be considered.
 -m,--merge             Whether to merge identical sequences (optional, default: false).
 -o,--output <arg>      Path to write the output. If not provided, the directory of the input storage is used. If a directory is provided, files are
                        created there. If a file is provided, its parent directory is used.
 -s,--samples <arg>     One or multiple sample identifiers to retrieve sequences for (optional).
 -v,--variable          Whether to only consider variable positions (optional, default: false).
```
</details>

### 🌐 Web Interface

MUSIAL is also available via a web interface at https://musial-tuevis.cs.uni-tuebingen.de/ currently running version `v2.3.10`.

### Build

MUSIAL `v2.4` is built with `JDK 21.0.6` and `Gradle 9.1.0`. If you want to compile the source code, run `gradle clean build` in the root directory of the project. The JavaDoc of the software is available at [https://integrative-transcriptomics.github.io/MUSIAL/javadoc/](https://integrative-transcriptomics.github.io/MUSIAL/javadoc/).

### Need Help?

- Detailed information about the software can be found in the repository's [Wiki](https://github.com/Integrative-Transcriptomics/MUSIAL/wiki). 
- Found an issue or have a feature request? Feel free to [Open a GitHub issue](https://github.com/Integrative-Transcriptomics/MUSIAL/issues/new).