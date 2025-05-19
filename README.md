<p align="center">
  <img width="256px" height="256px" src="https://github.com/Integrative-Transcriptomics/MUSIAL/blob/v2.4/media/logo_lucid_512.png">
</p>

# MUSIAL

**MUSIAL** (MUlti Sample varIant AnaLysis) is a Java command-line tool designed to analyze and summarize single nucleotide variants (SNVs) and insertions/deletions (indels) across multiple prokaryotic samples.
The software aggregates and analyzes variant calls from multiple samples of a prokaryotic species and provides an interface to generate comprehensive statistics and alignments at the genome, gene and protein level.
MUSIAL enables a comprehensive assessment of variability within a species at the genome, gene and protein level, providing insights into, for example, conserved and variable regions, diversity at the gene level and common proteoforms among samples.

## ✨ Features

- **Integrates SnpEff and other Sequence Ontology compliant annotations** to help interpret variants.
- **Projection to genomic features (genes) facilitates allele- and proteoform-specific information** that supports the characterization of individual samples.
- **VCF based sequence reconstruction** at nucleotide and protein sequence level and tabular reports on sample, feature and variant statistics.

## 🔧 Usage

An executable `jar` file (`Java 21`) is available from the [Releases](https://github.com/Integrative-Transcriptomics/MUSIAL/releases) section.
MUSIAL operates on a modular, task-based architecture that is primarily initiated by the `build` task, which creates a JSON file (_storage_) as its primary output; this is then used as input for all other tasks.

The general CLI usage is `java -jar MUSIAL-v2.4.0.jar <task>`, whereby the following tasks are available:

<details>
<summary><code>build</code> - Build a local database file (storage) in JSON format from variant calls; the mandatory input for other tasks.</summary>

```
Command line arguments of task build

 -C,--configuration <arg>   Path to a JSON file specifying the build task parameter configuration for MUSIAL.
```
</details>

<details>
<summary><code>expand</code> - Expand an existing storage file from variant call format (VCF) files.</summary>

```
Command line arguments of task expand

 -I,--storage <arg>    Path to a .json(.gz) file generated with the build task of MUSIAL.
 -m,--vcfMeta <arg>    Path to a .tsv or .csv file specifying sample annotations.
 -o,--output <arg>     Path to write the output file (default: overwrite input file).
 -p,--preview          Only report on novel entries without writing the updated storage.
 -V,--vcfInput <arg>   List of file or directory paths. All files must be in VCF format.
```
</details>

<details>
<summary><code>view</code> - View the content (features, samples or variants) and their attributes, of a MUSIAL storage file.</summary>

```
Command line arguments of task view

 -C,--content <arg>   One of sample, allele, call, variant, type, feature.
 -f,--filter <arg>    List of feature-, sample names, and/or positions for which the output is to be filtered (default: no filters). Entries may be
                      ignored depending on the content.
 -I,--storage <arg>   Path to a .json(.gz) file generated with the build task of MUSIAL.
 -o,--output <arg>    Path to directory or file to write the output to (default: stdout).
```
</details>

<details>
<summary><code>sequence</code> - Export FASTA format sequences of features from a MUSIAL storage file.</summary>

```
Command line arguments of task sequence

 -c,--content <arg>    One of `nt` or `aa` (default: `nt`).
 -F,--features <arg>   List of feature names to export data for. Non-coding features are skipped if `content` is `aa`.
 -I,--input <arg>      Path to a .json(.gz) file generated with the build task of MUSIAL.
 -k,--conserved        Export conserved sites.
 -m,--merge            Export sequences per allele or proteoform instead of per sample.
 -o,--output <arg>     Path to a directory to write the output files to (default: parent of input).
 -r,--reference        Include the reference sequence within the export.
 -s,--samples <arg>    List of sample names to restrict the sequence export to.
 -x,--strip            Strip all gap characters from the exported sequences.
```
</details>

---

Further details on the use of the software and internal workflows can be found in the repository [Wiki](https://github.com/Integrative-Transcriptomics/MUSIAL/wiki).

## 🌐 Web Interface

To provide user-friendly access to its functionalities, MUSIAL is available via a web interface at https://musial-tuevis.cs.uni-tuebingen.de/; currently running version `v2.3.8`.

## 🔨 Build

MUSIAL `v2.4` is built with `JDK 21.0.6` and `Gradle 8.2.1`. If you want to compile the source code, run `gradle clean build` in the root directory of the project. The JavaDoc of the software is available at [https://integrative-transcriptomics.github.io/MUSIAL/javadoc/](https://integrative-transcriptomics.github.io/MUSIAL/javadoc/).

## 🙋 Need Help?

- 🎓 Detailed information about the software can be found in the repository's [Wiki](https://github.com/Integrative-Transcriptomics/MUSIAL/wiki). 
- 🐛 Found an issue or have a feature request? Feel free to [Open a GitHub issue](https://github.com/Integrative-Transcriptomics/MUSIAL/issues/new).