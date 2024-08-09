<p align="center">
    <img src="media/icon.ico" width="128" alt="MUSIAL Icon"/>
</p>

## **MUlti Sample varIant AnaLysis**

MUSIAL computes prokaryotic genome, gene and protein sequence alignments from variant call datasets (using both SNVs and
InDels) derived from multiple samples of one species. It allows to assess the variability of a species on the genome,
gene, and protein sequence level, e.g., which positions in the species exhibit large or no variability, identify genes
with large variability, which samples agree on which proteoforms, etc. The project's architecture is organized in
tasks (you can find more information below). It provides multiple data export options for FASTA as well as tabular
formats.

---

## Dependencies and Building

The current version of MUSIAL `v2.3.8` is built with `JDK 17.0.8` and `Gradle 8.2.1`. A precompiled, executable **jar**
file is deposited at the **Releases** section of this repository. If you want to re-build MUSIAL,
type `gradle clean build` in the projects root directory. The **jar** file is then contained in the **/releases**
directory.

---

## Access via Web

MUSIAL is also accessible as a web-application at https://musial-tuevis.cs.uni-tuebingen.de/ (currently version `v2.8` is active), supporting not only the
computational capabilities of the software, but extending it with various user interactions and visualizations. The
source code of the web-application as well as a separate usage instruction will be available via the `web` branch of
this repository shortly.

---

## Workflow

We will soon add a description of MUSIAL's internal workflow.

---

## Usage

The compiled **jar** can be run from any command line terminal. MUSIAL implements distinct tasks; Thereby, the `build`
task is the starting point for all other tasks. The tasks currently available are listed below. Detailed usage
descriptions can be found in the collapsible sections.

```
 build            Build a local database/MUSIAL storage (opt. gzip compressed) JSON format file from variant calls.
 view_features    View annotated features from a MUSIAL storage file in a tabular format.
 view_samples     View annotated samples from a MUSIAL storage file in a tabular format.
 view_variants    View annotated variants from a MUSIAL storage file in a tabular format.
 export_table     Export variants from a MUSIAL storage file as a matrix-like .tsv file.
 export_sequence  Generate sequences in .fasta format from a MUSIAL storage file.
```

<details>
<summary><code>build</code></summary>

```
usage: java -jar MUSIAL-v2.3.3.jar build -C <arg> [-u] [-w <arg>]
command line arguments of task build
 -C,--configuration <arg>   Path to a .json file specifying the build configuration for MUSIAL. Please visit
                            https://github.com/Integrative-Transcriptomics/MUSIAL for a detailed explanation on how to specify the MUSIAL build
                            configuration file.
 -u,--uncompressed <arg>    Do not compress the storage file (Default: compressed).
 -w,--workdir <arg>         Path to a temporary working directory. By default './tmp/' is used.
```

The `build` task constitutes the first step for all analysis conducted with MUSIAL. The command uses the single
parameter `-c, --configuration` as input, pointing to the path of a **JSON** format build configuration file. The output
of the `build` task is a binary **brotli** compressed JSON file (`MUSIAL storage` file) that is used as the main input
for all other tasks.

#### Structure of the `build` Configuration File

The structure of the build configuration file follows a
strict [JSON schema](https://github.com/Integrative-Transcriptomics/MUSIAL/blob/v2.3/src/main/resources/musial_build_configuration.schema.json)
. The distinct properties and their meaning for running MUSIAL are described in a more human-readable form in the
following.

```
{
    "minimalCoverage": <Number>                 | Positive integer, variant call positions with read coverage below this value will be marked as ambiguous/rejected.
    "minimalFrequency": <Number>                | Float between 0.0 and 1.0, variant calls with read coverage relative to total position coverage below this value will be marked as ambiguous/rejected.
    "excludedPositions": {                      | Optional: Positions of contigs to exclude from the analysis.
        "<ContigName>": [<Number>, ... ]        | <ContigName> has to match the name of any sequence in 'referenceSequenceFile', <Number>s have to be any positions (1-based index) on that sequence; Unmatched entries are ignored.
    },
    "excludedVariants": {                       | Optional: Explicit variants to exclude from the analysis.
        "<ContigName>": [<Number>:<Var>, ... ]  | <ContigName> has to match the name of any sequence in 'referenceSequenceFile', <Number>s have to be any positions (1-based index) on that sequence; <Var> is interpreted as the explicit alternative (ALT) content of input .vcf files to ignore; Unmatched entries are ignored.
    },
    "referenceSequenceFile": "<FilePath>",      | Absolute or relative (to the working directory Java is run from) path to a .fasta|.fa|.fna file.
    "referenceFeaturesFile": "<FilePath>",      |                                                                   ... to a .gff3 file.
    "output": "<FilePath>",                     |                                                                   ... to store the output of the task. If the specified value does not end with .br, .br is appended at the end.
    "samples": {                                | Collection of samples, each sample is defined by one .vcf file.
        "<Name>": {                             | Any string value, used as internal name of the sample.
            "vcfFile": "<FilePath>",            | Absolute or relative (to the working directory Java is run from) path to a .vcf file.
            "annotations": {                    | Meta-information associated with the sample.
                "<Key>": "<Value>",             | <Key> and <Value> can be any strings.
                ...
            }
        },
        ...
    },
    "samplesDirectory": "<DirectoryPath>",      | Absolute or relative (to the working directory Java is run from) path to a directory. MUSIAL will collect all .vcf files in this directory as samples without annotations. The base name of the files are used as sample names.
    "features": {                               | Collection of features, each feature is considered an interval on any contig specified in 'referenceSequenceFile'.
        "<Name>": {                             | Any string value, used as internal name of the feature.
            "match_<AttributeKey>": "<Value>",  | <AttributeKey> has to match any attribute key in 'referenceFeaturesFile', <Value> is the value to match this feature for.
            "coding": true|false,               | Whether the feature is considered as a protein coding gene or not.
            "annotations": {                    | Meta-information associated with the feature.
                "<Key>": "<Value>",             | <Key> and <Value> can be any strings.
                ...
            }
        },
        ...
    }
}
```

#### Matching Features from a .gff File

The supposedly most complicated step is the definition of the features to be analyzed. This process can be explained
easily using an example. Imagine the following excerpt from a **.gff** file.

```
...
1		2	3	4	5	6	7	8	9
Contig1	Genbank	gene	239394	241367	.	+	.	ID=gene_0230;Name=priA;gbkey=Gene;gene=priA;gene_biotype=protein_coding;locus_tag=G_0230
Contig1	Genbank	gene	241468	242118	.	+	.	ID=gene_0231;Name=G_0231;gbkey=Gene;gene_biotype=protein_coding;locus_tag=G_0231
Contig1	Genbank	gene	242365	242910	.	+	.	ID=gene_0233;Name=G_0233;gbkey=Gene;gene_biotype=protein_coding;locus_tag=G_0233
Contig1	Genbank	gene	243045	243117	.	+	.	ID=gene_t0013;Name=G_t0013;gbkey=Gene;gene_biotype=tRNA;locus_tag=G_t0013
```

The 9th column represents the attribute column of the format. The corresponding entries in the build configuration file
for the first and last feature could look as follows:

```
"features": {
    "priA": {
        "match_Name": "priA",
        "coding": true,
        "annotations": {
            "biotype": "protein_coding"
        }
    },
    "G_t0013": {
        "match_locus_tag": "G_t0013",
        "coding": false,
        "annotations": {
            "biotype": "tRNA"
        }
    }
}
```

**Important:** If a genome-wide analysis is to be performed with MUSIAL, all contigs in the 'referenceSequenceFile' must be specified as separate features. If these are not already listed in the **.gff** file, they can be added manually by specifying the correct contig name, start and end attribute and an attribute to be matched:

```
1	2	3	4	5	6	7	8	9
Contig1	Custom	region	1	1139633	.	+	.	Name=genome
```

#### ❗ Input Restrictions

- MUSIAL requires the input reference sequence and all variant call format files to be indexed. If missing, the respective **.fai** and **.tbi** files are generated automatically.
- MUSIAL utilizes the `biojava GFF3Reader` to process **GFF** files:
  - The library is unable to parse files ending with .gff, so ensure that your GFF files use the .gff3 extension.
  - Please ensure that the GFF file does contain comment lines only at the start and no data other than the expected feature annotations are stored (many GFF files store sequence information in addition).
  - If contig names/**FASTA** headers are numbers, i.e., _>1_, _>2_, ... an index error will likely be thrown, as the value is interpreted as the index of the sequence in the 0-based index list of all sequences.
- Please ensure, that the contig names in the reference sequence, reference feature and variant call files match.
- Currently, only single sample **.vcf** files are supported, i.e., only one genotype per variant context.
- Complex InDel processing is made available by re-aligning the respective sequence content from the **.vcf** file, i.e., entries like `16333	ATTCA	GTTA` are split into `16333	A	G` and `16335	TC	T-`. **Note:** The outcome of this process may not be identical to the results of other alignment or mapping software, can lead to mixed substitution and InDel information and, thus, lead to somewhat ambiguous results. **We highly recommend to use variant call information without complex InDels**.

</details>

<details>
<summary><code>view_features</code></summary>

```
usage: java -jar MUSIAL-v2.3.3.jar view_features [-f <arg>] -I <arg> [-o <arg>]
command line arguments of task view_features
 -f,--features <arg>   Explicit space separated list of features to view (Default: all).
 -I,--storage <arg>    Path to a .json file generated with the build task of MUSIAL to view.
 -o,--output <arg>     Path to a file to write the output to (Default: stdout).
```

The output of the `view_features` task will look something like:

```
name	location	start	end	strand	number_of_alleles	number_of_proteoforms	number_of_substitutions	number_of_insertions	variable_positions	number_of_deletions	number_of_ambiguous	Annotation1	Annotation2
Gene1	Contig1		159684	160421	-	0			0			0			0			0			0			0			a		null
Gene2	Contig1		157943	159430	+	4			1			21			0			2.2177			12			4			b		1
```

#### Column Descriptions

- **name** The internal name of the feature.
- **location** The contig of the reference sequence this feature is located on.
- **start** The 1-based start position of the feature on the reference sequence in forward direction.
- **end** The 1-based end position of the feature on the reference sequence in forward direction.
- **strand** The strand (+/forward, -/reverse) of the feature.
- **number_of_alleles** Different nucleotide sequences of this feature across all samples.
- **number_of_proteoforms** Different aminoacid sequences of this (protein coding) feature across all samples.
- **number_of_substitutions** Nucleotide substitutions on this feature across all samples.
- **number_of_insertions** Nucleotide insertions (single base resolution) on this feature across all samples.
- **number_of_deletions** Nucleotide deletions (single base resolution) on this feature across all samples.
- **number_of_ambiguous** Ambiguous positions on this feature across all samples.
- **variable_positions** The percentage of variable nucleotide positions in percent relative to the feature length.
  Ambiguous calls are not counted.
- **Custom Annotations** The value for a user-defined annotation for this feature. All annotations of all viewed
  features are displayed as separate columns.
- ! All missing values are replaced with _null_.

</details>

<details>
<summary><code>view_samples</code></summary>

```
usage: java -jar MUSIAL-v2.3.3.jar view_samples -I <arg> [-o <arg>] [-s <arg>]
command line arguments of task view_samples
 -I,--storage <arg>   Path to a .json file generated with the BUILD task of MUSIAL to view.
 -o,--output <arg>    Path to a file to write the output to (Default: stdout).
 -s,--samples <arg>   Explicit space separated list of samples to view (Default: all).
```

The output of the `view_samples` task will look something like:

```
name	number_of_substitutions	number_of_insertions	number_of_deletions	number_of_ambiguous	allele_Gene1	proteoform_Gene1	allele_Gene2	proteoform_Gene2	Annotation1	Annotation2
Sample1	1			0			6			0			reference	reference		A1.s1.i0.d6.a0	P1.s56.i0.d4.a2.t0	a		1
Sample2	19			0			6			3			reference	reference		A3.s19.i0.d6.a3	P1.s56.i0.d4.a2.t0	a		null
```

#### Column Descriptions

- **name** The internal name of the sample.
- **number_of_substitutions** Nucleotide substitutions of this sample across all features.
- **number_of_insertions** Nucleotide insertions (single base resolution) of this sample across all features.
- **number_of_deletions** Nucleotide deletions (single base resolution) of this sample across all features.
- **number_of_ambiguous** Ambiguous positions of this sample across all features
- **allele\_[Feature]** The internal name of the assigned allele of this sample for each feature.
- **proteoform\_[Feature]** The internal name of the assigned proteoform of this sample for each feature.
- **Custom Annotations** The value for a user-defined annotation for this sample. All annotations of all viewed samples are displayed as separate columns.
- ! All missing values are replaced with _null_.

</details>

<details>
<summary><code>view_variants</code></summary>

```
usage: java -jar MUSIAL-v2.3.3.jar view_variants [-c <arg>] [-f <arg>] -I <arg> [-o <arg>] [-s <arg>]
command line arguments of task view_variants
 -c,--content <arg>    Sets the content type of viewed variants; One of `nucleotide` or `aminoacid` (Default: nucleotide).
 -f,--features <arg>   Explicit space separated list of features to restrict variants to (Default: all).
 -I,--storage <arg>    Path to a .json file generated with the BUILD task of MUSIAL to view.
 -o,--output <arg>     Path to a file to write the output to (Default: stdout).
 -s,--samples <arg>    Explicit space separated list of samples to restrict variants to (Default: all).
```

The output of the `view_variants -c nucleotide` task will look something like:

```
position	reference_content	alternate_content	feature		occurrence		type		frequency	snpeff_[ANN]
158915		C			T			Gene1		Sample1,Sample2		substitution	50		...
158916		A			C			Gene2		Sample1,Sample2		substitution	50		...
158930		CCTTCTT			C------			Gene2		Sample3			deletion	25		...
```

#### Column Descriptions

- **position** The position of the variant on the 1-based reference sequence (not relative to the feature).
- **reference_content** The reference sequence base content.
- **alternate_content** The alternative base content.
- **feature** The feature this variant is located on.
- **occurrence** Comma separated list of samples that yield this variant.
- **type** The type of the variant determined by MUSIAL (one of _substitution_, _insertion_, or _deletion_, with an optional \_ambiguous\_\_ prefix).
- **frequency** The frequency in percent of this variant relative to all samples.
- **snpeff\_[ANN]** MUSIAL conducts SnpEff annotation of all unambiguous nucleotide variant calls and extracts the added **ANN** fields as annotation fields. The full list of SnpEff **ANN** fields can be found [here](https://pcingola.github.io/SnpEff/snpeff/inputoutput/#ann-field-vcf-output-files).

#### Differences `view_variants -c aminoacid`

- The SnpEff annotation is omitted.
- The position is 1-based relative to the primary sequence of the protein.

</details>

<details>
<summary><code>export_table</code></summary>

```
usage: java -jar MUSIAL-v2.3.3.jar export_table [-c <arg>] -F <arg> -I <arg> -O <arg> [-s <arg>]
command line arguments of task export_table
 -c,--content <arg>   Sets the content type of viewed variants; One of `nucleotide` or `aminoacid` (Default: nucleotide).
 -F,--feature <arg>   The feature for which variants should be exported.
 -I,--storage <arg>   Path to a .json file generated with the BUILD task of MUSIAL to view.
 -O,--output <arg>    Path to a file to write the output to.
 -s,--samples <arg>   Explicit space separated list of samples to restrict variants to (Default: all).
```

The `export_table` task is used to create a complete overview of the variant calls of a subset (or all) analyzed samples with respect to a single feature. It can best be regarded as the short version of a multi-sample **.vcf** file. The output of the `export_table -c nucleotide` task will look something like:

```
Position	Reference	Sample1		Sample2		⋯
158193		G		.		1;72;G:0,A:72	⋯
158271		A		.		1;29;A:0,C:29	⋯
⋮		⋮		⋮		⋮		⋱
```

**Position** reflects the corresponding position in the reference sequence (1-based) and **Reference** the base in the reference sequence at this position. All subsequent columns describe the potential variant call at this position per sample or `.`, if no information about the position was present in the sample's **.vcf** file.
For `export_table -c nucleotide` each cell is described by `<CallIndex> ; <TotalCoverage> ; <A1>:<A1Coverage>,<A2>:<A2Coverage>,...`, where

- `<CallIndex>` The index of the called/most frequent allele with respect to the third field `<A1>:<A1Coverage>,...`.
- `<TotalCoverage>` The total read coverage at the position for the sample (extracted from the respective **.vcf** file).
- `<A1>:<A1Coverage>,<A2>:<A2Coverage>,...` A `,` separated list of alternative contents and their respective read support. The first entry at index 0 is the reference content.

#### Differences `export_table -c aminoacid`

- The position is 1-based relative to the primary sequence of the protein.
- Each cell is either the alternative content of the sample or `.`. This is due to the fact that only the most common allele per sample is used to derive information about the proteoform (see workflow).

</details>

<details>
<summary><code>export_sequence</code></summary>

```
usage: java -jar MUSIAL-v2.3.3.jar export_sequence [-a] [-c <arg>] -F <arg> -I <arg> [-k] [-m] -O <arg> [-r <arg>] [-s <arg>]
command line arguments of task export_sequence
 -a,--aligned           Whether to align sequences (Default: No alignment).
 -c,--content <arg>     Sets the content type of viewed variants; One of `nucleotide` or `aminoacid` (Default: nucleotide).
 -F,--feature <arg>     The feature for which variants should be exported.
 -I,--storage <arg>     Path to a .json file generated with the BUILD task of MUSIAL to view.
 -k,--conserved         Export conserved sites (Default: Only variantInformation sites).
 -m,--merge             Whether to merge samples by alleles and proteoforms (Default: No merging).
 -O,--output <arg>      Path to a file to write the output to.
 -r,--reference <arg>   Path to a .fasta file yielding the reference sequences with which the specified MUSIAL storage file was built. If the file is
                        not indexed, this wil be done automatically. This option is only required for `content=nucleotide` and `conserved`.
 -s,--samples <arg>     Explicit space separated list of samples to restrict variants to (Default: all).
```

The `export_sequence` task is used to create **FASTA** format sequence data from the variant calls of a subset (or all) analyzed samples with respect to a single feature.

</details>

<details>
<summary>The <code>MUSIAL storage</code> structure</summary>

#### This is not essential for using MUSIAL!

We will soon add a description of MUSIAL's internal storage structure.

</details>

---

## Open Developments and Limitations

#### Processing Efficiency

- When developing MUSIAL, we took care to ensure an efficient computing and storage structure. However, we do not use dedicated index files, which has some limitations:
  - We refrain from using large (eukaryotic) data sets (at least the possibility of processing them has not yet been validated).

#### Other OMICS Integration

- MUSIAL allows to integrate nucleotide variants to the protein level (i.e., the inference of SAVs, aminoacid InDels and proteoforms). A corresponding integration into the RNA level is currently not possible.
- We are looking for a method to predict the effects of SAVs and amino acid InDels on the integrity and function of proteins.

#### Upcoming Features

- Provision and description of an example data set.
- Benchmarking.

---

## Contact

If you encounter any problems when using MUSIAL, please feel free to open an issue or contact me via `simon.hackl@uni-tuebingen.de`
