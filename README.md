# **MUSIAL - MUlti Sample varIant AnaLysis**

`Contact:   simon.hackl@uni-tuebingen.de`

---

## **Description**:

![MUSIAL Logo](media/logo_inverse.png)

MUSIAL computes prokaryotic genome, gene and protein sequence alignments from variant call datasets (using both SNVs and Indels) derived from multiple samples of one species. It allows to assess the variability of a species on the genome, gene, and protein sequence level, e.g., which positions in the species exhibit large or no variability, identify genes with large variability, which samples agree on which proteoforms, etc. The project's architecture is organized in modules (you can find more information below). It provides multiple data export options for FASTA as well as tabular formats.

---

## **Dependencies and Building**:

- JDK 15+
- Gradle 7+

A precompiled **.jar** can be found at **Releases**. MUSIAL program can be built with gradle (https://gradle.org). For that just type `gradle clean build` in the projects root directory. The **.jar** file is then contained in the `/releases` directory.

---

## **Access via Web**:

As an alternative to using MUSIAL locally, you can access the software via a webserver at https://tuevis.cs.uni-tuebingen.de/. The source code of the web-plattform is available at https://github.com/Integrative-Transcriptomics/MUSIAL-WEB.

---

## **Usage**:

The compiled **.jar** can be run from any command line terminal. MUSIAL implements distinct modules to run tasks. The common input for each pipeline is a `JSON` file which describes which modules should be executed with which parameters. The exact structure of the MUSIAL configuration files is defined in the [MUSIAL configuration JSON schema](https://github.com/Integrative-Transcriptomics/MUSIAL/blob/v2.1/MUSIAL_CONFIGURATION.schema.json). All modules except the `BUILD` module use a database like local `JSON` file (the output of the `BUILD` module) as input. The structure of this file is defined in the [MUSIAL output JSON schema](https://github.com/Integrative-Transcriptomics/MUSIAL/blob/v2.1/MUSIAL_BUILD_OUTPUT.schema.json). Prettified versions of the schemas can be accessed from the `./musial_*_schema_pretty/` directories in this repository.

Currently available modules are:

```
MUSIAL v2.1 | License: GPL-3.0 License | Contact: simon.hackl@uni-tuebingen.de
usage: java -jar MUSIAL-v2.1.jar -c <arg> [-k] [-s]
command line arguments:

-c,--configuration <arg> Path to a .json file specifying the run configurations for one or more MUSIAL modules. Please visit https://github.com/Integrative-Transcriptomics/MUSIAL for a detailed explanation on how to specify MUSIAL configuration files.

Available modules are:
BUILD: Generate a new variants dictionary .json file from multiple .vcf files wrt. one reference (.fasta + .gff).
EXTRACT: Extract stored nucleotide or aminoacid variants of specified samples and features as a .tsv (tabular overview) or .fasta (sequences) file.

-k,--compress If set, any output files will be compressed.

-s,--silent If set, no console output will be generated.
```
