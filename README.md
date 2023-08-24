# **MUSIAL - MUlti Sample varIant AnaLysis**

`Contact:   simon.hackl@uni-tuebingen.de`

---

## **Description**:

![MUSIAL Logo](media/logo.png)

MUSIAL computes prokaryotic genome, gene and protein sequence alignments from variant call datasets (using both SNVs and Indels) derived from multiple samples of one species.
It allows to assess the variability of a species on the genome, gene, and protein sequence level, e.g., which positions in the species exhibit large or no variability, identify genes with large variability, which samples agree on which proteoforms, etc.
The project's architecture is organized in tasks (you can find more information below).
It provides multiple data export options for FASTA as well as tabular formats.

---

## **Dependencies and Building**:

- JDK 17.0.8+
- Gradle 8.2.1

A precompiled **.jar** can be found at **Releases**.
The MUSIAL program can be built with gradle (https://gradle.org).
For that just type `gradle clean build` in the projects root directory.
The **.jar** file is then contained in the `/releases` directory.

---

## **Access via Web**:

As an alternative to using MUSIAL locally, you can access the software via a webserver at https://tuevis.cs.uni-tuebingen.de/.
The source code of the web-platform is available at https://github.com/Integrative-Transcriptomics/MUSIAL-WEB.

---

## **Usage**:

The compiled **.jar** can be run from any command line terminal.
MUSIAL implements distinct tasks; thereby, the first step of the pipeline is constituted by the build task, constructing a local JSON database of variant call inforamtion used by the other tasks.
The build task uses as special `JSON` format configuration file as input.
A schema of this file is defined at [MUSIAL build configuration JSON schema](https://github.com/Integrative-Transcriptomics/MUSIAL/blob/v2.2/musial_build_configuration_schema/schema_doc.html).

Currently available tasks are:

```
[MUSIAL v2.2 | License: GPL-3.0 License | Contact: simon.hackl@uni-tuebingen.de]
usage: java -jar MUSIAL-v2.2.jar
MUSIAL provides distinct tasks to computes prokaryotic genome, gene and protein sequence alignments from variant call datasets from multiple samples
of one species. Available tasks are:

 build            Build a local .json database (MUSIAL storage) from variant calls.
 view_features    View annotated features from a built MUSIAL storage in a tabular format.
 view_samples     View annotated samples from a built MUSIAL storage in a tabular format.
 view_variants    View annotated variants from a built MUSIAL storage in a tabular format.
 export_table     Export variants from a built MUSIAL storage into a matrix-like .tsv file.
 export_sequence  Generate sequences in .fasta format from a built MUSIAL storage.


Call `java -jar MUSIAL-v2.2.jar <task> [-h|--help]` for more information.
```

We have added a (zip compressed) example dataset based on _M. tuberculosis_ samples you can run MUSIAL with at the examples directory. The origin and output as well as the configuration specification of the example dataset are described inside the respective README file.