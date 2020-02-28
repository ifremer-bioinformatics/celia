[![Licence](https://img.shields.io/badge/licence-Affero_GPL_3.0-orange.svg)]()
[![Version](https://img.shields.io/badge/version-beta-red.svg)]()
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.07.0-blue.svg)](https://www.nextflow.io/)
[![Install](https://img.shields.io/badge/install-SeBiMER_gitlab-brightgreen.svg)](https://gitlab.ifremer.fr/bioinfo/CELIA)

# Introduction

CELIA (automatiC gEnome assembLy marIne prokAryotes) is a bioinformatics workflow used to automate the genome assemblies of prokaryotes from Illumina data

For now, the workflow allows:
- to make a quality control of the data using FastQC and MultiQC
- to assemble the genome using Unicycler
- to remove contaminants (vectors, adpatators or PhiX) using homemade script and NCBI UniVec database
- to evaluate the assembly quality using BUSCO, mapping coverage (Bowtie2 + MoseDepth)
- to compute ANI scores using dedicated database

## Quick start

1. Clone the current gitlab repertory in your working directory on DATARMOR

```
git clone https://gitlab.ifremer.fr/bioinfo/celia.git
```

2. Once on DATARMOR, edit the celia/config/params.config file with your analysis parameters

3. Add a directory with your data (paired fastq.gz files)

4. Run the analysis

```
./RunCELIA.sh
```

## Workflow process

<img width="700" height="636" src="./CELIA_workflow.jpg">

## Further development

- Identify and isolate assembled plasmids
- Structural and functional annotation (Prokka)
- Generate HTML report

## License and Credits
MAEVA is released under the GNU Affero General Public License, Version 3.0. AGPL

It is developped by Alexandre Cormier, bioinformatics engineer at the bioinformatics service of IFREMER (SeBiMER).

-- (c) 2020 - SeBiMER, Ifremer
