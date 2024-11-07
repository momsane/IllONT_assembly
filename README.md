# Metabee

This repository contains all the scripts and conda environments to assemble high quality genomes from Oxford Nanopore and Illumina reads.

## Installation

Requirements:
- conda v23.10.0
- mamba v1.5.6\
The following tools installed in their own conda environment:
- snakemake v7.32.4
- R v.4.3.2
  - r-ggplot2 v3.5.0
  - r-tidyverse v2.0.0
- seqkit v.2.6.1
- fastqc v0.11.8
- multiqc v1.6 
- bbmap v39.01
- nanoq v0.10.0
- filtlong v0.2.1
- hybracter v.0.9.0
- bandage v0.8.1 

- checkm v1.2.2 with the database downloaded, decompressed and dearchived. The path to this database must be put in the environment yaml file at the end, so that it is stored in an environment variable.
- gtdb-tk v2.4.0 with the database downloaded and decompressed. The path to this database must be put in the environment yaml file at the end, so that it is stored in an environment variable.

All environment YAML files can be found in /envs.

## Usage

Note: the resources directive in each rule of the snakefile is written for execution on the slurm cluster of UNIL.

1. Create a directory with the following sub-directories:

```
.
│
└─── benchmarks
│   
└─── data
│
└─── logs
│   
└─── results
│
└─── workflow
```

2. Clone this repository in /workflow.
3. Download and uncompress all required databases (CheckM, GTDB-Tk, DRAM).
4. Set up the metadata.tsv file in /config with the following columns:
- sample
- Illumina_avail: 1 if Illumina reads are available, 0 if not
- estimated_genome_size in bp
5. Upload raw Illumina and ONT reads to /data and add the paths to the snakefile.
6. Activate your snakemake conda environment and run the snakefile.

If snakemake raises issues with conflicting versions when creating the conda environments, it helps to (temporarily) set the channel priority to flexible.

Check points requiring manual input:
- after running fastqc, add column 'Illumina_adapter' to metadata.tsv file indicating which adapter was found in each sample
- if hybracter cannot assemble contigs, it will exit without producing an output, which will cause snakemake to consider it as failed. You need to check the log file.
- if hybracter cannot assemble a genome of size somewhat close to the estimated_genome_size, it will flag the output as incomplete.
- after running hybracter, add column 'genome_assembled' to metadata.tsv file indicating if hybracter could assemble a genome (1) or not (0) 
- in rule all, comment the line requesting the hybracter output for all samples, otherwise snakemake will keep trying to run hybracter on the samples for which contigs cannot be assembled
- after running the genomes QC until CheckM, add column 'include_downstream' to indicate if you wish to pursue the analysis of the genome (1) or not (0). If the genome is less than 95% complete (CheckM) it might not be worth it.

## Contributions

