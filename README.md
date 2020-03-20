# MetaKraken

Snakemake pipeline for profiling composition of microbial communities from metagenomic shotgun sequencing data using Kraken2 and Bracken.

## Overview

Input:

Quality processed paired-end fastq files from shotgun metagenome sequencing.

Output:

* Table of microbial species and their relative abundance for each sample, output/merged_abundance_table.txt
* Heatmap of abundance results, output/abundance_heatmap_species.png
* Profile plot pdfs generated per sample or a user defined metadata varaible

## Pipeline summary

### Steps


1. Kraken2 is used to generate profiles of microbial clades and their abundances.
2. Bracken is used to re-estimate reads assigned to species.
3. Kreprot2mpa is used to generate profiles in mpa format which is similar to one used by MetaphlAn2
4. Normalized mpa files to generate percentages of each clade
5. Merge all sampel profiles into one table.
6. Generate a heatmap with N most abundant species.
7. Generate barplots with taxonomic information
 

## Installation
To use this pipeline, navigate to your project directory and clone this repository into that directory using the following command:

```
git clone https://github.com/SycuroLab/metakraken.git metakraken
```

Note: you need to have **conda** and **snakemake** installed in order to run this. To install conda, see the instructions [here](https://github.com/ucvm/synergy/wiki). 

To install snakemake using conda, run the following line:

```
conda install -c bioconda -c conda-forge snakemake
```

See the snakemake installation [webpage](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for further details. 

## Config file

All the parameters required to run this pipeline are specified in a config file, written in yaml. See/modify the provided example file with your custom parameters, called config.yaml. 

### Important parameters

#### General 

1. list_files: a file containing the sample names
2. path: path to trimmed fastq files
3. for: Suffix for forward reads
4. rev:  suffix for reverse reads

#### Kraken + Bracken related

1. kraken_db: Location of the Kraken database to use.
2. level: Level for bracken taxa (defualt 'S', Option:'D', 'P', 'C', 'O', 'F', 'G', 'S')
3. redreadlen: Read length required for Bracken ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)
    

4. threshold: specifies the minimum number of reads required for a classification at the specified rank

#### Barplots related (used by plot_profile.R)

1. variableX:  X-axis variable to make barplots at different taxa level
2. variableFacet: Facet variable to make barplots at different taxa level
3. topN: Top N most abundant families/genra/species to be plotted by the R script
4. metadata: comma delimited (csv) metadata file with Sample Names same as used in list_files file.

## Running Instructions

Test the pipeline by running snakemake -np. This command prints out the commands to be run without actually running them.

To run the pipeline on the Synergy compute cluster, enter the following command from the project directory:

snakemake --cluster-config cluster.json --cluster 'bsub -n {cluster.n} -R {cluster.resources} -W {cluster.walllim} -We {cluster.time} -M {cluster.maxmem} -oo {cluster.output} -e {cluster.error}' --jobs 500 --use-conda
The above command submits jobs to Synergy, one for each sample and step of the QC pipeline. Note: the file cluster.json contains the parameters for the LSF job submission system that Synergy uses. In most cases, this file should not be modified.

Results and log files
Snakemake will create a directory for the results of the pipeline as well as a directory for log files. Log files of each step of the pipeline will be written to the logs directory.

# Addition information 

Information about  Kraken2 can be found on https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual
Information about  Bracken can be found on https://ccb.jhu.edu/software/bracken/

## How to build Databases:

### Kraken:
At the moment standard Kraken2 db is downloaded and prepared to be used with this snakemake file. It is saved as kraken2db_2020_02_26 in  /gpfs/snyder_work/shared/lsycuro_labshare/dbs/kraken2db_2020_02_26

The command used to build it is

```
kraken2-build --standard --threads THREADS --db DBNAME
```

You can also build custom databases. Please read the kraken2 manual for more information.


### Bracken:

Bracken files are built using Kraken2 db. It requres Kmer length (defualt=35) and read length of the paried-end reads. You might have to re-run bracke-build if your read length has not already been build

```
bracken-build -d KRAKEN_DB -t THREADS -k KMER_LEN -l READ_LEN
```


Note: Kraken2 is working with nullarbor conda file but not with kraken2 conda file
