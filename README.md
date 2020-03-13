# MetaKraken

Snakemake pipeline for profiling composition of microbial communities from metagenomic shotgun sequencing data using MetaPhlAn2.

## Overview

### Input:

Quality processed paired-end fastq files from shotgun metagenome sequencing.

### Output:

Table of microbial species and their relative abundance for each sample, output/merged_abundance_table.txt
Heatmap of abundance results, output/abundance_heatmap_species.png


## Pipeline summary

Steps


1. Kraken2 is used to generate profiles of microbial clades and their abundances.
2. Bracken is used to re-estimate reads assigned to species.
3. Kreprot2mpa is used to generate profiles in mpa format which is similar to one used by Metaphlan
4. Normalized mpa files to generate percentages of each clade
5. Merge all sampel profiles into one table.
6. Generate a heatmap with 25 most abundant species.
7. Generate barplots with taxonomic information
 

Installation
To use this pipeline, navigate to your project directory and clone this repository into that directory using the following command:

git clone https://github.com/SycuroLab/metaphlan.git metaphlan
Note: you need to have conda and snakemake installed in order to run this. To install conda, see the instructions here.

To install snakemake using conda, run the following line:

conda install -c bioconda -c conda-forge snakemake
See the snakemake installation webpage for further details.

Config file
All the parameters required to run this pipeline are specified in a config file, written in yaml. See/modify the provided example file with your custom parameters, called config.yaml. This is the only file that should be modified before running the pipeline. Make sure to follow the syntax in the example file in terms of when to use quotations around parameters.

Data and list of files
Specify the full path to the directory that contains your data files in the config file. You also need to have a list of sample names which contains the names of the samples to run the pipeline on, one sample per line. You can run this pipeline on any number or subset of your samples. Sample names should include everything up to the R1/R2 (or 1/2) part of the file names of the raw fastq files. Specify the path and name of your list in the config file.

Description of parameters
Parameter	Description	Example
list_files	Full path and name of your sample list.	"/home/aschick/project/list_files.txt"
path	Location of input files.	"/home/aschick/project/data/filtered/"
for	Suffix of forward reads.	"_R1_filtered.fastq"
rev	Suffix of reverse reads.	"_R2_filtered.fastq"
Running the pipeline on Synergy

Test the pipeline by running snakemake -np. This command prints out the commands to be run without actually running them.

To run the pipeline on the Synergy compute cluster, enter the following command from the project directory:

snakemake --cluster-config cluster.json --cluster 'bsub -n {cluster.n} -R {cluster.resources} -W {cluster.walllim} -We {cluster.time} -M {cluster.maxmem} -oo {cluster.output} -e {cluster.error}' --jobs 500 --use-conda
The above command submits jobs to Synergy, one for each sample and step of the QC pipeline. Note: the file cluster.json contains the parameters for the LSF job submission system that Synergy uses. In most cases, this file should not be modified.

Results and log files
Snakemake will create a directory for the results of the pipeline as well as a directory for log files. Log files of each step of the pipeline will be written to the logs directory.

# Notes and Options

Information about  Kraken2 can be found on https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual
Information about  Bracken can be found on https://ccb.jhu.edu/software/bracken/

Databases:

Kraken:
At the moment standard Kraken2 db is downloaded and prepared to be used with this snakemake file. It is saved as kraken2db_2020_02_26.

The command used to build it is

kraken2-build --standard --threads 24 --db $DBNAME

Bracken:

Bracken database is built using Kraken2 db. It requres Kmer length (defualt=35) and read length of the paried-end reads. You might have to re-run bracke-build if your read length has not already been build


bracken-build -d ${KRAKEN_DB} -t ${THREADS} -k ${KMER_LEN} -l ${READ_LEN}



Options:

Following are the possible options for Kraken2 



Usage: kraken2 [options] <filename(s)>

Options:
  --db NAME               Name for Kraken 2 DB
                          (default: none)
  --threads NUM           Number of threads (default: 1)
  --quick                 Quick operation (use first hit or hits)
  --unclassified-out FILENAME
                          Print unclassified sequences to filename
  --classified-out FILENAME
                          Print classified sequences to filename
  --output FILENAME       Print output to filename (default: stdout); "-" will
                          suppress normal output
  --confidence FLOAT      Confidence score threshold (default: 0.0); must be
                          in [0, 1].
  --minimum-base-quality NUM
                          Minimum base quality used in classification (def: 0,
                          only effective with FASTQ input).
  --report FILENAME       Print a report with aggregrate counts/clade to file
  --use-mpa-style         With --report, format report output like Kraken 1's
                          kraken-mpa-report
  --report-zero-counts    With --report, report counts for ALL taxa, even if
                          counts are zero
  --memory-mapping        Avoids loading database into RAM
  --paired                The filenames provided have paired-end reads
  --use-names             Print scientific names instead of just taxids
  --gzip-compressed       Input files are compressed with gzip
  --bzip2-compressed      Input files are compressed with bzip2
  --help                  Print this message
  --version               Print version information







bracken-build -d kraken2db_2020_02_26/ -t 32 -k 35 -l 100

bracken is build at the moment with recommended kmer size of 35 and readlength of the test data was about 100

bracken -d kraken2db_2020_02_26 -i vaginal_sample_default -l S -o vaginal_sample_default.bracken


/export/home/hramay/miniconda3/envs/kraken_bracken_final/bin/bracken: illegal option -- h
Usage: bracken -d MY_DB -i INPUT -o OUTPUT -r READ_LEN -l LEVEL -t THRESHOLD
  MY_DB          location of Kraken database
  INPUT          Kraken REPORT file to use for abundance estimation
  OUTPUT         file name for Bracken default output
  READ_LEN       read length to get all classifications for (default: 100)
  LEVEL          level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
  THRESHOLD      number of reads required PRIOR to abundance estimation to perform reestimation (default: 0)


Kraken2 is working with nullarbor conda file but not with kraken2 conda file
