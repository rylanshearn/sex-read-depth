# sex-read-depth
Takes paired end short illumina reads in fastq format from a male and female and computes ratio of coverage along a reference genome

## About
This is a collection of bash scripts that make use of `bowtie2`, `samtools`, and `awk` in order to map reads to a reference genome, and analyse the depth ratio along this genome. These tasks are performed by the following scripts:
### bowtieindex.sh
This takes a `.fasta` reference file and indexes it to enable bowtie2 to align reads to it.
### readscount.sh
This counts the number of reads for a sample (two `.fastq.gz` input files) and stores the number in an output file to be used later by `samtoolsdepth.sh`
### alignbowtie.sh
This aligns reads from two `fastq.gz` files onto a bowtie2 indexed reference genome
### samtoolssort.sh
This runs a few different commands that convert the `.sam` alignment file to `.bam`, then sorts the alignment file
### samtoolsdepth.sh
This calculates the read depth at each position along the reference genome. This is piped to awk, which corrects for the number of reads and calculates the coverage ratio between two samples on the fly
### samprofile.sh
generates polymorphism profile file from a .sam alignment
### profilesum.sh
This takes two snp profiles as output from sam2pro, and combines them in a single summary file
### *.pbs
Each of the above scripts have `.pbs` scripts that communicate with the Torque job submission system. This is the place to go if you want to change some of the resource allocation for each job)
### coveragemaster.pbs
This script launches all of the scripts for male:female read depth analyses through their equivalent `.pbs` scripts, specifying appropriate job dependencies and some resources that are more likely to need changing between datasets. This enables you to run analyses for male and female at the same time.
### snpdens.pbs
This script launches all of the scripts for male:female SNP diversity analyses through their equivalent `.pbs` scripts, starting from read alignment and counting. You can turn off the early stages of the pipeline if they have already been completed by coveragemaster.pbs, just follow the instructions in snpdens.pbs


## Requirements
This has only been tested with:
  * bowtie2 version 2.2.7
  * samtools version 0.1.19-96b5f2294a

Use these versions if you don't want any issues

## Installation
1. Clone the repository:
```
git clone https://github.com/rylanshearn/sex-read-depth
```
2. The path to the following must be added to `~/.profile`:
  * `/path/to/sex-coverage-cluster/`
  * `/path/to/bowtie2-2.2.7/`
  * `/path/to/samtools/`
  * `/path/to/sam2pro/`

## Usage
usage instructions are given in `coveragemaster.pbs`, and also in each `.sh` script

