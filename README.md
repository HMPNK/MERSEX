# MERSEX
## SuperFast Kmer based / reference free sex marker search for wgs datasets
### Prerequisites
Install KMC, IDBA, GNU parallel, pigz and mawk
### Server requirements
currently parameters are set for a 96 CPU / 256+ GB RAM machine, processing 104 fish genomes (7-14Gbp reads per individual) had a peak mem usage of approximately 128 Gb and a runtime of ~4h.
### INSTALL BY CONDA/MAMBA
```sh
#create conda/mamba profile and install tools
mamba create -n MERSEX -c bioconda kmc pigz mawk parallel idba rapidgzip

#activate profile
mamba activate MERSEX

#get the repository
git clone https://github.com/HMPNK/MERSEX.git
cd MERSEX
cd src
bash COMPILE.sh #(you may have to install g++ and libraries first)
cd ..

#Copy or link your female data in directory "females" and your male data in directory "males"
#make sure that each individual has its ID separated with "_" in the filename
#and the filename contains only one "_": for example F1_1.fq.gz, F2_1.fq.gz, M1_xyz.fq.gz, M2_abc.fq.gz M2_def.fq.gz
#In this example data of the files "M2_abc.fq.gz" and "M2_def.fq.gz" will be merged while Kmer counting, because both are from sample "M2"
#also only use gzipped fastq (*.fq.gz suffix mandatory!)
#a minimum of 10 males and 10 females is recommended, more always makes sense... 

#run the script
bash MERSEX.sh

#PS: The core dump of the idba_ud assembler can be ignored.
This is expected behavior, given that we do not have paired end reads, which the assembler expects in the last step after contigs have already been build.
```
### Latest changes v0.3:
* streamlined pipeline, faster, less disk IO
* C++ tool kmc_merge for fast and memory efficient DB merging (up to hundreds)
* C++ tool fisher_mt for sex specific significance calculations and extraction of most significant kmers
* Assembly of significant kmers into contigs
