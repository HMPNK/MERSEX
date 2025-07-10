# MERSEX
## SuperFast Kmer based sex marker search for wgs datasets
### Prerequisites
Install KMC, GNU parallel, pigz and mawk
### INSTALL BY CONDA/MAMBA
```sh
#create conda/mamba profile and install tools
mamba create -n MERSEX -c bioconda kmc pigz mawk parallel

#activate profile
mamba activate MERSEX

#get the repository
git clone https://github.com/HMPNK/MERSEX.git
cd MERSEX

#Copy or link your female data in directory "females" and your male data in directory "males"
#make sure that each individual has its ID separated with "_" from in the filename
#and the filename contains only one "_": for example F1_1.fq.gz, F2_1.fq.gz, M1_xyz.fq.gz, M2_abc.fq.gz M2_def.fq.gz
#data of the files M2_abc.fq.gz M2_def.fq.gz will be merged while Kmer counting, because both are from sample "M2"
#also only use gzipped fastq (*.fq.gz suffix mandatory!) 

#run the script
bash MERSEX-v0.1.sh

```
