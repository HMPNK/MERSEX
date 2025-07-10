# MERSEX
## SuperFast Kmer based sex marker search for wgs datasets
### Prerequisites
Install KMC, GNU parallel, pigz and mawk
### INSTALL BY CONDA/MAMBA
```sh
mamba create -n MERSEX -c bioconda kmc pigz mawk parallel
#get the repository
git clone https://github.com/HMPNK/MERSEX.git
#Copy or link your female data in directory "females" and your male data in directory "males"
#make sure that each individual has its ID separated with "_" from in the filename
#and the filename contains only one "_": for example F1_1.fq.gz, F2_1.fq.gz
#also only use gzipped fastq (*.fq.gz suffix mandatory!) 

#run the script
bash MERSEX-v0.1.sh

```
