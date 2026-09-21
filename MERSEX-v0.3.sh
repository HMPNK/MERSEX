#!/usr/bin/env bash
set -e
set -o pipefail
set -u

##A PIPELINE FOR FINDING SEX-specific kmers ##
#NOT for POOLSEX data! Needs wgs sequencing from individuals!
#needs at least 6-10x whole genome coverage per sample
#Set LOWCOV=1 for datasets with less than 6x whole genome coverage
#read files MUST be sorted by sex and go to folder "males" and folder "females"
#read files naming scheme MUST BE <sample-id>_1.fq.gz and <sample-id>_2.fq.gz !!!

#settings (currently optimized for server with =128 GB RAM, >=1TB FREE DISK, >=12 CPU THREADS)
export LOWCOV=0      #Set to 1 for low sequencing coverage data (1-6x) per individual (genome skimming approaches), may increase runtime and disk usage!

export KMCTHREADS=6  #kmc kmer counting of fastq files, number of threads
export KMCJOBS=2     #kmc kmer counting of fastq files, number of kmc jobs running in parallel
                     #Do not use more than 4 parallel jobs! IO is limiting! Instead use more kmc threads!
export KMCMEM=60     #max memory per KMC job

export MERGE=6       #reading threads for kmc_merge (the number of true CPU cores is optimal)

export KPOS=1        #Kmer counts equal or larger are treated as valid match (kmer present in sample) in fisher test
export KNEG=0        #Kmer counts equal or lower are treated as no match (kmer absent in sample) in fisher test
export PVAL=0.01     #P-value cut off, values equal or smaller than this go to signifcant sex differences table
export PJOBS=6       #number of threads for creating significant sex differences table

#creating list of readfiles per sample id (for low coverage data we just use each file twice to meet the minimum kmer criteria of 2!)
if [ $LOWCOV -eq 0 ]
then
        echo "USING DEFAULT MODE"
        find males/| grep fq.gz$ | sort -V | awk '{split($1,a,/[/_]/);print $1 > "male-"a[2]".list"}'
        find females/| grep fq.gz$ | sort -V | awk '{split($1,a,/[/_]/);print $1 > "female-"a[2]".list"}'
else
        echo "USING LOW COVERAGE MODE"
        find males/| grep fq.gz$ | sort -V | awk '{split($1,a,/[/_]/);print $1"\n"$1 > "male-"a[2]".list"}'
        find females/| grep fq.gz$ | sort -V | awk '{split($1,a,/[/_]/);print $1"\n"$1 > "female-"a[2]".list"}'
fi

#CREATING KMC RUNs per sample
ls *.list | awk -v cpu=$KMCTHREADS -v mem=$KMCMEM '{gsub(".list","");print "kmc -sm"mem" -t"cpu" -k27 -r @"$1".list "$1" ./kmc_tmp_dir/ > "$1".kmc.log 2>&1"}' > KMC_BATCH.sh

#executing KMC runs
mkdir ./kmc_tmp_dir
#not more than 4 parallel jobs! IO is limiting! Instead more kmc threads (i.e. -t20 )
cat KMC_BATCH.sh | nohup parallel -j $KMCJOBS

#dump Dbs to single table

export KMCDBS="$(find | cut -f 2 -d'/' | grep kmc_pre$ | sed "s/.kmc_pre//g" |sort -V | awk '{printf $1" "}' )"
export FCOUNT="$(find | cut -f 2 -d'/' | grep kmc_pre$ | grep ^female | wc -l )"
export MCOUNT="$(find | cut -f 2 -d'/' | grep kmc_pre$ | grep ^male | wc -l )"

echo "Merging KMC databases into table:"
echo $KMCDBS
echo "Number of female samples: "$FCOUNT
echo "Number of   male samples: "$MCOUNT

./kmc_merge -b 131072 -S 131072 -L 4097 -t $MERGE -z -o KMCMERGED.tsv.gz $KMCDBS

#Extract significant sex differentiating kmers:
pigz -dc KMCMERGED.tsv.gz | ./fisher_mt -f $FCOUNT -H -p $PVAL -t $PJOBS --cmin1 $KPOS --cmin2 $KNEG | pigz -c > KMCMERGED_pval$PVAL.tsv.gz
