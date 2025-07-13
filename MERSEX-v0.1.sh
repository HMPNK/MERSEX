##A PIPELINE FOR FINDING SEX-specific kmers ##
#NOT for POOLSEX data! Needs sequencing from individuals!
#needs at least 7-10x whole genome coverage per sample (not tested on lower so far)
#read files MUST be sorted by sex and go to folder "males" and folder "females"
#read files naming scheme MUST BE <sample-id>_1.fq.gz and <sample-id>_2.fq.gz

#creating list of readfiles per sample id
find males/| grep fq.gz$ | sort -V | awk '{split($1,a,/[/_]/);;print $1 > a[2]".list"}'
find females/| grep fq.gz$ | sort -V | awk '{split($1,a,/[/_]/);;print $1 > a[2]".list"}'

#CREATING KMC RUNs per sample
ls *.list | awk '{print "kmc -t20 -k27 -r @"$1" "$1".res ./kmc_tmp_dir/ > "$1".res.log 2>&1"}' > KMC_BATCH.sh

#executing KMC runs
mkdir ./kmc_tmp_dir
#not more than 4 parallel jobs! IO is limiting! Instead more kmc threads (i.e. -t20 )
cat KMC_BATCH.sh | nohup parallel -j 4

#create operations file for union of all kmers
find males/ females/ | grep fq.gz | sort -V | cut -f 1 -d '_' | uniq | awk 'BEGIN{print "INPUT:";}{n=split($1,a,"/");if(a[1]=="males"){i++;print "m"i" = " a[2]".list.res";text=text"m"i" + "} else if(a[1]=="females"){j++;print "f"j" = " a[2]".list.res";text=text"f"j" + "}} END{print "OUTPUT:\nallunion = "substr(text,1,length(text)-2)"\nOUTPUT_PARAMS:\n-cs1000000"}' > union-op-def.file

#execute union (may use lots of memory here)
kmc_tools complex union-op-def.file

#remove counts from union
kmc_tools transform allunion compact allunion-nocounts

#replace counts of all kmers with those of single samples (no match kmers will have count "1" !!!)
#and dump all counts to gzipped files (8 threads in parallel work!)
#the combination of these steps saves disk space (increasing parallel threads will use more disk for temporary kmc databases!)
ls *list | awk '{n=$1;gsub(".list","",n);print "kmc_tools simple allunion-nocounts -ci0 "$1".res -ci0 union all_"n" -ocright; kmc_dump all_"n" /dev/stdout | cut -f 2 | pigz -c > all_"n".counts.gz; rm all_"n"*kmc_*"}' | nohup parallel -j 8

#remove old KMC-DBs
rm *res* -f

#dump all kmers for once
kmc_dump allunion /dev/stdout | cut -f 1 | pigz -c > all.kmers.gz

#free disk
rm -f *.kmc_suf *.kmc_pre

#dump to table
#make fifos
awk '{if($1!="INPUT:" && $1!="OUTPUT:"){print};if($1=="OUTPUT:"){exit;}}' union-op-def.file | awk 'BEGIN{printf "mkfifo "} {printf $1" "}' | bash
#gunzip to fifo's
awk '{if($1!="INPUT:" && $1!="OUTPUT:"){print};if($1=="OUTPUT:"){exit;}}' union-op-def.file | awk '{gsub(".list.res",".counts.gz",$3);print "gunzip -c all_"$3" > "$1" &"}' | bash
mkfifo start
gunzip -c all.kmers.gz > start &

#final table merge, Kmer-positive female/male counting and compress
export FCOUNT="$(awk '{if($1!="INPUT:" && $1!="OUTPUT:"){print};if($1=="OUTPUT:"){exit;}}' union-op-def.file | grep ^f | wc -l)"
awk '{if($1!="INPUT:" && $1!="OUTPUT:"){print};if($1=="OUTPUT:"){exit;}}' union-op-def.file | awk 'BEGIN{text="start";} {text=text" "$1;} END{print "paste "text; }' | bash | parallel --keep-order -j 16 -l 1000000 --pipe 'mawk -v cmin=2 -v fcount=$FCOUNT -f count.samples.awk - ' | pigz -c > FINAL-table.tsv.gz

#clean-up
#remove FIFOs
rm -f start f? f?? f??? m? m?? m???
#remove lists and related
rm -f *list* nohup.out
#remove tmp-dir
rm kmc_tmp_dir/ -rf
#remove temporary count-files
rm -f *.counts.gz all.kmers.gz

