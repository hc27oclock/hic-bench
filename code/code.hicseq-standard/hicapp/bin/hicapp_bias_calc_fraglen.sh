####
# calculate fragment length bias for Hi-C data
####

#### set path of the pipeline
bindir=`dirname $0`


#### get pars
rebed=$1 #rebed=/michorlab/hjwu/code/hic_pipeline/ref/hg19/Digest_hg19_MboI_None.clear
chrom=$2 #chrom=/michorlab/hjwu/code/hic_pipeline/ref/hg19/hg19.genome
res=$3 #res=1000000
tag=$4 #tag=1m


#### calculate fragment based bias score
if [ ! -e $rebed.frag.fl.sortBed ]
then
	#### get fragment length
	awk '{print $1"\t"$2"\t"$3"\t"$1"_"int(($2+$3)/2)"\t"$3-$2}' $rebed |sortBed > $rebed.frag.fl

	#### sort files
	bedtools sort -i $rebed.frag.fl > $rebed.frag.fl.sortBed

	#### clear intermediate files
	rm $rebed.frag.fl
else
	echo "#*_*# <$rebed.frag.fl.sortBed> already exists...Skip~~~"
fi


#### calculate resolution based average bias score
if [ ! -e $rebed.$tag.fl ]
then
	## make window
	if [ ! -e $chrom.$tag.sortBed ]
	then
		bedtools makewindows -g $chrom -w $res |awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"\t"$1"_"$2}' > $chrom.$tag
    bedtools sort -i $chrom.$tag > $chrom.$tag.sortBed
	else
		echo "#*_*# <$chrom.$tag> already exists...Skip~~~"
	fi
	## average by window (NA: the bin doesn't include any fragment)
	# fragment length
	bedtools map -a $chrom.$tag.sortBed -b $rebed.frag.fl.sortBed -c 5 -o mean -null NA > $rebed.$tag.fl
	# make the nature chrom order
	awk 'NR==FNR {a[$1"_"$2"_"$3]=$0} NR>FNR {print a[$1"_"$2"_"$3]}' $rebed.$tag.fl $chrom.$tag > $$; mv $$ $rebed.$tag.fl
	# format
	cut -f5- $rebed.$tag.fl > $$; mv $$ $rebed.$tag.fl

else
	echo "#*_*# <$rebed.$tag.fl> already exists...Skip~~~"
fi









