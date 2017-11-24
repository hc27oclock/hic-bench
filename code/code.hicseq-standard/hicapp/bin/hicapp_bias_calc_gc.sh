####
# calculate GC bias for Hi-C data
####

#### set path of the pipeline
bindir=`dirname $0`


#### get pars
rebed=$1 #rebed=/michorlab/hjwu/code/hic_pipeline/ref/hg19/Digest_hg19_MboI_None.clear 
chrom=$2 #chrom=/michorlab/hjwu/code/hic_pipeline/ref/hg19/hg19.genome
fasta=$3 #fasta=/michorlab/hjwu/code/hic_pipeline/ref/hg19/hg19.fa
res=$4 #res=1000000
tag=$5 #tag=1m
gcflank=$6 #gcflank=200


#### calculate fragment based bias score
if [ ! -e $rebed.frag.gc.sortBed ]
then
	#### GC content
  # calculate fragment bed file
  awk '{print $1"\t"$2"\t"$3"\t"$1"_"int(($2+$3)/2)}' $rebed |sortBed > $rebed.frag
	# get the enzyme cutsites
	awk '{print $1"\t"$2"\t"$2+1"\t"$1"_"$2"\n"$1"\t"$3"\t"$3+1"\t"$1"_"$3}' $rebed |uniq > $rebed.cutsite
	# get the N bp flanking region of cutsites
	bedtools flank -i $rebed.cutsite -g $chrom -b $gcflank |awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2}' |awk '$3>$2' |sortBed |uniq > $rebed.cutsite.flank${gcflank}bp
	# calculate GC content for each cutsite
	bedtools nuc -fi $fasta -bed $rebed.cutsite.flank${gcflank}bp > $rebed.cutsite.flank${gcflank}bp.nuc
	# calculate GC content for each fragment
	bedtools map -a $rebed.frag -b $rebed.cutsite.flank${gcflank}bp.nuc -c 6 -o mean -null 0 > $rebed.frag.gc

	#### sort files
	bedtools sort -i $rebed.frag.gc > $rebed.frag.gc.sortBed

	#### clear intermediate files
	rm $rebed.frag $rebed.cutsite $rebed.cutsite.flank${gcflank}bp $rebed.cutsite.flank${gcflank}bp.nuc $rebed.frag.gc

else
	echo "#*_*# <$rebed.frag.gc.sortBed> already exists...Skip~~~"
fi


#### calculate resolution based average bias score
if [ ! -e $rebed.$tag.gc ]
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
	# effective length
	bedtools map -a $chrom.$tag.sortBed -b $rebed.frag.gc.sortBed -c 5 -o mean -null NA > $rebed.$tag.gc
	# make the nature chrom order
	awk 'NR==FNR {a[$1"_"$2"_"$3]=$0} NR>FNR {print a[$1"_"$2"_"$3]}' $rebed.$tag.gc $chrom.$tag > $$; mv $$ $rebed.$tag.gc
	# format
	cut -f5- $rebed.$tag.gc > $$; mv $$ $rebed.$tag.gc
 
else
	echo "#*_*# <$rebed.$tag.gc> already exists...Skip~~~"
fi









