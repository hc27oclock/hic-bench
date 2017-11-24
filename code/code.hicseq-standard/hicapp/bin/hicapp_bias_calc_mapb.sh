####
# calculate mappability bias for Hi-C data
####

#### set path of the pipeline
bindir=`dirname $0`


#### get pars
rebed=$1 #rebed=/michorlab/hjwu/code/hic_pipeline/ref/hg19/Digest_hg19_MboI_None.clear 
chrom=$2 #chrom=/michorlab/hjwu/code/hic_pipeline/ref/hg19/hg19.genome
res=$3 #res=1000000
tag=$4 #tag=1m
mapb=$5 #mapb=$bdir/ref/hg19/mappability100mer.bw
mapbflank=$6 #mapbflank=500


#### calculate fragment based bias score
if [ ! -e $rebed.frag.mapb.sortBed ]
then
	#### get the mappability score
  # calculate fragment bed file
  awk '{print $1"\t"$2"\t"$3"\t"$1"_"int(($2+$3)/2)}' $rebed |sortBed > $rebed.frag
  # get the enzyme cutsites
	awk '{print $1"\t"$2"\t"$2+1"\t"$1"_"$2"\n"$1"\t"$3"\t"$3+1"\t"$1"_"$3}' $rebed |uniq > $rebed.cutsite
	# get the N bp flanking region of cutsites
	bedtools flank -i $rebed.cutsite -g $chrom -b $mapbflank |awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2}' |awk '$3>$2' |sortBed |uniq > $rebed.cutsite.flank${mapbflank}bp
	# split large bed file into pieces to fit bigWigAverageOverBed
	mkdir $rebed.cutsite.flank${mapbflank}bp.cut
	split -a 5 -d -l 2900 $rebed.cutsite.flank${mapbflank}bp $rebed.cutsite.flank${mapbflank}bp.cut/cut
	# loop files to calculate mappability for each cutsite
	for files in $rebed.cutsite.flank${mapbflank}bp.cut/cut*; do
		$bindir/bigWigAverageOverBed -bedOut=$files.map100mer.bed $mapb $files $files.map100mer
	done
	# merge all pieces into one
	cat $rebed.cutsite.flank${mapbflank}bp.cut/cut*.map100mer.bed > $rebed.cutsite.flank${mapbflank}bp.map
	# calculate mappability for each fragment
	bedtools map -a $rebed.frag -b $rebed.cutsite.flank${mapbflank}bp.map -c 5 -o mean -null 0 > $rebed.frag.mapb
	# clean tmp files
	rm -rf $rebed.cutsite.flank${mapbflank}bp.cut

	#### sort files
	bedtools sort -i $rebed.frag.mapb > $rebed.frag.mapb.sortBed

	#### clear intermediate files
	rm $rebed.frag $rebed.cutsite $rebed.cutsite.flank${mapbflank}bp $rebed.cutsite.flank${mapbflank}bp.map $rebed.frag.mapb

else
	echo "#*_*# <$rebed.frag.mapb.sortBed> already exists...Skip~~~"
fi


#### calculate resolution based average bias score
if [ ! -e $rebed.$tag.mapb ]
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
	bedtools map -a $chrom.$tag.sortBed -b $rebed.frag.mapb.sortBed -c 5 -o mean -null NA > $rebed.$tag.mapb
	# make the nature chrom order
	awk 'NR==FNR {a[$1"_"$2"_"$3]=$0} NR>FNR {print a[$1"_"$2"_"$3]}' $rebed.$tag.mapb $chrom.$tag > $$; mv $$ $rebed.$tag.mapb
	# format
	cut -f5- $rebed.$tag.mapb > $$; mv $$ $rebed.$tag.mapb

else
	echo "#*_*# <$rebed.$tag.mapb> already exists...Skip~~~"
fi









