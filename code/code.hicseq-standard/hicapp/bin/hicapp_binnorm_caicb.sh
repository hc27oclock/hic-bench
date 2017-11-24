#!/bin/bash
####
# caICB correction on Hi-C data
####

## set path of the pipeline
bindir=`dirname $0`

## set par
pair=$1
icb=$2
K=$3 #parameter K in caICB method
res=$4 #resolution
outiqr=10 #outiqr
ld=$(expr $res + 1)
hd=$(expr $K \* $res)


if [ ! -e $pair.selected.allpairs.$K ]
then
	## pre-select contact pairs to reduce memory useage for R
	sed 's/_/\t/g' $pair |awk -v ld=$ld -v hd=$hd 'function abs(v) {return (v<0)?-v:v} $1==$3 && abs($4-$2)>=ld && abs($4-$2)<=hd {print $1"_"$2"\t"$3"_"$4"\t"$5"\t"abs($4-$2)}' > $pair.selected.$K

	## get all possible pairs, and add stats into the pair file
	$bindir/hicapp_binnorm_getallpairs.pl $pair.selected.$K $icb $ld $hd > $pair.selected.allpairs.$K

	## clean tmp
	rm $pair.selected.$K
else
	echo "#*_*# <$pair.selected.allpairs.$K> already exists...Skip~~~"
fi


## caICB correction
if [ ! -e $pair.selected.allpairs.$K.adj ]
then
	Rscript $bindir/hicapp_binnorm_caICB.R $pair.selected.allpairs.$K $icb $outiqr
else
	echo "#*_*# <$pair.selected.allpairs.$K.adj> already exists...Skip~~~"
fi



