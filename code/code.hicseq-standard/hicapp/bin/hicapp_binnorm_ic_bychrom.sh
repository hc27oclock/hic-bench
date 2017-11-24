#!/bin/bash
####
# run IC normalization by chromosomes
####


#### set path of the pipeline
bindir=`dirname $0`


#### set pars
matdir=$1 #bychrom matrix dir
cmd=$2 #ic or ic_mes or ic_mep
ncpu=$3 #number cpu
iter=$4 #number of iteration
out=$5 #output file name

#### loop chrom to run IC normalization
for files in $matdir/*.mat; do
	## IC normalization
	$bindir/hicapp_binnorm_ic.sh $files $cmd $ncpu $iter
done


#### merge IC bias vector for all chromosome
cat $matdir/*.icb > $out

