#!/bin/tcsh

source ./inputs/params/params.tcsh

set tool = hicmatrix

set chrom_excluded = 'chr[MY]'

set d = 500000

set hicmatrix_params = ( \
--min-lambda=0.0 --max-lambda=1.0 --n-lambda=6 --gamma=0 \
--preprocess=none \
--method=ratio \
--distance=`echo $d/$bin_size | bc` \
--distance2=`echo $d/$bin_size | bc` \
--skip-distance=0 \
--flank-dist=`echo $d/$bin_size | bc` \
--slope=1.2 \
--fdr=0.10 \
--track-dist=`echo 2000000/$bin_size | bc` \
--presentation=none \
)


