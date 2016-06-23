#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                                                      # excluded chromosomes
set diff_domains_params = ( \
--min-lambda=0.0 --max-lambda=1.0 --n-lambda=6 --gamma=0 \
--preprocess=none \
--method=ratio \
--distance=`echo 1000000/$bin_size | bc` \
--distance2=`echo 1000000/$bin_size | bc` \
--skip-distance=2 \
--flank-dist=`echo 500000/$bin_size | bc` \
--slope=1.1 \
--fdr=0.01 \
--z1=1.0 \
--z2=0.5 \
--track-dist=`echo 2000000/$bin_size | bc` \
)


