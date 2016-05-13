#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                                                      # excluded chromosomes

set boundary_scores_params = ( \
--min-lambda=0.0 --max-lambda=1.0 --n-lambda=6 --gamma=0 \
--preprocess=log2 \
--distance=`echo 500000/$bin_size | bc` \
--distance2=`echo 500000/$bin_size | bc` \
--skip-distance=0 \
--flank-dist=`echo 500000/$bin_size | bc` \
--tolerance=0.01 \
)


