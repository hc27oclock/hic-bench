#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                                                      # excluded chromosomes

set boundary_scores_params = ( \
--min-lambda=0.0 --max-lambda=1.0 --n-lambda=6 --gamma=0 \
--preprocess=none \
--distance=`echo 2000000/$bin_size | bc` \
--distance2=1 \
--skip-distance=2 \
--flank-dist=`echo 500000/$bin_size | bc` \
)


