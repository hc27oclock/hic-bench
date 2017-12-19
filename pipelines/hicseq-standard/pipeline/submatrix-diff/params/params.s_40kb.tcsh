#!/bin/tcsh

source ./inputs/params/params.tcsh

set submatrix_diff_params = (                  \
--pseudo-count=10                              \
--square-size=`echo 400000/$bin_size | bc`     \
--max-dist=`echo 4000000/$bin_size | bc`       \
--logfc-cutoff=0.2                             \
--qval-cutoff=1e-4                             \
)



