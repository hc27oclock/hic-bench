#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MYX]'       # excluded chromosomes

set max_dist = `echo 2500000/$bin_size | bc`      # number of bins (max distance = 2.5Mb)

set compare_params = "--max-dist=$max_dist --n-dist=5 --min-lambda=0.0 --max-lambda=1.0 --n-lambda=6 --gamma=0"           # lambda/gamma are only used if estimation was done with max-lambda=Inf

set prep = none

