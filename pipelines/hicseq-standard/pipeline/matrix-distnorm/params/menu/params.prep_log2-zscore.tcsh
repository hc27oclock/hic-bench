#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'       # excluded chromosomes

set preprocess = log2-zscore

set max_dist = `echo 10000000/$bin_size | bc`          # number of bins (max distance = 10Mb)
set distnorm_params = "--max-dist=$max_dist"


