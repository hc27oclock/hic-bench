#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'       # excluded chromosomes

set loop_params = "--bin-size=$bin_size --lambda-id=1 --min-count=20 --min-zscore=2.0 --min-distance=40000 --max-distance=10000000"        # parameters for identifying significant interactions


