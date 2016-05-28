#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'       # excluded chromosomes

set diff_params = "--min-count=50 --min-zscore=2.0 --min-distance=40000"


