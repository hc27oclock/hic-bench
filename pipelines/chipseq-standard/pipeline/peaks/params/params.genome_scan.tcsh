#!/bin/tcsh

# load the custom tcsh environment
source ./inputs/params/params.tcsh

# HPC modules to load; this adds/removes these programs in the PATH
module unload gcc
module unload samtools
module load samtools/1.3
module load homer/v4.6

set use_input = 'true'

set ref = "$genome_dir/genome.bed"

set caller_params = "-s 5000 -d 1000"

set min_fc="1"

set annot_params = "annotate2 -i --upstream-dist 100000 --downstream-dist 100000 --proximal-dist 1000 $genome_dir/gene-name.bed"
