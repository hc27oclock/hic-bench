#!/bin/tcsh

source ./inputs/params/params.tcsh

# unload potentially conflicting modules
module unload gcc
module unload python

# macs is part of python/2.7.3 module (used to be macs/2.0.10.20131216)
module load python/2.7.3

set extsize = `./code/read-sample-sheet.tcsh $sheet "$objects" fragmentation-size`
set extsize = `echo $extsize | tools-vectors m -n 0`
set caller_params = "--nomodel --extsize=$extsize --qvalue 0.05"
set use_input = 'true'
set annot_params = "annotate2 -i --upstream-dist 100000 --downstream-dist 100000 --proximal-dist 1000 $genome_dir/gene-name.bed"
