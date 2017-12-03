#!/bin/tcsh

source ./inputs/params/params.tcsh

set aligner = gtools
set genome = `./code/read-sample-sheet.tcsh $sheet $object genome`
set genome_index = inputs/genomes/$genome/bowtie2.index/genome
set align_params = "--min-len 30 --len-diff 5"

