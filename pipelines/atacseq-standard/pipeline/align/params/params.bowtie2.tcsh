#!/bin/tcsh

source ./inputs/params/params.tcsh

set aligner = bowtie2
set genome = `./code/read-sample-sheet.tcsh $sheet $object genome`
set genome_index = inputs/genomes/$genome/bowtie2.index/genome
set align_params = "--local -x $genome_index --minins 50 --maxins 2000 --no-mixed --no-discordant "
set mapq = 30
set release = inputs/release
set excluding_regions = "~gongy05/references/$genome/blacklist.bed"
set excluding_chrom = "chrM"
