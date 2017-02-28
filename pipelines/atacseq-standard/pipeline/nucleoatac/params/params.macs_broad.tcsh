#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload gcc
module unload python
# macs now part of python module
# module load macs/2.0.10.20131216
module load python/2.7.3

set nucleoatac = ""
set nucleoatac_path = '/ifs/home/gongy05/utilities/NucleoATAC/bin/nucleoatac'
set shift_dist = "100"
set genome_fa = $genome_dir/bowtie2.index/genome.fa
