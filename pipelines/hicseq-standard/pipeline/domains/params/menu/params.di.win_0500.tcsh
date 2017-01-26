#!/bin/tcsh

#source ./inputs/params/params.tcsh

# Load required modules
module unload r
module unload java
module load java/1.6
module unload gcc
module load gcc/4.4
#module unload matlab
module load matlab/R2013a

module list

set tool = di
set domaincallpath = "/ifs/home/cl3011/ROTATION_3/Resources/Software/domaincall_software"
set chrom_excluded = 'chr[MY]'
set window_size = 500000 
set min = 2
set prob = 0.99
set faipath = "/ifs/home/cl3011/hicbench-new-analysis/inputs/genomes/hg19/bowtie2.index/genome.fa.fai"
set bin_size = 40000
