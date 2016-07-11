#!/bin/tcsh

# unload tools that may cause conflicts with other tools
module unload gcc
module unload samtools
module unload java
module unload r
module unload python

# load all tools
module load kentutils/329
module load samtools/1.2.1
module load bedtools/2.22.0
module load java/1.7
picard-tools/1.88
module load bowtie2
# module load macs/2.0.10.20131216
module load r/3.2.0

# set sample sheet
set sheet = inputs/sample-sheet.tsv

