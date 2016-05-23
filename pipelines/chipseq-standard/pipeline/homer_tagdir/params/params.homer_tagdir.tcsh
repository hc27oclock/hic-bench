#!/bin/tcsh

# load the custom tcsh environment
source ./inputs/params/params.tcsh

# HPC modules to load; this adds/removes these programs in the PATH
module unload gcc
module load samtools/1.3
module load homer/v4.6

# some old exaple params from peaks step..
# set extsize = `./code/read-sample-sheet.tcsh $sheet "$objects" fragmentation-size`
# set extsize = `echo $extsize | tools-vectors m -n 0`

# what does this do??
set use_input = 'true'

