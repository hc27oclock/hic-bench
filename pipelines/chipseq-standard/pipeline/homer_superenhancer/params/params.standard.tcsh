#!/bin/tcsh

# load the custom tcsh environment
source ./inputs/params/params.tcsh

# HPC modules to load; this adds/removes these programs in the PATH
module unload gcc
module load samtools/1.3
module load homer/v4.6


# whether to call peaks against ChIP-Seq INPUT sample files specified in the sample sheet (aka the control samples)
set use_input = 'true'

# http://homer.salk.edu/homer/ngs/peaks.html


# other params to use go here
set extra_params = ""
# e.g.
# set extra_params = "-size <#> -F <#>"
