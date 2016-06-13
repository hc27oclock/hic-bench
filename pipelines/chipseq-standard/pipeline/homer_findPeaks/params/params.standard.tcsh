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

# the style of peak calling to use
# set peak_style = 'factor' 
set peak_style = 'histone' 

# -style <option> (Specialized options for specific analysis strategies)
#                        factor (transcription factor ChIP-Seq, uses -center, output: peaks.txt,  default)
#                        histone (histone modification ChIP-Seq, region based, uses -region -size 500 -L 0, regions.txt)
#                        groseq (de novo transcript identification from GroSeq data, transcripts.txt)
#                        tss (TSS identification from 5' RNA sequencing, tss.txt)
#                        dnase (Hypersensitivity [crawford style (nicking)], peaks.txt)
#                        super (Super Enhancers, superEnhancers.txt)
#                        mC (Cytosine methylation (BS-seq/methylC-seq), regions.txt)

# other params to use go here
set extra_params = ""
# e.g.
# set extra_params = "-size <#> -F <#>"