#!/bin/tcsh

# load the custom tcsh environment
source ./inputs/params/params.tcsh

# HPC modules to load; this adds/removes these programs in the PATH
module unload gcc
module unload samtools
module load samtools/1.3
module load homer/v4.6

# some old exaple params from peaks step..
# set extsize = `./code/read-sample-sheet.tcsh $sheet "$objects" fragmentation-size`
# set extsize = `echo $extsize | tools-vectors m -n 0`

# what does this do??
set use_input = 'true'



# -style <option> (Specialized options for specific analysis strategies)
#                        factor (transcription factor ChIP-Seq, uses -center, output: peaks.txt,  default)
#                        histone (histone modification ChIP-Seq, region based, uses -region -size 500 -L 0, regions.txt)
#                        groseq (de novo transcript identification from GroSeq data, transcripts.txt)
#                        tss (TSS identification from 5' RNA sequencing, tss.txt)
#                        dnase (Hypersensitivity [crawford style (nicking)], peaks.txt)
#                        super (Super Enhancers, superEnhancers.txt)
#                        mC (Cytosine methylation (BS-seq/methylC-seq), regions.txt)

if (`echo $objects | tr ' ' '\n' | tr '-' '\n' | grep -iEc '^H2AZ|^H[234]K[0-9]' | sort | uniq` == 1) then
  set peak_style = 'histone'   # histone
else
  set peak_style = 'factor'
endif
set caller_params = $peak_style

# other params to use go here
set extra_params = ""



set annot_params = "annotate2 -i --upstream-dist 100000 --downstream-dist 100000 --proximal-dist 1000 $genome_dir/gene-name.bed"
