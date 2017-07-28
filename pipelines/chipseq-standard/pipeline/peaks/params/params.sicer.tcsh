#!/bin/tcsh

source ./inputs/params/params.tcsh

module load python/2.7.11-gcc-4.7.0

set extsize = `./code/read-sample-sheet.tcsh $sheet "$objects" fragmentation-size`
set extsize = `echo $extsize | tools-vectors m -n 0`
if (`echo $objects | tr ' ' '\n' | tr '-' '\n' | grep -iEc '^H2AZ|^H[234]K[0-9]' | sort | uniq` == 1) then
  set win_size = '200'	# histone
else
  set win_size = '100'	# TF
endif
set caller_params = "1 $win_size $extsize 0.74 `echo 3\*$win_size | bc` 0.01" 
# caller_params:
# [genome] [redundancy threshold] [windwo size (bp)] [fragment size (bp)] [effective genome fraction] [gap size (bp)] [FDR]
# 	- genome: mm8, mm9, rn4, hg18, hg19, sacCer1, dm2, dm3, prmbe, tair8. If the desired reference genome is not in the list, user can easily add customized reference genome information in the file GenomeData.py under / lib. For example, if a new reference genome “NewGenome” contains two chromosomes “chr1” and “chrX,” with length 100 and 200, respectively, user will need to:
# 		a. Add this entry to dictionary “species_chroms”: “NewGenome”: NewGenome_chroms
# 		b. Add this entry to dictionary “species_chrom_lengths”: “NewGenome”: NewGenome_chrom_length
# 		c. Add this list to GenomeData.py: NewGenome_chroms=[“chr1,” “chrX”]
# 		d. Add this dictionary to GenomeData.py: NewGenome_chrom_lengths={“chr1”: 100, “chrX”: 200} 
# 	- redundancy threshold: number of redundant reads kept for analysis. In this example it is 1. Redundant reads refer to reads with exactly the same genomic location and orientation. For typical ChIP-Seq datasets, this is likely due to PCR amplification artifact. To remove this potential bias, we generally recommend removing the redundancy and retaining only 1 read for each set of redundant reads.
# 	- windwo size (bp): resolution of SICER algorithm, the width (in bps) of window in comparing ChIP with control library. A window too narrow will exaggerate local fluctuation in each window, while a window too large will cause over-smoothing of data and lose resolution. For histone modifications, one can use 200 bp (a number approximately the length of a single nucleosome and a linker), for transcription factors, one can use 50-100 bps.
# 	- fragment size (bp): the average size (in bps) of ChIP fragment, used to assign a ChIP read to the center of the DNA fragment. Typical sonication outputs ChIP fragment of 150–300 bps long. Recommended FRAGMENT_SIZE=150, means the shift is 75.
# 	- effective genome fraction: Effective Genome as fraction of the genome size. It depends on read length, 0.8, which is recommended for single-end ChIP-Seq with 50 bp read length. It can be found or computed from Uniqueome.
# 	- gap size (bp): the gap size (in bps) allowed in SICER filtering. needs to be multiples (multiply by 3 work better) of window size. 600 (200*3) as recommended for extended histone modification marks like H3K27me3. For more careful consideration, users can plot the aggregate score of all significant islands (output: *.scoreisland) as a function of g. If the aggregate score reaches a maximum inside of the range of g explored, the gap size corresponding to the highest aggregate score would be a good choice.
# 	- FDR: cut off
# 		a. for without control, this is E-value (use 100 if you are trying to use 0.01 in FDR) (poorly-documented)
set use_input = 'true'
set annot_params = "annotate2 -i --upstream-dist 100000 --downstream-dist 100000 --proximal-dist 1000 $genome_dir/gene-name.bed"

