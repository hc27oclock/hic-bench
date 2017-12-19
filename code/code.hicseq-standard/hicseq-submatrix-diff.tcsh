#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-submatrix-diff.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT1 OBJECT2
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set object1 = $4
set object2 = $5

set objects = ($object1 $object2)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# code dir
set code_dir = code/hicseq-submatrix-diff-scripts

# perform Hi-C submatrix differential analysis
Rscript $code_dir/find-hic-squares.r -v -o $outdir --bin-size $bin_size $submatrix_diff_params $branch/$object1 $branch/$object2
cat $outdir/out.tsv | scripts-skipn 1 | tr '\t' ' ' | sed 's/ /|/' | sed 's/ /|/' | sed 's/ /|/' | sed 's/ /|/' | sed 's/ /\t/' | gtools-overlaps bin -v -i --print-labels --print-regions coding-tss-bins.bed | tr '|' '\t' >! $outdir/out1.tsv


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


