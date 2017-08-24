#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-optimal-lambda.tcsh OUTDIR PARAM-SCRIPT COMPARE-MATRICES-BRANCH [OBJECTS]
##

if (($#argv < 3) || ($#argv > 4)) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# check
if ($#objects > 1) then
  scripts-send2err "Error: only 1 input object allowed for this operation."
  exit
endif

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

scripts-send2err "Calculating optimal lambda..."
set inpfile = $branch/$objects[1]/compare.tsv
Rscript code/run-optimal-lambda.r $outdir $inpfile $min_improvement >! $outdir/results.tsv
cat $outdir/results.tsv | awk -v p=$max_pvalue '$4<=p' | tail -1 >! $outdir/optimal.tsv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




