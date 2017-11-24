#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-diff.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT1 OBJECT2
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
set codedir = code/hicseq-domains-diff-scripts

# first, determine the TADs to be used
set domains1 = $domains_branch/$object1/domains.k=001.bed
set domains2 = $domains_branch/$object2/domains.k=001.bed

# find consistent TADs between samples
perl $codedir/find-consistent-domains.pl $domains1 $domains2 $max_boundary_dist $bin_size $outdir/domains1.tsv $outdir/domains2.tsv

# perform Hi-C fold-change analysis
$codedir/run_comparison.sh $codedir/differential_tad_activity.r $branch/$object1 $branch/$object2 $outdir/domains1.tsv $outdir/domains2.tsv $min_tad_size $max_tad_size $is_normalize $centrotelo_file $bin_size $outdir

# create figures
R --no-save $outdir/final_results.tsv $gene_name FALSE FALSE \
        $object1 $object2 $bin_size $min_tad_size $outdir/final_results < $codedir/differential_tad_activity_expression.r



# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


