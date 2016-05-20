#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-annotations-stats.tcsh OUTPUT-DIR PARAM-SCRIPT ANNOTATIONS-BRANCH OBJECTS
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

if ($#objects != 1) then
  scripts-send2err "Error: hicseq-annotations-stats.tcsh requires exactly one input object."
  exit 1
endif

# compute loci counts and frequencies
set n_bins = `cat $genome_dir/genome.bed | gtools-regions win -s $bin_size -d $bin_size | wc -l`
set loci_reg = $branch/$objects[1]/loci.reg
cat $loci_reg | cut -d' ' -f3 | sed 's/$/\/'$bin_size'/' | bc | paste - $loci_reg | cut -d' ' -f1 | sort -u | cut -f2 | sort | uniq -c | sed "s/^ */$n_bins /" | tools-cols 2 0 1 | sed 's/ /\t/' | tools-vectors div -n 6 >! $outdir/freq.tsv

# compute stats
set score_column = 3    # raw unprocessed score (directly from input matrix; no scaling, no dist-norm)
Rscript ./code/hicseq-annotations-enrichments.r $outdir $branch/$objects[1]/table.annotated.tsv $outdir/freq.tsv $nbest $score_column

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


