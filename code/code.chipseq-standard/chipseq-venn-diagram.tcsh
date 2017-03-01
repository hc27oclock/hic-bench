#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: chipseq-venn-diagram.tcsh OUTPUT-DIR PARAMETER-SCRIPT PEAKS-BRANCH [OBJECTS]
##

# process command-line inputs
if (($#argv < 3) || ($#argv > 4)) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

# inputs
set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# if samples is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
scripts-send2err "Setting parameters..."
source $params
scripts-send2err "-- Parameters: "


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------


# extract group level peak calling result
set group_branch = `echo $branch | sed -E 's/peaks.by_sample/peaks.by_group/'`
set group_name = `basename $outdir`

# create merged peaks reference
set peaks = `echo $objects | tr ' ' '\n' | awk -v b=$branch '{print b"/"$1"/peaks.bed"}'`

# plot annotations

scripts-send2err "-- CMD: "

scripts-send2err "Rscript --vanilla code/chipseq-venn-diagram.r -o $outdir $peaks $group_branch/$group_name/peaks.bed"
#Rscript --vanilla code/chipseq-venn-diagram.r -o $outdir $peaks $group_branch/$group_name/peaks.bed

foreach file ($bed_files)
  set filename = `basename $file:r`.pdf
  scripts-send2err "Rscript --vanilla code/chipseq-venn-diagram.r -o $outdir -f $filename $peaks $group_branch/$group_name/peaks.bed $file"
  Rscript --vanilla code/chipseq-venn-diagram.r -o $outdir -f $filename $peaks $group_branch/$group_name/peaks.bed $file
end
 
# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


