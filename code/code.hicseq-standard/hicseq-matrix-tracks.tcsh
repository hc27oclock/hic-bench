#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-tracks.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)          # TODO: allow multiple samples, so that this operation can be run by-group

# test number of input objects
set object = $objects[1]
if ($#objects != 1) then
  scripts-send2err "Error: this operation allows only one input object!"
  exit 1
endif

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir enzyme bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# distance-normalize matrix for each chromosome
set inpdir = $branch/$object
echo -n "" >! $outdir/track.washu.tsv
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv`)
  scripts-send2err "Processing $mat..."

  set f = $inpdir/$mat
  
  # inverse-rotate matrix if necessary
  if (`cat $f | head -1 | grep 'chr.*chr' | wc -l` == 1) then
    set fnew = $f																											# "full" matrix
  else
    set fnew = $outdir/tmp-matrix.tsv
    Rscript ./code/hic-matrix.r rotate -v --row-labels -o $fnew $f		# inverse-rotate matrix
  endif

  # convert  
  cat $fnew | gtools-hic convert -v --col-labels -t '	' >> $outdir/track.washu.tsv
  
  # cleanup
  rm -f $outdir/tmp-matrix.tsv
end

# Compress and index
bgzip $outdir/track.washu.tsv
tabix -p bed $outdir/track.washu.tsv.gz

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


