#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-compare-matrices.tcsh OUTDIR PARAMS BRANCH OBJECT1 OBJECT2
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

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$object1 $object2" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# compare matrices (one chromosome at a time)
set jid = ()
foreach f (`cd $branch/$object1; ls -1 matrix.*.tsv matrix.*.RData | grep -vwE "$chrom_excluded"`)
  set chr = `echo $f | cut -d'.' -f2`
  scripts-send2err "Processing matrix $f..."
  if ((-e $branch/$object1/$f) && (-e $branch/$object2/$f)) then
    set jpref = $outdir/__jdata/job.$chr
    scripts-create-path $jpref
    set chr_size = `cat $genome_dir/genome.bed | grep "^$chr	" | gtools-regions n | cut -f2`
    set n_bins = `echo $chr_size/$bin_size+1 | bc`
    set mem = `echo "100*2.0*$n_bins*$n_bins/1000000000+5" | bc`G                 # TODO: for estimated RData matrices, memory should take into account the number of lambdas...
    scripts-send2err "-- requested memory = $mem"
    set jid = ($jid `scripts-qsub-run $jpref 1 $mem Rscript ./code/hic-matrix.r compare -v -o $outdir/$chr $compare_params $branch/$object1/$f $branch/$object2/$f`)
  endif
end
scripts-send2err "Waiting for all jobs [$jid] to complete..."
scripts-qsub-wait "$jid"

# Collect all correlation coefficients along with sample and lambda info
set comp = `basename $outdir`
set methods = `cd $outdir; ls -1 *.cor.*.tsv | sed 's/[^.]\+\.cor\.//' | sed 's/\.tsv$//' | sort -u`
set header = `cd $outdir; cat *.cor.*.tsv | head -1 | cut -f2-`
foreach method ($methods)
  echo "SAMPLE-1 SAMPLE-2 COMPARISON METHOD CHROMOSOME LAMBDA $header" | tr ' ' '\t' >! $outdir/cor.$method.tsv
  foreach f ($outdir/*.cor.$method.tsv)
    set chr = `basename $f | cut -d'.' -f1`
    cat $f | scripts-skipn 1 | sed "s/^/$object1\t$object2\t$comp\t$method\t$chr\t/" >> $outdir/cor.$method.tsv
  end
  rm -f $outdir/*.cor.$method.tsv
end 

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




