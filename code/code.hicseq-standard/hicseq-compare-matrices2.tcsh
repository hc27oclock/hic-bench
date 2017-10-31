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
    # setup job
    set jpref = $outdir/__jdata/job.$chr
    scripts-create-path $jpref
    set chr_size = `cat $genome_dir/genome.bed | grep "^$chr	" | gtools-regions n | cut -f2`
    set n_bins = `echo $chr_size/$bin_size+1 | bc`
    set mem = `echo "100*2.0*$n_bins*$n_bins/1000000000+5" | bc`G
    scripts-send2err "-- requested memory = $mem"
    
    # compare matrices
    set jid = ($jid `scripts-qsub-run $jpref 1 $mem ./code/run-hicrep.tcsh $outdir/out.$chr $branch/$object1/$f $branch/$object2/$f $bin_size $max_dist $prep`)
  endif
end
scripts-send2err "Waiting for all jobs [$jid] to complete..."
scripts-qsub-wait "$jid"

# Collect comparisons of all chromosomes into one single file
echo "SAMPLE-1 SAMPLE-2 CHROMOSOME LAMBDA VALUE" | tr ' ' '\t' >! $outdir/compare.tsv
cat $outdir/out.*/compare.tsv >> $outdir/compare.tsv
# Clean-up
rm -rf $outdir/out.*

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




