#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-boundary-diff.tcsh OUTDIR PARAMS BRANCH OBJECT1 OBJECT2
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

foreach f (`cd $branch/$object1; ls -1 matrix.*.tsv matrix.*.RData | grep -vwE "$chrom_excluded"`)
  set chr = `echo $f | cut -d'.' -f2`
  scripts-send2err "Processing matrix $f..."
  if ((-e $branch/$object1/$f) && (-e $branch/$object2/$f)) then
    Rscript ./code/hic-matrix.r bdiff -v -o $outdir/$chr --row-labels $diff_domains_params $branch/$object1/$f $branch/$object2/$f             # TODO: use scripts-qsub-run to assign memory
  endif
end

# Collect diff-domains from all chromosomes
set K = `ls -1 $outdir/*/bdiff.k=*.tsv | sed 's/.*\/bdiff.k=//' | cut -d'.' -f1 | sort -u`
foreach k ($K)
  foreach t (bdiff.k=$k.tsv all_data.k=$k.tsv)
    cat $outdir/*/$t | head -1 >! $outdir/$t
    foreach tt ($outdir/*/$t)
      cat $tt | sed '1d' >> $outdir/$t
    end
  end
  set score_col = `head -1 $outdir/bdiff.k=$k.tsv | tr '\t' '\n' | grep -n '^ratio-zdiff$' | cut -d':' -f1`
  cat $outdir/bdiff.k=$k.tsv | sed '1d' | cut -f1,$score_col | awk '$2>0' | sed 's/:/\t/' | sed 's/-/\t/' | gtools-regions shiftp -5p -1 -3p 0 >! $outdir/boundary_gain.k=$k.bed
  cat $outdir/bdiff.k=$k.tsv | sed '1d' | cut -f1,$score_col | awk '$2<0' | sed 's/:/\t/' | sed 's/-/\t/' | gtools-regions shiftp -5p -1 -3p 0 >! $outdir/boundary_loss.k=$k.bed
end

# cleanup
rm -rf $outdir/*/*.RData $outdir/*/*.tsv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




