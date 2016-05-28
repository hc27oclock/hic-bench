#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-interactions.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT(S)
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
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir enzyme bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# create table of interactions for each chromosome
set inpdirs = `echo $objects | tr ' ' '\n' | awk -v v=$branch '{print v "/" $1}'`
set objects_list = `echo $objects | tr ' ' ','`
set mat_list = `cd $inpdirs[1]; ls -1 matrix.*.tsv matrix.*.RData | grep -vwE "$chrom_excluded"`
set n_reads =
foreach inpdir ($inpdirs)
  if (-e $inpdir/stats.tsv) then
    set n_reads = ($n_reads `cat $inpdir/stats.tsv | grep '^ds-accepted-intra	' | cut -f2`)
  else
    set n_reads = ($n_reads 1)
  endif
end
set n_reads = `echo $n_reads | tr ' ' ','`
scripts-send2err "number of reads for scaling = $n_reads"
set jid = 
foreach mat ($mat_list)
  scripts-send2err "Processing input matrix $mat..."
  set mat_files = `echo $inpdirs | tr ' ' '\n' | awk -v v=$mat '{print $1 "/" v}'`
  set outmatdir = `echo $mat | sed 's/\.tsv$//' | sed 's/\.RData$//'`
  set a = `echo "3*$#inpdirs" | bc`
  set mem = `./code/calc-matrix-memory.tcsh $inpdirs[1]/$mat $a 5`
  scripts-send2err "requested memory = $mem"
  set jdata = $outdir/__jdata
  scripts-create-path $jdata/
  set jpref = $jdata/job.$outmatdir
  set jid = ($jid `scripts-qsub-run $jpref 1 $mem Rscript ./code/hic-matrix.r loops -v -o $outdir/$outmatdir -L $objects_list --n-reads=$n_reads $loop_params $mat_files`)
end
scripts-qsub-wait "$jid"

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


