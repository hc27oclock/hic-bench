#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-estimated.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT(S)
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

# setup
./code/create-ignored-loci.tcsh $genome_dir $bin_size >! $outdir/ignored_loci.txt

module unload r
module load r/3.3.0
# run estimation
set inpdir = $branch/$object
#set matrix_sum = `cat $inpdir/matrix.chr*.tsv | cut -f2- | grep -v ^chr | tr '\t' ' ' | tools-vectors sum -v | tools-matrix csum`
#set scale = `echo 1000000000/$matrix_sum | bc -l`
set scale = 1.0
set jid = ()
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv | grep -vwE "$chrom_excluded"`)
  scripts-send2err "Processing input matrix $mat..."
  set mem1 = `./code/calc-matrix-memory.tcsh $inpdir/$mat 2.0 40 | sed 's/G$//'` 
  set mem2 = `./code/calc-matrix-memory.tcsh $inpdir/$mat 1.0 0 | sed 's/G$//'`
  set mem = `echo "$mem1+$mem2" | bc`G 
  scripts-send2err "requested memory = $mem"
  set outmat = `echo $mat | sed 's/.tsv$/.RData/'`
  set jpref = $outdir/__jdata/job.$outmat
  scripts-create-path $jpref
  set jid = ($jid `scripts-qsub-run $jpref 1 $mem ./code/run-fastMeanFilter.tcsh $prep $smooth_params $inpdir/$mat $outdir/$outmat`)
end
scripts-qsub-wait "$jid"

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


