#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-ic.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT(S)
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
if (! $?NSLOTS) then
  set threads = 16
else
  set threads = $NSLOTS
endif
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# run iteractive correction
set inpdir = $branch/$object
set jid = 
set matrices = `cd $inpdir; ls -1 matrix.*.tsv | grep -vwE "$chrom_excluded"`
foreach mat ($matrices)
  scripts-send2err "Processing $mat..."
  
  # estimate required memory
  set mem = `./code/calc-matrix-memory.tcsh $inpdir/$mat 1.5 5`
  scripts-send2err "requested memory = $mem"

  # run IC
#  set jpref = $outdir/__jdata/job.`echo $mat | sed 's/\.[^.]\+$//'`
#  scripts-create-path $jpref
#  set jid = ($jid `scripts-qsub-run $jpref 1 $mem ./code/hicseq-matrix-cnv-melt.tcsh $outdir $params $branch $object $mat`)
  ./code/hicseq-matrix-cnv-melt.tcsh $outdir $params $branch $object $mat
end

# wait until all jobs are completed
scripts-send2err "Waiting until all jobs are completed..."
scripts-qsub-wait "$jid"

# hicapp normalization
scripts-send2err "Normalizing matrices..."
cat $outdir/matrix.chr*.tmp | sed 's/_/\t/2' | awk '$3!=0' | tr '_' '\t'  | awk '{print $1"_"$2-1"\t"$3"_"$4-1"\t"$5}' >! $outdir/matrix.pair
set p = `pwd`
cd $outdir
set tmp = `mktemp` 
cut -f1,3 $p/$genome_dir/genome.bed >! $tmp
$p/code/hicapp/hicapp_caICB --threads $threads -s matrix.pair -r $bin_size -c $tmp
awk 'NR==FNR {a[$1]=$2} NR>FNR {B=a[$1]*a[$2]; if (B>0) print $0"\t"$3/B}' matrix.bychrom.caicb matrix.pair | tr '_' '\t'  | awk '{print $1"_"$2+1"_"$3"_"$4+1"\t"$6}' >! matrix.pair.norm
cd $p

# unmelt the matrices
set matrices = `cd $outdir; ls -1 matrix.*.tmp | grep -vwE "$chrom_excluded"`
foreach mat ($matrices)
  scripts-send2err "Processing $mat..."
  
  # estimate required memory
  set mem = `./code/calc-matrix-memory.tcsh $outdir/$mat 1.5 5`
  scripts-send2err "requested memory = $mem"
  #set jpref = $outdir/__jdata/job.`echo $mat | sed 's/\.[^.]\+$//'`
  #scripts-create-path $jpref
  #set jid = ($jid `scripts-qsub-run $jpref 1 $mem ./code/hicseq-matrix-cnv-unmelt.tcsh $outdir $params $branch $object $mat`)
  ./code/hicseq-matrix-cnv-unmelt.tcsh $outdir $params $branch $object $mat
end

# wait until all jobs are completed
scripts-send2err "Waiting until all jobs are completed..."
scripts-qsub-wait "$jid"

# clean up
rm -rf $tmp $outdir/matrix.bychrom.* $outdir/matrix.pair.selected.* $outdir/matrix.pair $outdir/matrix.pair.norm

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


