#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-ic-mirny.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT(S) MATRIX
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)          # TODO: allow multiple samples, so that this operation can be run by-group
set mat = $5

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

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# convert to IC input format
scripts-send2err "Converting to melted matrix format..."
set inpdir = $branch/$object
set x = $outdir/$mat.tmp
  
tools-table -x0 $inpdir/$mat | tr ':-' '\t' | cut -f1,2,4,5,7 | sed 's/\t/_/' | sed 's/\t/_/' | sed 's/\t/_/' >! $x


# this part is to do hicapp CNV correction individually for different chromosome, the result is a little bit different than doing it all together as default practice.
#scripts-send2err "Normalizing matrices..."
#cat $x | sed 's/_/\t/2' | awk '$3!=0' | tr '_' '\t'  | awk '{print $1"_"$2-1"\t"$3"_"$4-1"\t"$5}' >! $x.pair
#set p = `pwd`
#cd $outdir
#set y = `basename $x`
#set tmp = `mktemp`
#cut -f1,3 $p/$genome_dir/genome.bed >! $tmp
#$p/code/hicapp/hicapp_caICB --threads $threads -s $y.pair -r $bin_size -c $tmp
#awk 'NR==FNR {a[$1]=$2} NR>FNR {B=a[$1]*a[$2]; if (B>0) print $0"\t"$3/B}' $y.bychrom.caicb $y.pair | tr '_' '\t'  | awk '{print $1"_"$2+1"_"$3"_"$4+1"\t"$6}' >! $y.pair.norm
#cd $p
#join -t '	' --nocheck-order -e 0 -o 0,2.2 -1 1 -2 1 -a1 $x $x.pair.norm | sed 's/_/\t/g' | awk -v f=$bin_size '{print $1":"$2"-"$2+f-1"\t"$3":"$4"-"$4+f-1"\t"$5}' | sed 's/:/\t/' | sed 's/-/\t/' | gtools-regions bounds -g $genome_dir/genome.bed | sed 's/\t/:/' | sed 's/\t/-/' | tools-cols -t 1 0 2 | sed 's/:/\t/' | sed 's/-/\t/' | gtools-regions bounds -g $genome_dir/genome.bed | sed 's/\t/:/' | sed 's/\t/-/' | tools-cols -t 1 0 2 | tools-table -c | tr ' ' '\t' | sed 's/\t$//'  >! $x.tsv
# clean up
#rm -rf $tmp $x.bychrom.* $x.pair.selected.* $x.pair $x.pair.norm


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


