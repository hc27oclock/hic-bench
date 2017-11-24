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

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# convert to IC input format
scripts-send2err "Converting to unmelted matrix format..."
set input = $outdir/$mat
set output = `echo $input | sed 's/.tmp//g'`
set chr = `basename $input | tr '.' '\t' | cut -f2 | awk '{print $1"_"}'`
set tmp = `mktemp`
set tmp1 = `mktemp`

grep $chr $outdir/matrix.pair.norm | sort -k1,1 >! $tmp
sort -k1,1 $input >! $tmp1

join -t '	' -e 0 -o 0,2.2 -1 1 -2 1 -a1 $tmp1 $tmp | sed 's/_/\t/g' | sort -k2,2n -k4,4n | awk -v f=$bin_size '{print $1":"$2"-"$2+f-1"\t"$3":"$4"-"$4+f-1"\t"$5}' | sed 's/:/\t/' | sed 's/-/\t/' | sed 's/\t/@/4' | gtools-regions bounds -g $genome_dir/genome.bed | sed 's/@/\t/' | sed 's/\t/:/' | sed 's/\t/-/' | tools-cols -t 1 0 2 | sed 's/:/\t/' | sed 's/-/\t/' | sed 's/\t/@/4' | gtools-regions bounds -g $genome_dir/genome.bed | sed 's/@/\t/' | sed 's/\t/:/' | sed 's/\t/-/' | tools-cols -t 1 0 2 | tools-table -c | tr ' ' '\t' | sed 's/\t$//'  >! $output


# clean up
rm -rf $input $tmp $tmp1

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


