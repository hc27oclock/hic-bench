#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: run-hicrep.tcsh OUTPUT-DIR MATRIX-1 MATRIX-2 BIN-SIZE MAX-DIST
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set mat1 = $2
set mat2 = $3
set bin_size = $4
set max_dist = $5

set max_dist_bite = `echo ${max_dist}\*${bin_size} | bc`

# create output directory
mkdir -p $outdir

# extract matrices from RData file
set tmpdir = $outdir
Rscript ./code/hic-matrix.r matrices -v --reverse-log2 -o $tmpdir/mat1 $mat1
Rscript ./code/hic-matrix.r matrices -v --reverse-log2 -o $tmpdir/mat2 $mat2

wc -l $tmpdir/mat?/*

# compare matrices
module unload r
module load r/3.3.0

rm -rf $outdir/compare.tsv
touch $outdir/compare.tsv
foreach f (`ls $tmpdir/mat1/* | xargs -n1 basename`)
  set k = `echo $f | tr "=." "\t" | cut -f3`
  set sample1 = `echo $outdir | xargs dirname | xargs basename | sed -E "s/\./\t/g" | cut -f1`
  set sample2 = `echo $outdir | xargs dirname | xargs basename | sed -E "s/\./\t/g" | cut -f2`
  set chr = `echo $outdir | xargs basename | sed -E "s/\./\t/g" | cut -f2`
  set cor = `Rscript code/run-hicrep.r -b $bin_size -s 0 -r inf -m $max_dist_bite  $tmpdir/mat1/$f $tmpdir/mat2/$f`
  echo "$sample1 $sample2 $chr $k $cor" | tr ' ' '\t' >> $outdir/compare.tsv
# echo "Rscript code/run-hicrep.r -b $bin_size -s 0 -r inf -m $max_dist_bite  $tmpdir/mat1/$f $tmpdir/mat2/$f"
end
