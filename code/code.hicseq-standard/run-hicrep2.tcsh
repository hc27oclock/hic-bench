#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: run-hicrep.tcsh OUTPUT-DIR MATRIX-1 MATRIX-2 BIN-SIZE MAX-DISTANCES SMOOTH-PARAMS
##

if ($#argv != 6) then
  grep '^##' $0
  exit
endif

set outdir = $1
set mat1 = $2
set mat2 = $3
set bin_size = $4
set max_distances = ($5)
set smooth_params = ($6)

# create output directory
mkdir -p $outdir

# extract matrices from RData file
set ext = `echo $mat1 | sed 's/.*\.//'`
set tmpdir = $outdir      #TODO: set this to $TMP
if ($ext == 'tsv') then
  mkdir -p $tmpdir/mat1 $tmpdir/mat2
  cp $mat1 $tmpdir/mat1/`basename $mat1`
  cp $mat2 $tmpdir/mat2/`basename $mat2`
else
  Rscript ./code/hic-matrix.r matrices -v --reverse-log2 -o $tmpdir/mat1 $mat1
  Rscript ./code/hic-matrix.r matrices -v --reverse-log2 -o $tmpdir/mat2 $mat2
endif

# compare matrices
module unload r
module load r/3.3.0

rm -rf $outdir/compare.tsv
touch $outdir/compare.tsv
foreach max_dist ($max_distances)
  foreach smooth ($smooth_params)
    foreach f (`ls $tmpdir/mat1/* | xargs -n1 basename`)
      if (`echo $f | grep -c '='` == 0) then
        set k = 001
      else
        set k = `echo $f | tr "=." "\t" | cut -f3`
      endif
      set sample1 = `echo $outdir | xargs dirname | xargs basename | sed -E "s/\./\t/g" | cut -f1`
      set sample2 = `echo $outdir | xargs dirname | xargs basename | sed -E "s/\./\t/g" | cut -f2`
      set chr = `echo $outdir | sed 's/.*\.//'`
      set cor = `Rscript code/run-hicrep.r -b $bin_size -s $smooth -r inf -m $max_dist $tmpdir/mat1/$f $tmpdir/mat2/$f`
      echo "$sample1 $sample2 $chr $k $cor $max_dist $smooth" | tr ' ' '\t' >> $outdir/compare.tsv
    end
  end
end

rm -rf $tmpdir/mat1 $tmpdir/mat2

