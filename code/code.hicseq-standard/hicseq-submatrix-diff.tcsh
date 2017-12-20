#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-submatrix-diff.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT1 OBJECT2
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

set objects = ($object1 $object2)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# code dir
set code_dir = code/hicseq-submatrix-diff-scripts

# perform Hi-C submatrix differential analysis
Rscript $code_dir/find-hic-squares.r -v -o $outdir --bin-size $bin_size $submatrix_diff_params $branch/$object1 $branch/$object2

# generated annotated bed file
set t = $outdir/tmp
mkdir -p $t
cat $genome_dir/gene-name.bed | gtools-regions reg | gtools-regions pos -op 5p | sed 's/^/tss:/' >! $t/annot.bed
#cat $genome_dir/gene-name.bed | gtools-regions reg | sed 's/^/gbody:/' >> $t/annot.bed
set columns = `head -1 $outdir/out.tsv | tr '\t' '\n' | grep -n '^locus[12]$' | cut -d':' -f1`
cat $outdir/out.tsv | scripts-skipn 1 | cut -f`echo $columns | tr ' ' ','` | tr '\t' '\n' | sort -u | sed 's/^/_\t/' | gtools-overlaps overlap -i -label $t/annot.bed | cut -d':' -f2- | tools-cols -t 1 0 | tr ' ' '_' | sort | tools-mergeuniq -merge -t ',' >! $t/bins.tsv
cat $outdir/out.tsv | scripts-skipn 1 | cut -f`echo $columns | tr ' ' ','` | tr '\t' '\n' | tr ' ' '_' | sort -u | join -t'	' -a1 -e 'N/A' -o 1.1 2.2 - $t/bins.tsv | sed 's/_/ /' | sed 's/_/ /' | sed 's/_/ /' | tools-cols -t 1 0 | gtools-regions bed | cut -f-4 | sort -k1,1 -k2,2n >! $t/bins-annotated.bed

# generate annotated table
cat $outdir/out.tsv | scripts-skipn 1 | tr '\t' ' ' | sed 's/ /|/' | sed 's/ /|/' | sed 's/ /|/' | sed 's/ /|/' | sed 's/ /\t/' | gtools-overlaps bin -v -i --print-labels --print-regions $t/bins-annotated.bed | tr '|' '\t' | sort -k2,2g >! $outdir/results.tsv

# cleanup
rm -rf $t

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


