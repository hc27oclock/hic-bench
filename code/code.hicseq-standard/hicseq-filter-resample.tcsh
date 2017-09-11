#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hic-filter-resample.tcsh OUTPUT-DIR PARAM-SCRIPT ALIGNMENT-BRANCH OBJECT(S)
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir enzyme"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# resample
set stats_files = `echo $objects | tr ' ' '\n' | awk -v b=$branch '{print b"/"$1"/stats.tsv"}'`
set n_total = `cat $stats_files | grep '^ds-accepted-intra	' | cut -f2 | tools-matrix csum -ll`
scripts-send2err "total reads = $n_total"
Rscript ./code/create-random-vector.r $n_reads $n_total 1 >! $outdir/v
set reg_files = `echo $objects | tr ' ' '\n' | awk -v b=$branch '{print b"/"$1"/filtered.reg.gz"}'`
cat $reg_files | gunzip | awk '$2==$6' | paste $outdir/v - | grep '^1' | cut -f2- | gzip >! $outdir/filtered.reg.gz

# create stats file
set n_used_reads = `cat $outdir/v | grep ^1 | wc -l`
cat $branch/$objects[1]/stats.tsv | cut -f1 | sed 's/$/ 0 0/' | sed "s/^read-pairs .*/read-pairs $n_used_reads 1/" | sed "s/^ds-accepted-intra .*/ds-accepted-intra $n_used_reads 1/" | tr ' ' '\t' >! $outdir/stats.tsv
rm -f $outdir/v

# save variables
source ./code/code.main/scripts-save-job-vars

scripts-send2err "Done."




