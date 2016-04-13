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

# test number of input objects
set object = $objects[1]
if ($#objects != 1) then
  scripts-send2err "Error: this operation allows only one input object!"
  exit 1
endif

# run parameter script
source $params

# indentify genome directory
set genome = `./code/read-sample-sheet.tcsh $sheet "$objects" genome`
set enzyme = `./code/read-sample-sheet.tcsh $sheet "$objects" enzyme`
set genome_dir = inputs/genomes/$genome

# create path
scripts-create-path $outdir/

# resample
set n_total = `cat $branch/$object/stats.tsv | grep '^ds-accepted-intra	' | cut -f2`
Rscript ./code/create-random-vector.r $n_reads $n_total 1 >! $outdir/v
cat $branch/$object/filtered.reg.gz | gunzip | awk '$2==$6' | paste $outdir/v - | grep '^1' | cut -f2- | gzip >! $outdir/filtered.reg.gz

# create stats file
set n_used_reads = `cat $outdir/v | grep ^1 | wc -l`
cat $branch/$object/stats.tsv | cut -f1 | sed 's/$/ 0 0/' | sed "s/^read-pairs .*/read-pairs $n_used_reads 1/" | sed "s/^ds-accepted-intra .*/ds-accepted-intra $n_used_reads 1/" | tr ' ' '\t' >! $outdir/stats.tsv
rm -f $outdir/v

# save variables
set >! $outdir/job.vars.tsv

scripts-send2err "Done."




