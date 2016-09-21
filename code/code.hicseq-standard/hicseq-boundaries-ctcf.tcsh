#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hic-boundaries-ctcf.tcsh OUTPUT-DIR PARAM-SCRIPT ALIGNMENT-BRANCH OBJECT(S)
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

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# cell type
set celltype = `./code/read-sample-sheet.tcsh $sheet "$objects" celltype`

# Concatenate blacklisted regions
cat $black_lists | gtools-regions bed >! $outdir/blacklist.bed

# boundary-CTCF overlaps
set domaindir = $branch/$object
set pref = domains
set kappas = `ls -1 $domaindir/domains.k=*.bed | sed "s/.*$pref\.//" | sed 's/\.bed$//' | sort -u`
foreach k ($kappas)
  set domains = $domaindir/$pref.$k.bed
  set boundaries = $outdir/boundaries.$k.bed
  ( gtools-regions pos -op 5p $domains ; gtools-regions pos -op 3p $domains ) | gtools-regions shiftp -5p -$flank_dist -3p +$flank_dist | gtools-overlaps subset -inv $outdir/blacklist.bed | scripts-sortbed >! $boundaries
  cat $peaks_dir/$celltype-CTCF/peak-scores.bed | gtools-overlaps subset -i $boundaries >! $outdir/peaks.bed
end

# save variables
set >! $outdir/job.vars.tsv

scripts-send2err "Done."




