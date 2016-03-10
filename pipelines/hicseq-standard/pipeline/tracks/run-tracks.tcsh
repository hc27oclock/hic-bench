#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-tracks.tcsh [--dry-run]
##

#% This step generates Hi-C tracks for the WashU Genome Browser. 
#% TABLES: 
#% FIGURES: 

# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit
endif

set opt = "$1"

# setup
set op = tracks
set inpdirs = "inpdirs/*"
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 4,20G
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/hicseq-$op.tcsh"

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "group" 1

# run
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run



