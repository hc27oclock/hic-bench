#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: run-fastMeanFilter.tcsh PREP SMOOTH-PARAMS MATRIX-1 OUTPUT
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set prep = $1
set smoothing = "$2"
set mat1 = $3
set outdir = $4

module unload r
module load r/3.3.0

Rscript ./code/run-fastMeanFilter.r -s $smoothing -p $prep -o $outdir $mat1


