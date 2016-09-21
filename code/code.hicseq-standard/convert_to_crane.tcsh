#!/bin/tcsh

##
## USAGE: convert_to_crane.tcsh INPUT-MATRIX OUTPUT-MATRIX 
##

if ($#argv != 2) then
  grep '^##' $0
  exit
endif

set input_matrix = $1
set output_matrix = $2

# Convert
cat $input_matrix | sed 1d | sed 's/:/	/' | sed 's/-/	/' | awk -F"\t" 'BEGIN {OFS = FS} {$3=$3+1; print}' | sed 's/	/:/' | sed 's/	/-/' | sed 's/^/hg19\|/' | awk '{print "bin"NR"|"$0}'>! $input_matrix.x1.tmp
cut -f1 $input_matrix.x1.tmp | tr '\n' '	' | sed 's/^/	/' >! $input_matrix.x2.tmp
awk '{print}' $input_matrix.x2.tmp $input_matrix.x1.tmp >! $output_matrix
rm -rf $input_matrix.*.tmp
