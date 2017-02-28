#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # custom shell environment

##
## USAGE: read-sample-sheet2.tcsh SAMPLE-SHEET-TSV OBJECT-NAME(S) FEATURE-NAME(S)
##

if ($#argv != 3) then
  grep '^##' $0
  exit 1
endif

set sheet = $1
set objects = ($2)
set features = ($3)

set features2 = `echo $features | tr ' ' '\n' | sed 's/^/^/' | sed 's/$/$/' | tr '\n' '|' | sed 's/|$//'`
set columns = `head -1 $sheet | tr '\t' '\n' | grep -nE "$features2" | cut -d':' -f1 | tr '\n' ',' | sed 's/,$//'`

if ("$objects" == "*") then
  cat $sheet | sed '1d' | cut -f$col
else
  set t = `scripts-create-temp`
  echo $objects | tr ' ' '\n' | sort -u >! $t
  cat $sheet | sed '1d' | cut -f$columns | sort | join -t '	' $t -
  rm -f $t
endif






