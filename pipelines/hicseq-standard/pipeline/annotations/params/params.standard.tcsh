#!/bin/tcsh

source ./inputs/params/params.tcsh

set genes_bed = $genome_dir/gene.bed                         # gene BED6 file for annotation of interactions
set cell_type = `echo $objects[1] | cut -d'-' -f1 | cut -d'_' -f1`
if (! -e inputs/data.external/$cell_type) then
  set loci_bed = ()
else
  set loci_bed = `find inputs/data.external/$cell_type -maxdepth 1 -name '*.bed'`
endif


