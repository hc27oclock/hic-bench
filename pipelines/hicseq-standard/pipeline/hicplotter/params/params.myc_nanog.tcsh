#!/bin/tcsh

source ./inputs/params/params.tcsh

# HiCplotter path
set hicplotter_path = ./code/HiCPlotter2.py

# create bedgraphs for boundary scores
set bscores_branch = ../boundary-scores/results/boundary-scores.by_sample.prep_none/`echo $branch | sed 's/.*results\///'`
set cell_type = `./code/code.main/read-sample-sheet.tcsh $sheet "$objects" cell_type`
set methods = (inter DI ratio)
set kappas = `cd $bscores_branch/$objects[1]; ls -1 all_scores.k=*.tsv | cut -d'.' -f2`
set bedgraphs = ()
set bedgraph_labels = ()
foreach kappa ($kappas)
  set f = $bscores_branch/$objects[1]/all_scores.$kappa.tsv
  set bfiles = ()
  set blabels = ()
  
  # add bscores for each method
  foreach m ($methods)
    set k = `head -1 $f | tr '\t' '\n' | grep -n "^$m"'$' | cut -d':' -f1`             # get column number for method
    cat $f | sed '1d' | cut -f1,$k | sed 's/:/\t/' | sed 's/-/\t/' | grep -v '	NA$' >! $outdir/bscores.$kappa.$m.bedGraph               # convert to bedgraph, remove NA values
    set bfiles = ($bfiles ../bscores.$kappa.$m.bedGraph)
    set blabels = ($blabels $m)
  end

  # add CTCF ChIP-seq
  if (-e inputs/data.external/$cell_type/CTCF.bedGraph) then
    set bfiles = ($bfiles `readlink -f inputs/data.external/$cell_type/CTCF.bedGraph`)
    set blabels = ($blabels CTCF)
  endif

  # add to the lists
  set bedgraphs = ($bedgraphs `echo $bfiles | tr ' ' ','`)
  set bedgraph_labels = ($bedgraph_labels `echo $blabels | tr ' ' ','`)
end

# regions to plot
set regions = (chr8:125000000-133000000 chr12:3940000-11950000)
set tiles = "params/regions.bed"
set tiles_labels = "regions"
set highlight = 0
set highlight_bed = ""
set fileheader = 0         # Either 1 or 0 (header / no header)
set insulation_score = 0   # Either 1 or 0 (include insulation index or not)



