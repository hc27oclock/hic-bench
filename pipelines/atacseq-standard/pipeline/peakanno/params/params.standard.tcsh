#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload r
module load r/3.3.0

set promoter_proximal = 3000      # Extending promoter upstream and downstream by nt
set include_input = 'false'       # include ChIP inputs
