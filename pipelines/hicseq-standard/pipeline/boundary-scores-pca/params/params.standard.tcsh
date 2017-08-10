#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MYX]'              # excluded chromosomes

set group_var = 'cell-type'                      # grouping variable (from sample sheet) to be used for color assignment)

set pca_params = '--show-text --use-short-names --plain'

