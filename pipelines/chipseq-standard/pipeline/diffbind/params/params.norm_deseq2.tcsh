#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload r
# module load r/3.2.3
module load r/3.3.0
module unload zlib

set diffbind_factor = 'group'
set diffbind_blocking_factor = ''
set norm_method = 'DESEQ2'
set external_table =                      # additional data (in TSV format) to be merged with diffbind results


