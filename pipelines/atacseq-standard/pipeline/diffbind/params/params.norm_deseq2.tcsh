#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload r
# module load r/3.2.3
module load r/3.3.0
module unload zlib

set diffbind_factor = 'group'
set diffbind_blocking_factor = ''
set norm_method = 'DESEQ2'


