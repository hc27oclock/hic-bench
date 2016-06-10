#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload zlib
module unload r
# module load r/3.2.3
module load r/3.3.0

set diffbind_factor = 'group'
set diffbind_blocking_factor = ''


