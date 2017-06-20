#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload r
module unload java
module load r/3.3.0
module unload zlib

set diffbind_factor = 'group'
set diffbind_blocking_factor = ''


