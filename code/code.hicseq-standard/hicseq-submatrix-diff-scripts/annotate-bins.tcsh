#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: annotate-bins.tcsh OUTPUT-DIR INPUT-TABLE GENOME-DIR
##

if ($#argv != 3) then
  grep '^##' $0
  exit
endif

set outdir = $1
set inp = $2
set genome_dir = $3

# generated annotated bed file
set t = `scripts-create-temp $outdir`
cat $genome_dir/gene-name.bed | gtools-regions reg | gtools-regions pos -op 5p | sed 's/^/tss:/' >! $t.annot.bed

# annotate bins
set columns = `head -1 $inp | tr '\t' '\n' | grep -n '^locus[12]$' | cut -d':' -f1`
set columns = `echo $columns | tr ' ' ','` 
cat $inp | scripts-skipn 1 | cut -f$columns | tr '\t' '\n' | sort -u | tools-cols -t 0 0 | sed 's/ /_/' | sed 's/ /_/' | sed 's/ /_/' | gtools-overlaps overlap -i -label $t.annot.bed | cut -f1 | sed 's/:/\t/' | sort | tools-mergeuniq -merge -t ',' >! $t.x2

# prepare output table for merging
cat $inp | cut -f$columns | tr ' ' '_' | paste - $inp | scripts-skipn 1 | sort >! $t.x1

# merge bin pairs with annotations
cat $t.x1 | cut -f-2 | join -t '	' -a1 -e 'N/A' -o 1.1 1.2 2.2 - $t.x2 | sort -k2,2 | join -t '	' -a1 -e 'N/A' -1 2 -o 1.1 1.2 1.3 2.2 - $t.x2 | sed 's/\t/|/' | sort >! $t.x3

# generate final annotated table
echo `head -1 $inp` locus1-annotations locus2-annotations | tr ' ' '\t'
cat $t.x1 | sed 's/\t/|/' | sort | join -t'	' - $t.x3 | cut -f2- | sort -k2,2g 

# cleanup
rm -f ${t}.* $t

# done
scripts-send2err "Done."


