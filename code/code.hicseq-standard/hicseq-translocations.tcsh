#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-translocations.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

set filtered_reads = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/filtered.reg.gz"}'`
foreach r ($res)
  set y = $outdir/out.res_${r}kb
  cat $filtered_reads | gunzip | awk '$2!=$6' | gtools-regions sort | gtools-hic bin -v -g $genome_dir/genome.bed --bin-size ${r}000 | sort | uniq -c | awk -v c=$min_reads '$1>=c' | sort -k1,1rg >! $y.0
  cat $y.0 | sed 's/^ *//' | tr '\t:' ' ' | cut -d' ' -f1,2,4 | tools-cols 1 2 0 | sed 's/ /|/' | sort | tools-mergeuniq -merge >! $y.1
  cat $y.1 | tools-vectors n | paste - $y.1 | cut -f2- | sort -rn >! $y.2
end

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


