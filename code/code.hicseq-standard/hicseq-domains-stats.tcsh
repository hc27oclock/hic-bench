#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-domain-stats.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
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
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------
# Isolate all the kappas
set kappas = `ls -1 $branch/*/domains.*.bed | sed 's/.*\///' | cut -d'.' -f2 | sort -u`

# Call the Rscript that plots the domain numbers and
# domain distribution for each one of the kappas (lambdas)
# separetely
foreach kappa ($kappas)
	echo $kappa
        Rscript ./code/plot-domains-stats.r $outdir $branch $kappa
end

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


