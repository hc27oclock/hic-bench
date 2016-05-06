#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-matrix-pca.tcsh OUTDIR PARAM-SCRIPT BOUNDARY-SCORES-BRANCH [OBJECTS]
##

if (($#argv < 3) || ($#argv > 4)) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# if objects is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# Check for number of objects
if ($#objects < 2) then
	scripts-send2err "Error: more than one input objects are required."
	exit 1
endif 

# Generate PCA plots
./code/read-sample-sheet.tcsh $sheet "$objects" $group_var yes | awk '{print $2":"$1}' >! $outdir/labels.tsv
set matrices = `cd $branch/$objects[1]; ls -1 matrix.*.tsv | grep -vwE "$chrom_excluded"`
foreach object ($objects)
  scripts-send2err "Processing $object..."
  echo -n "" >! $outdir/data.$object.tsv
  foreach mat ($matrices)
    cat $branch/$object/$mat | scripts-skipn 1 | cut -f2- | tr -s '\t' '\n' >> $outdir/data.$object.tsv
  end
end
paste $outdir/data.*.tsv >! $outdir/matrix.tsv
rm -f $outdir/data.*.tsv
set metric = max
set nbest = 100000
cat $outdir/matrix.tsv | tr '\t' ' ' | tools-vectors $metric >! $outdir/matrix.$metric.tsv
set cutoff = `cat $outdir/matrix.$metric.tsv | sort -rg | head -$nbest | tail -1`
echo "$objects" | tr ' ' '\t' >! $outdir/matrix.filtered.tsv
cat $outdir/matrix.$metric.tsv | tools-vectors test -g -e -c $cutoff | paste - $outdir/matrix.tsv | grep ^1 | cut -f2- | tools-rows -number >> $outdir/matrix.filtered.tsv
scripts-perform-pca.r -v -o $outdir -L $outdir/labels.tsv --show-text --use-short-names $outdir/matrix.filtered.tsv
rm -f $outdir/matrix.tsv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




