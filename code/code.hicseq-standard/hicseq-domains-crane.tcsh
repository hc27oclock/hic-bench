#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-domains-crane.tcsh OUTDIR PARAMS BRANCH OBJECTS
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

if ($#objects != 1) then
  scripts-send2err "Error: hicseq-domains-crane.tcsh requires exactly one input object."
  exit 1
endif

# create ignored loci file
./code/create-ignored-loci.tcsh $genome_dir $bin_size >! $outdir/ignored_loci.txt

# call domains for each chromosome
set inpdir = $branch/$objects[1]
set workdir = $outdir/work
mkdir -p $workdir

set est_matrices = `cd $inpdir; ls -1 matrix.*.tsv matrix.*.RData | grep -vwE "$chrom_excluded"`
foreach est_mat ($est_matrices)
  scripts-send2err "Processing matrix $est_mat..."
  set chr = `echo $est_mat | cut -d'.' -f2`

  # extract matrices from RData file
  if (`echo $est_mat | grep -c '\.RData$'` == 1) then
    Rscript ./code/hic-matrix.r matrices -v $hicmatrix_params -o $workdir/tmp $inpdir/$est_mat
  else
    mkdir -p $workdir/tmp
    cat $inpdir/$est_mat >! $workdir/tmp/matrix.k=001.tsv
  endif
  
  # convert matrices and run crane
  foreach mat ($workdir/tmp/matrix.*.tsv)
    set pref = `basename $mat .tsv | sed 's/^matrix\.//'`.$chr
    set inpmat = $pref.matrix.txt
    ./code/convert_to_crane.tcsh $mat $workdir/$inpmat
    set p = `pwd`
    set cranepath_abs = `readlink -f $cranepath`
    cd $workdir
    perl $cranepath_abs/matrix2insulation.pl -i $inpmat -is $inssqr -ids $idspan -im $insmode -nt $noise_thr -bmoe $bmoerr -v
    rm -f $inpmat
    cd $p
  end
  rm -rf $workdir/tmp
end

# Collect domains from all chromosomes
set K = `ls -1 $outdir/*/*k=*.bed | sed "s/.*k=//" | cut -d'.' -f1 | sort -u`
foreach k ($K)
  cat $outdir/*/k=$k.*.bed | grep -v track | gtools-regions shift -start -1 -stop -1 | gtools-regions inv -g $genome_dir/genome.bed | cut -f-3 | sed 's/^chr//g' | sed 's/^X	/23	/g' | sort -k1,1n -k2,2n | sed 's/^23	/X	/g' | sed 's/^/chr/g' >! $outdir/domains.k=$k.bed
end

# remove unused files
rm -rf $workdir

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
