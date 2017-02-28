#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-virtual-4c.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# check number of objects
if ($#objects > 1) then
  send2err "Error: hicseq-virtual-4c currently allows only 1 input object."
  exit 1
endif

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# Create the working directory
if ($?TMP) then
  set tempdir = $TMP
else
  set tempdir = $outdir
endif
set workdir = $tempdir/work
mkdir -p $workdir

# create input matrices
set object1 = $objects[1]
set chrom_included = `cat $loci | cut -f1 | sort -u`
foreach chr ($chrom_included)
  scripts-send2err "Processing chromosome $chr..."
  
  # number of reads
  set scaling = 1.0     # default
  if ($scale == yes) then
    if (-e $branch/$object1/stats.tsv) then
      set n_reads = `cat $branch/$object1/stats.tsv | grep '^ds-accepted-intra	' | cut -f2`
      set scaling = `echo $n_reads/1000000000 | bc -l`
    endif
  endif
  scripts-send2err "- scaling factor = $scaling"
  
  # extract matrices
  if (-e $branch/$object1/matrix.$chr.RData) then
    Rscript ./code/hic-matrix.r matrices -v -o $workdir/tmp $branch/$object1/matrix.$chr.RData
  else if (-e $branch/$object1/matrix.$chr.tsv) then
    mkdir -p $workdir/tmp
    cat $branch/$object1/matrix.$chr.tsv >! $workdir/tmp/matrix.k=001.tsv
  endif
  
  # process matrices
  foreach k (`cd $workdir/tmp; ls -1 matrix.*.tsv | cut -d'.' -f2`)
    set f = $workdir/tmp/matrix.$k.tsv
    if (`cat $f | head -1 | grep 'chr.*chr' | wc -l` == 1) then
      set fnew = $f																											# "full" matrix
    else
      set fnew = $workdir/matrix.tsv
      Rscript ./code/hic-matrix.r rotate -v --row-labels -o $fnew $f		# inverse-rotate matrix
    endif
    set baits = `cut -f1 $fnew | scripts-skipn 1 | tools-cols -t 0 0 | sed 's/:/\t/' | sed 's/-/\t/' | gtools-overlaps overlap -i -label -t '|' $loci | cut -f4`
    foreach b ($baits)
      set blocus = `echo $b | cut -d'|' -f1`
      set bname = `echo $b | cut -d'|' -f2 | tr ':' '_'`
      set col = `cat $fnew | grep -n "^$blocus	" | cut -d':' -f1`
      cat $fnew | cut -f1,$col | scripts-skipn 1 | awk '$2>0' | sed 's/$/	'$scaling'/' | awk '{print $1"\t"$2/$3}' | sed 's/:/\t/' | sed 's/-/\t/' >! $outdir/$bname.bedgraph
    end
  end
  
end



# Cleanup
if ($tempdir == $outdir) then
  rm -rf $workdir
endif
	
# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


