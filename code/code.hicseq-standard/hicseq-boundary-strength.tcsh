#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-boundary-strength.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
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

# genome information (for TSS annotations)
set release = release     # TODO: $genome_dir???

# identify all input directories
set domains_dir = $branch
set matrix_branch = `echo $branch | sed 's/.*results\/domains\.[^/]\+.//'`
set matrix_task = `echo $matrix_branch | cut -d'.' -f1`
set matrix_dir = ../$matrix_task/results/$matrix_branch
set grouping = `echo $branch | sed 's/.*results\/domains\.//' | cut -d'.' -f1`
set bscores_dir = ../boundary-scores/results/boundary-scores.$grouping.$bscore_params/$matrix_branch

# filter boundaries near centromeres and telomeres!
if ($filter_centrotelo == true) then
  scripts-send2err "Filtering boundaries near centromeres and telomeres..."
  cat $genome_dir/centrotelo.bed | gtools-regions shiftp -5p -2000000 -3p +2000000 >! $outdir/filter.bed
else
  touch $outdir/filter.bed
endif
  
# analyze all samples
scripts-send2err "Classifying boundaries..."
echo -n "" >! $outdir/eval-SE.tsv
if ($n_boundaries == 0) then 
  set auto_nboundaries = true
else
  set auto_nboundaries = false
endif
foreach sample ($objects)
  scripts-send2err "$sample"

  # make output dir
  mkdir -p $outdir/$sample

  # obtain boundary scores    
  if (! -e $domains_dir/$sample/score.$optimal.bedgraph) then
    cat $bscores_dir/$sample/all_scores.$optimal.tsv | cut -f1,9 | grep -wv NA | scripts-skipn 1 | tr ':-' '\t' >! $outdir/scores.bedgraph    # TODO: this will use hic-ratio scores as default...
  else
    scripts-send2err "Using $domains_dir/$sample/score.$optimal.bedgraph..."
    cat $domains_dir/$sample/domains.$optimal.bed | gtools-regions inv -g $genome_dir/genome.bed | awk '{print $1,$2,$3,$1":"$2":"$3}' | gtools-overlaps overlap -i -label -t ' ' $domains_dir/$sample/score.$optimal.bedgraph | cut -f4,5 | sort | tools-mergeuniq -merge | tools-vectors min -n 6 | awk '{print $1,-1*$2}' | tr -s ': ' '\t' >! $outdir/scores.bedgraph
  endif

  # identify lost boundaries
  cat $outdir/scores.bedgraph | gtools-overlaps subset -i -inv $domains_dir/$sample/domains.$optimal.bed | gtools-overlaps subset -inv $outdir/filter.bed | sort -k4,4g >! $outdir/ranked-boundaries.bedgraph
  if ($auto_nboundaries == true) then
    set N = `cat $outdir/ranked-boundaries.bedgraph | wc -l`
    set n_boundaries = `echo $N/$n_categories | bc`
  endif
  scripts-send2err "N = $N, nboundaries = $n_boundaries"
  set N = `echo $n_categories\*$n_boundaries | bc`
  scripts-send2err "new-N = $N"
  cat $outdir/ranked-boundaries.bedgraph | tail -$N | gtools-regions shift -start -1 | gtools-regions shift -start -$flank_boundaries -stop +$flank_boundaries | split -l $n_boundaries - $outdir/$sample/blost.
  
  # clean-up
  rm -f $outdir/scores.bedgraph
  
  # rename files
  set k = 1
  foreach b ($outdir/$sample/blost.*)
    set bnew = `echo $b | sed 's/\.[a-z]\+$'"/.0$k.bed/"`
    mv $b $bnew
    @ k ++
  end

  # evaluate proximity to super-enhancers
  set p = `pwd`
  cd $outdir/$sample
  ln -s $p/code
  set cell_type = `echo $sample | cut -d'-' -f1`
  ( grep . blost.0* | cut -f-3 | tr ':' '\t' | sed 's/^blost.//' | tools-cols -t 1 2 3 0 | sed 's/.bed$//' ; cat $p/superenhancer/superenhancer/$cell_type-H3K27ac.bed | sed 's/$/\tSE/' ) | scripts-sortbed >! blost-SE.bed
  gtools-regions gdist blost-SE.bed | grep SE | cut -f-2 | tr '\t' '\n' | sort | uniq -c | grep -v SE | awk '{print $2,$1}' | sort >! b1
  cat blost-SE.bed | grep -v SE | cut -f4 | sort | uniq -c | awk '{print $2,$1}' | sort | join - b1 | sed 's/ /\t/' | tools-vectors -div -n 6 | grep -v KBM7-untreated-MboI | awk -v var=$sample '{print var"\t"$0}' >> $p/$outdir/eval-SE.tsv
  cd $p
end

# overlap with CTCF
foreach filter (noncoding tss all)
  scripts-send2err "Evaluating CTCF, filtered by $filter..."
  ./code/hicseq-boundary-strength-eval.tcsh $outdir $release CTCF $filter | perl -F/\\t/ -lane '@a=split(" ",$F[1]);$sum=0;foreach $ele (@a){$sum=$sum+$ele};$mean=$sum/@a;for($i=1;$i<=@a;$i++){$final=$a[$i-1]/$mean;print "$F[0]\tb$i\t$final"}' >! $outdir/eval-CTCF-$filter.tsv
end

# overlap with repeat elements
#./code/hicseq-boundary-strength-eval.tcsh $outdir $genome_dir repeats >! $outdir/eval-repeats.tsv

# create matrix for heatmap of boundary strength categories across all samples
set p = `pwd`
cd $outdir
ln -s $p/code
cat */blost.* | scripts-sortbed | gtools-regions link -d 100000 >! x1
cat x1 | gtools-regions n | paste - x1 | awk '$2>=40000' | cut -f3- >! x11
grep . */blost.* | tr ':' '\t' | tools-cols -t 1 2 3 0 | gtools-regions center | gtools-regions pos -op 5p >! x2
cat x11 | gtools-overlaps overlap -i -label -t '|' x2 | sort -u | sed 's/_|//' | cut -f-4 | sed 's/\/blost./\t/' | sed 's/.bed$//' | sort -k1,1 -k2,2g -k3,3 | sed 's/\t/:/' | sed 's/\t/:/' | sed 's/\t/|/' | tools-mergeuniq -merge | tools-vectors min -n 0 | tr '|' '\t' | tools-table -c -n 0 >! bscore-matrix.tsv
(printf "locus"; cat bscore-matrix.tsv | tr ' ' '\t' | sed -E 's/\t$//' | sed -E 's/-untreated//g' | sed 's/-DMSO//g' | sed 's/rep//g' | cut -f`head -1 bscore-matrix.tsv | tr ' \t' '\n' | cat -n | grep -v "NcoI\|DpnII" | cut -f1 | tr ' \t\n' ' ' | sed -E 's/ +/,/g' | sed -E 's/^,|,$//g'`) >! x1
head -1 x1 >! bscore-matrix.tsv
cat x1 | tr '\t' ' ' | tr ':' '\t' | tr ' ' ':' | sed 's/:/\t/' | scripts-skipn 1 | gtools-overlaps subset -inv $p/$genome_dir/centrotelo.bed | tr '\t' ' ' | tr ':' '\t' | tr ' ' ':' | sed 's/:/\t/3' >> bscore-matrix.tsv
rm -f x1 x2 x11

Rscript code/hicseq-boundary-strength-publication.r `ls eval-CTCF-* | tr '\n' ',' | sed 's/,$//g'` eval-SE.tsv bscore-matrix.tsv
rm -rf Rplots.pdf

cd $p

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


