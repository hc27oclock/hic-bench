#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-integrate-icgc.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
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

# collect boundary strength data
(cd $branch; grep . all-samples/*/blost.*.bed) | sed 's/:/\t/' | tools-cols -t 1 2 3 0 | scripts-sortbed >! $outdir/boundary_strength_data_sorted.bed

# setup
set release = release     # TODO: $genome_dir???
set mpath = `pwd`
cd $outdir
ln -s $mpath/code
ln -s $mpath/$release
ln -s $mpath/data.local/ICGC/unique_SEs_from_SEA_sorted.bed
ln -s $mpath/data.local/ICGC/unique_enhancers_from_FANTOM_sorted.bed
ln -s $mpath/data.local/ICGC/ICGC_deletions.bed
ln -s $mpath/data.local/ICGC/ICGC_duplications.bed
ln -s $mpath/data.local/human-protein-atlas/gene-class.tsv
set log = log
echo -n '' >! $log

# gene classes
cat release/gene-name.bed | gtools-regions reg | sort | join -t'	' - release/gene-name.description | grep protein_coding: | cut -f-2 | gtools-regions bed | cut -f-6 | scripts-sortbed -i >! protein_coding.bed
cat release/gene-name.bed | gtools-regions reg | sort | join -t'	' - release/gene-name.description | grep -E '	protein_coding:|	lincRNA:|	miRNA:|	rRNA:' | cut -d':' -f1 | tools-cols -t 2 0 1 | sed 's/\t/:/' | gtools-regions bed | cut -f-6 | scripts-sortbed -i >! gene-name-biotype.bed
#cat release/gene-name.bed | gtools-regions reg | sort | join -t'	' - release/gene-name.description | cut -d':' -f1 | cols -t 2 0 1 | sed 's/\t/:/' | gtools-regions bed | cut -f-6 | scripts-sortbed -i >! gene-name-biotype.bed
cat release/gene-name.bed | gtools-regions reg | sort | join -t'	' - gene-class.tsv | tools-cols -t 2 0 1 | sed 's/\t/:/' | gtools-regions bed | cut -f-6 | scripts-sortbed -i >! gene-name-class.bed
 
foreach t (duplications deletions)
  cat ICGC_${t}.bed | awk -v max=$max_size '$3-$2<=max' | tools-rows -number | awk '{print $2"\t"$3"\t"$4"\t"$1"-"$5}' >! ICGC_${t}_uniqid.bed
  cat ICGC_${t}_uniqid.bed | gtools-regions inv -g release/genome.bed | cut -f-3 | tools-rows -number | tools-cols -t 1 2 3 0 >! ICGC_${t}_uniqid_inv.bed
end

# Find number of duplications and deletions
set nodups = `cat ICGC_duplications_uniqid.bed | wc -l`
set nodels = `cat ICGC_deletions_uniqid.bed | wc -l`

# find closest boundary for each SE (ignore overlaps)
bedtools closest -a unique_SEs_from_SEA_sorted.bed -b boundary_strength_data_sorted.bed -t all -io >! SE_closest_boundaries.tsv
cat SE_closest_boundaries.tsv | awk '{print $8"|"$4"\t"$1" +",$2,$3,$5" +",$6,$7}' | cut -d'/' -f2- | gtools-regions connect | gtools-regions bed | cut -f-4 >! SE_closest_boundaries.bed
echo "SEs and closest boundaries:" >> $log
cat SE_closest_boundaries.bed | cut -f4 | cut -d'|' -f1 | cut -d'/' -f2 | tools-countuniq -s | cut -f2 | sed 's/^/SE\t/' | tools-mergeuniq -merge | tools-vectors div -n 3 >> $log
echo >> $log

# find closest boundary for regular enhancer (ignore overlaps)
bedtools closest -a unique_enhancers_from_FANTOM_sorted.bed -b boundary_strength_data_sorted.bed -t all -io >! enhancer_closest_boundaries.tsv
cat enhancer_closest_boundaries.tsv | awk '{print $8"|"$4"\t"$1" +",$2,$3,$5" +",$6,$7}' | cut -d'/' -f2- | gtools-regions connect | gtools-regions bed | cut -f-4 >! enhancer_closest_boundaries.bed
echo "Enhancers and closest boundaries:" >> $log
cat enhancer_closest_boundaries.bed | cut -f4 | cut -d'|' -f1 | cut -d'/' -f2 | tools-countuniq -s | cut -f2 | sed 's/^/enhancer\t/' | tools-mergeuniq -merge | tools-vectors div -n 3 >> $log
echo >> $log

# find closest genes/boundaries
#bedtools closest -a boundary_strength_data_sorted.bed -b gene-name-biotype.bed -t all >! boundary_closest_gene.tsv
#echo "Boundaries and closest genes:" >> $log
#cat boundary_closest_gene.tsv | grep blost | cut -f4,8 | cut -d'/' -f3- | cut -d':' -f1 | tools-cols -t 1 0 | tools-countuniq -s | cut -f1,3 | tools-mergeuniq -merge | vectors div -n 3 >> $log
bedtools closest -a gene-name-biotype.bed -b boundary_strength_data_sorted.bed -t all >! gene_closest_boundary.tsv
echo "Genes and closest boundaries:" >> $log
cat gene_closest_boundary.tsv | grep blost | cut -f4,10 | tools-cols -t 1 0 | cut -d'/' -f3- | cut -d':' -f1 | tools-cols -t 1 0 | tools-countuniq -s | cut -f1,3 | tools-mergeuniq -merge | tools-vectors div -n 3 >> $log
echo >> $log
bedtools closest -b gene-name-class.bed -a boundary_strength_data_sorted.bed -t all >! boundary_closest_gene_class.tsv
echo "Genes and closest boundaries:" >> $log
cat boundary_closest_gene_class.tsv | grep blost | cut -f4,8 | cut -d'/' -f3- | cut -d':' -f1 | tools-cols -t 1 0 | tools-countuniq -s | cut -f1,3 | tools-mergeuniq -merge | tools-vectors div -n 3 >> $log
echo >> $log

foreach t (duplications deletions)
#  echo "${t} and boundaries"
#  cat boundary_strength_data_sorted.bed | gtools-overlaps overlap -label ICGC_${t}_uniqid.bed | cut -f4 | cut -d':' -f1 | cut -d'/' -f3 | sort | uniq -c

  echo "${t} and closest boundaries (inside/outside):" >> $log
  bedtools closest -a ICGC_${t}_uniqid.bed -b boundary_strength_data_sorted.bed -t all >! closest-outside-${t}.tsv
  cat closest-outside-${t}.tsv | cut -f8 | grep blost | cut -d'/' -f2- | tr '/' '\t' | tools-countuniq -s | sort | tools-table -c -n 0 | scripts-skipn 1 | sed 's/ $//' | tr -s ' ' '\t' >! closest-outside-${t}.boxplot.tsv
  cat closest-outside-${t}.boxplot.tsv | tools-matrix cstat | grep -E '^AVG|STD' | tools-matrix T >! closest-outside-${t}.stats.tsv
  cat closest-outside-${t}.tsv | cut -f8 | grep blost | cut -d':' -f1 | cut -d'/' -f3 | tools-countuniq -s | cut -f2 | sed "s/^/$t\t/" | tools-mergeuniq -merge | tools-vectors div -n 3 >> $log
  echo >> $log

#  echo "${t} and closest boundaries (inside):"
#  bedtools closest -a ICGC_${t}_uniqid_inv.bed -b boundary_strength_data_sorted.bed -t all -io >! closest-inside-${t}.tsv
#  cat closest-inside-${t}.tsv | cut -f8 | grep blost | cut -d':' -f1 | cut -d'/' -f3 | sort | uniq -c

  # allow some flexibility, don't require 100% overlap
  echo "co-${t} SEs and closest boundaries:" >> $log          
  bedtools intersect -a SE_closest_boundaries.bed -b ICGC_${t}_uniqid.bed -f $f -wa -wb >! co-${t}-SE.tsv
  cat co-${t}-SE.tsv | cut -f4 | tr '/' '\t' | cut -d'|' -f1 | tools-countuniq -s | sort | tools-table -c -n 0 | scripts-skipn 1 | sed 's/ $//' | tr -s ' ' '\t' >! co-${t}-SE.boxplot.tsv
  cat co-${t}-SE.boxplot.tsv | tools-matrix cstat | grep -E '^AVG|STD' | tools-matrix T >! co-${t}-SE.stats.tsv
  cat co-${t}-SE.tsv | cut -f4 | cut -d'|' -f1 | cut -d'/' -f2 | tools-countuniq -s | cut -f2 | sed "s/^/co-$t\t/" | tools-mergeuniq -merge | tools-vectors div -n 3 >> $log
  echo >> $log

  echo "co-${t} enhancers and closest boundaries:" >> $log
  bedtools intersect -a enhancer_closest_boundaries.bed -b ICGC_${t}_uniqid.bed -f $f -wa -wb >! co-${t}-enhancer.tsv
  cat co-${t}-enhancer.tsv | cut -f4 | tr '/' '\t' | cut -d'|' -f1 | tools-countuniq -s | sort | tools-table -c -n 0 | scripts-skipn 1 | sed 's/ $//' | tr -s ' ' '\t' >! co-${t}-enhancer.boxplot.tsv
  cat co-${t}-enhancer.boxplot.tsv | tools-matrix cstat | grep -E '^AVG|STD' | tools-matrix T >! co-${t}-enhancer.stats.tsv
  cat co-${t}-enhancer.tsv | cut -f4 | cut -d'|' -f1 | cut -d'/' -f2 | tools-countuniq -s | cut -f2 | sed "s/^/co-$t\t/" | tools-mergeuniq -merge | tools-vectors div -n 3 >> $log
  echo >> $log
   
  # genes
  cat co-${t}-SE.tsv | grep blost.05 | cut -f5-8 | gtools-overlaps overlap -i -label protein_coding.bed | sed 's/.*://' | sort | uniq -c | sort -rn >! co-${t}-SE-gene-freq.txt
  cat co-${t}-enhancer.tsv | grep blost.05 | cut -f5-8 | gtools-overlaps overlap -i -label protein_coding.bed | sed 's/.*://' | sort | uniq -c | sort -rn >! co-${t}-enhancer-gene-freq.txt
end
  
  # Get the other r/version
  module unload r
  module load r/3.3.2

  # Now pass to R what is needed to make calculations and create the graph
  Rscript code/calculate_fractions.r 'ICGC_duplications_uniqid.bed' 'ICGC_deletions_uniqid.bed' 

cd $mpath

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


