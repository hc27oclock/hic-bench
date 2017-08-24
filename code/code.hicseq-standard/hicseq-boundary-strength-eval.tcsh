#!/bin/tcsh

##
## USAGE: hicseq-boundary-strength-eval INP-DIRECTORY GENOME-DIR [CTCF|repeats] [filter={tss,noncoding,all}]
##

if ($#argv < 2) then
  grep '^##' $0
  exit
endif

set inpdir = $1
set genome_dir = $2
set regions = $3
set filter = $4

if ($regions == "") set regions = CTCF
set repeats = (Alu L1 MIR L2 Simple_repeat Low_complexity ERVL-MaLR hAT-Charlie ERV1 ERVL TcMar-Tigger CR1)

set p = `pwd`
set cell_types = (`cd $inpdir; ls -1 */blost.01.bed | cut -d'/' -f1 | cut -d'-' -f1 | sort -u`)
foreach c ($cell_types)
  set c1 = `echo $c | cut -d'-' -f1`
  if ($regions == "CTCF") then
    set bedfiles = CTCF:$p/CTCF_peaks/$c1-CTCF/peaks.bed
  else if ($regions == "repeats") then
    set bedfiles = 
    foreach repeat ($repeats)
      set bedfiles = ($bedfiles ${repeat}:$p/$repeat.bed)
    end
  else 
    echo "Error: provided region type not found!"
    exit
  endif

  foreach ref ($bedfiles)
    set peaks = `echo $ref | cut -d':' -f2-`
    set refname = `echo $ref | cut -d':' -f1`
    foreach bdir ($inpdir/$c-*)
      set sample = `echo $bdir | sed 's/.*\///'`
      cd $bdir
      echo -n "`basename $bdir`:$refname\t"
      foreach blost (blost.*)
        set lambda = `echo $blost | cut -d'.' -f2`
        set b_size = `cat $blost | gtools-regions n | cut -f2 | matrix csum -l`

        if ($filter == "noncoding") then
          set b_sum = `cat $peaks | cut -f-3,5 | gtools-overlaps subset -i -inv $p/$genome_dir/tss.bed | gtools-overlaps subset -i -inv $p/$genome_dir/gene.bed | gtools-overlaps subset -i $blost | cut -f4 | matrix csum -n 6`  # non-coding only
        else if ($filter == "tss") then
          set b_sum = `cat $peaks | cut -f-3,5 | gtools-overlaps subset -i $p/$genome_dir/tss.bed | gtools-overlaps subset -i $blost | cut -f4 | matrix csum -n 6`     # TSS only
        else 
          set b_sum = `cat $peaks | cut -f-3,5 | gtools-overlaps subset -i $blost | cut -f4 | matrix csum -n 6`      # all
        endif
        
        set b_score = `echo "1000000*$b_sum/$b_size" | bc -l | vectors format -n 3`               # sum of peak scores (or number of repeat elements) per boundary size (Mb)
        echo -n "$b_score "
      end
      echo ""
      cd $p
    end
  end
end




