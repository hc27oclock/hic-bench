##
# convert bam to hicsum format
##

#### par set
bam=$1
bname=${bam%.bam}
re=$2
out=$3

#### bam2bedpe
bamToBed -bedpe -i $bam > $bname.bedpe


#### get the fragment mid position
awk '{if($9=="+") print $1"\t"$2"\t"$2+1"\t"$9; else print $1"\t"$3-1"\t"$3"\t"$9;}' $bname.bedpe > $bname.read1
bedtools intersect -a $bname.read1 -b $re -wao |awk '{print int(($6+$7)/2)}' > $bname.frag1
awk '{if($10=="+") print $4"\t"$5"\t"$5+1"\t"$10; else print $4"\t"$6-1"\t"$6"\t"$10;}' $bname.bedpe > $bname.read2
bedtools intersect -a $bname.read2 -b $re -wao |awk '{print int(($6+$7)/2)}' > $bname.frag2


#### get final contact output
paste $bname.read1 $bname.frag1 $bname.read2 $bname.frag2 |cut -f1,2,4,5,6,7,9,10 > $out


#### clean files
rm -rf $bname.read1 $bname.frag1 $bname.read2 $bname.frag2 $bname.bedpe



