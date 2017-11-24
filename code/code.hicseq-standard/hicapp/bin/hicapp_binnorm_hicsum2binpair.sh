#!/bin/bash
##
# Generate bin based interaction contact pairs from hicsum file
##

## set par
input=$1
res=$2

## process
awk -v res=$res '{rd1=$1"\t"int($2/res)*res; rd2=$5"\t"int($6/res)*res; print rd1"\t"rd2}' $input | sort| uniq -c |awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' |sed 's/chrX/chr97/; s/chrY/chr98/; s/chrM/chr99/' |sort -t $'\t' -k1.4,1.5n -k2,2n -k3.4,3.5n -k4,4n |sed 's/chr97/chrX/; s/chr98/chrY/; s/chr99/chrM/' |awk '{print $1"_"$2"\t"$3"_"$4"\t"$5}'

