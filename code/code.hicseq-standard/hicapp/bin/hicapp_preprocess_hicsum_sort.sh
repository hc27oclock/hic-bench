#!/bin/bash
##
# sort hicsum file based on genomic coordinates
##

input=$1

sed 's/chrX/chr97/; s/chrY/chr98/; s/chrM/chr99/' $input |sort -t $'\t' -k1.4,1.5n -k2,2n -k5.4,5.5n -k6,6n |sed 's/chr97/chrX/; s/chr98/chrY/; s/chr99/chrM/'

