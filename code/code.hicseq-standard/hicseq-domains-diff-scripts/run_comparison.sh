#!/bin/bash

script=$1
folder_processed_sample_1=$2
folder_processed_sample_2=$3
tads_sample_1=$4
tads_sample_2=$5

min_tad_size=$6
max_tad_size=$7
is_normalize=$8
centrotelo_file=$9
#params=$6
bin_size=${10}
out_folder=${11}
log=$out_folder/log.log

#is_normalize=TRUE
#bin_size=40000

echo "Estimating differences between sample_1 and sample_2 Hi-C data"
echo "Estimating differences between sample_1 and sample_2 Hi-C data" > $log

out_all=$out_folder/final_results.tsv

if [ ! -d $out_folder ]; then
	mkdir $out_folder
fi

rm -rf $out_all
touch $out_all

for matrix_file in $folder_processed_sample_1/matrix.chr*tsv; do
	matrix_file_base=`basename $matrix_file`
	chr=`echo "$matrix_file_base" | cut -f2 -d. | cut -f1 -d_`
	echo "Run comparison between sample_1 and sample_2 on $matrix_file_base"
	echo "Run comparison between sample_1 and sample_2 on $matrix_file_base" >> $log
	if [ ! -f $folder_processed_sample_2/matrix.${chr}.tsv ]; then
		echo "Cannot find file $folder_processed_sample_2/matrix.${chr}.tsv, checking next chromosome"
		echo "Cannot find file $folder_processed_sample_2/matrix.${chr}.tsv, checking next chromosome" >> $log
		continue
	fi
	matrix_file_base_wo_ext="${matrix_file_base%.*}"
#	chr=`echo "$matrix_file_base" | cut -f2 -d.`
#	echo "matrix_file_base_wo_ext = $matrix_file_base_wo_ext"
#	exit
	n=`wc -l $tads_sample_1 | cut -f1 -d" "`
	echo "R --no-save $matrix_file $folder_processed_sample_2/matrix.${chr}.tsv $tads_sample_1 $tads_sample_2 $chr $n $is_normalize $bin_size $centrotelo_file $max_tad_size \
		$out_folder/$matrix_file_base_wo_ext < $script" >> $log
	echo "R --no-save $matrix_file $folder_processed_sample_2/matrix.${chr}.tsv $tads_sample_1 $tads_sample_2 $chr $n $is_normalize $bin_size $centrotelo_file $max_tad_size \
		$out_folder/$matrix_file_base_wo_ext < $script"
	R --no-save $matrix_file $folder_processed_sample_2/matrix.${chr}.tsv $tads_sample_1 $tads_sample_2 $chr $n $is_normalize $bin_size $centrotelo_file $max_tad_size \
		$out_folder/$matrix_file_base_wo_ext < $script

	if [ -f $out_folder/${matrix_file_base_wo_ext}_results-table.tsv ]; then
		sed '1d' $out_folder/${matrix_file_base_wo_ext}_results-table.tsv >> $out_all
	fi
done

if [ -f $out_folder/matrix.chr1_results-table.tsv ]; then
	head -n 1 $out_folder/matrix.chr1_results-table.tsv | cat - $out_all > $out_all.new
	mv -f $out_all.new $out_all
fi

if [ ! -f $out_all ]; then
	(>&2 echo "Error: Results file $out_all is not created")	
fi






