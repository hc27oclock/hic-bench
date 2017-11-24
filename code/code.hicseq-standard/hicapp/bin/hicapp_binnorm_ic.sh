####
# normalize matrix by IC (Hi-Corrector1.1 software used)
####

## set path of the pipeline
bindir=`dirname $0`
bdir=$bindir/..
icedir=$bdir/extrasoft/Hi-Corrector/bin

## set par
matrix=$1 #contact matrix
cmd=$2 #ic or ic_mes or ic_mep
ncpu=$3 #number cpu
iter=$4 #number of iteration
nrow=$[$(wc -l $matrix | awk '{print $1}')-1]

## start to normalize matrix
if [ ! -e $matrix.icb ]
then
	## remove diagonal of the matrix
	#awk -F $'\t' 'BEGIN {OFS = FS} NR==1 {print} NR>1 {$NR=0; print}' $matrix > $matrix.rmdiag

	## run Hi-Corrector
	if [ $cmd == "ic" ]
	then
		## run ic
		$icedir/ic $matrix $nrow $iter 1 1 ./$matrix.icb > ./$matrix.icb.log
	else
		## run in parallel
		/opt/openmpi/bin/mpirun -np $ncpu $icedir/ic_mep --inputFile=$matrix --numTask=$ncpu --memSizePerTask=1000 --numRows=$nrow --maxIteration=$iter --hasHeaderRow=1 --hasHeaderCol=1 --jobID=$matrix --outputFile=./$matrix.icb > ./$matrix.icb.log
	fi

	## add tag for the bias vector
	paste <(cut -f1 $matrix |awk 'NR>1') $matrix.icb > $matrix.icb.anno
	mv $matrix.icb.anno $matrix.icb

else	
	echo "#*_*# <$matrix.icb> already exists...Skip~~~"
fi













