# README for HiCapp v1.0.0 #
--------

HiCapp is a Hi-C analysis pipeline which can correct for the copy number bias in tumor Hi-C data using caICB correction algorithm. HiCapp receives Hi-C pair-end sequence reads to generate a corrected Hi-C map in interested binning resolutions. HiCapp currectly have three modules: ***hicapp_preprocess***, ***hicapp_binnorm*** and ***hicapp_bias***.

* ***hicapp_preprocess*** starts from mate-pair fastq files to generate cleaned contact pairs (named hicsum).
* ***hicapp_binnorm*** takes hicsum file to generate a corrected Hi-C map in certain binning resolution.
* ***hicapp_caICB*** takes binned Hi-C contact pair file to generate a corrected Hi-C map.
* ***hicapp_bias*** computes the four explicit biases in Hi-C data used for evaluation of normalization.

### How to download HiCapp? ###
--------
Download the pipeline using git clone or using the "download" link.

     git clone https://bitbucket.org/mthjwu/hicapp

### Software dependencies ###
--------
The pipeline has been tested in Linux system. 
The following softwares need user to install by themselves.

* R (https://www.r-project.org/)
* perl (https://www.perl.org/)
* bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* bedtools (http://bedtools.readthedocs.org/en/latest/)
* (optional) R-packages: CNTools 

The following softwares are already included in the pipeline, so don't need user to install by themselves.

* Hi-Corrector (http://zhoulab.usc.edu/Hi-Corrector/)
* HiCup (http://www.bioinformatics.babraham.ac.uk/projects/hicup/)

### Files required ###
--------
* Reference genome fasta file (e.g.: hg19.fa)
* bowtie2 index files of reference genome (e.g.: hg19.1.bt2, ...)
* Chromosome size file (e.g.: for hg19 https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes)
* (optional) Copy number segmentation file of the Hi-C analyzed sample (e.g.: K562.seg)
* (optional) CRG Mappability file (e.g.: for hg19 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig)

### How to run HiCapp pipeline from raw fastq data? ###
--------
Let's take K562 Hi-C data from study GSE63525 as an example. We want to use human genome version 19 (hg19) as the reference genome.

* Install R, perl, bowtie2, bedtools

* Prepare hg19.fa, bowtie2 index (hg19.1.bt2...) and chromosome size file (hg19.genome: https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes)

* Run mapping and preprocessing to get the hicsum file:

        hicapp_preprocess -1 K562_R1.fq.gz -2 K562_R2.fq.gz -g hg19.fa -e ^GATC,MboI

* Run binning and normalization on 1Mb resolution to get the bias vector in <K562.1m.bychrom.caicb> file:

        hicapp_binnorm -s K562.sorted.hicsum -r 1000000 -t 1m -c hg19.genome

* Then one can easily get the caICB normalized counts from <pair> file and <caicb> file (see "File format description" section):

        awk 'NR==FNR {a[$1]=$2} NR>FNR {B=a[$1]*a[$2]; if (B>0) print $0"\t"$3/B}' K562.1m.bychrom.caicb K562.1m.pair > K562.1m.pair.norm

* Calculate explicit biases on 1Mb resolution: (R-packages: CNTools is needed for copy number bias calculation; CRG Mappability file is needed for mappability calculation)

        hicapp_bias -g hg19.fa -e ^GATC,MboI -c hg19.genome -m mappability100mer.bw -r 1000000 -t 1m -z K562.seg


### How to perform caICB correction from binned raw Hi-C contact counts? ###
--------
Let's take K562 Hi-C data from study GSE63525 as an example. We want to use human genome version 19 (hg19) as the reference genome.

* Install R, perl

* Prepare chromosome size file (hg19.genome: https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes)

* Prepare binned Hi-C contact map in certain resolution (Let's say K562 in 1Mb binning resolution: K562.1m.pair, check the format of <pair> file in "File format description" section) 

* Run caICB correction to get the bias vector in <K562.1m.bychrom.caicb> file:

        hicapp_caICB -s K562.1m.pair -r 1000000 -c hg19.genome

* One can easily get the caICB normalized counts from <pair> file and <caicb> file (see "File format description" section):

        awk 'NR==FNR {a[$1]=$2} NR>FNR {B=a[$1]*a[$2]; if (B>0) print $0"\t"$3/B}' K562.1m.bychrom.caicb K562.1m.pair > K562.1m.pair.norm


### File format description ###
--------
* <hicsum> file has 8 columns (the file has no header, position is one-indexed):

        <read1 chromosome> <read1 position 5' end> <read1 strand> <read1 insilico enzyme cutting site> <read2 chromosome> <read2 position 5' end> <read2 strand> <read2 insilico enzyme cutting site>
        chr1	10242	+	8003	chr8	129590490	+	129587691
        chr1	13055	+	8003	chr8	130235011	-	130236586

* <pair> file is the binned Hi-C contact map in certain resolution, it has 3 columns (the file has no header, position is zero-indexed):

        <chr_start of anchor 1> <chr_start of anchor 2> <counts>
        chr1_0  chr1_0          12
        chr1_0  chr1_1000000    66
        chr1_0  chr1_2000000    19
        chr1_0  chr1_3000000    11
  
* <caicb> file contains the bias vector estimated by caicb algorithms, it has 2 columns (the file has no header, position is zero-indexed):

        <chr_start of bin> <bias vector>
        chr1_0          0.112787822292498
        chr1_1000000    0.690712993542956
        chr1_2000000    0.946608928048394
        chr1_3000000    0.995235490936498


### Who do I talk to? ###
--------
* Hua-Jun Wu (hjwu@jimmy.harvard.edu)