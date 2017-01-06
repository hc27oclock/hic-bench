# ChIP-Seq Standard Pipeline

User manual for the pipeline can be found here:

https://github.com/NYU-BFX/hic-bench_documentation/blob/master/HiC-Bench_manual.pdf.

ChIP-Seq pipeline usage follows the same structure and format as the HiC pipleline instructions listed there.


## Overview

```bash
# create project directory
~/hic-bench/code/code.main/pipeline-new-analysis chipseq-standard /path/to/<project_directory>

# set input files
cd /path/to/<project_directory>
./code/setup-sample-files.sh <fastq_source_dir>

# create sample sheet template from samples
cd inputs
./code/create-sample-sheet.tcsh <hg19|mm10> <fragment-size>

# !! FILL IN THE CONTROL AND SAMPLE GROUPINGS IN THE SAMPLE SHEET MANUALLY !!

# run the pipeline
cd /path/to/<project_directory>
./code.main/pipeline-execute <project_ID|project_directory_name> <your.email.goes.here@address.edu>

```


# Full Walkthrough
## 0. Create a new analysis

If a new analysis project has not already been created, do so with the following command:

```bash
~/hic-bench/code/code.main/pipeline-new-analysis chipseq-standard /path/to/<project_directory>
```

## 1. Set input files

There are two methods to set up your sample files. Pick one of the following:

### Automatic Method (preferred)

From the project directory directory, run:
```bash
./code/setup-sample-files.sh <fastq_source_dir>
```
This will scan the given directory and create symlinks to any found FASTQs in `./inputs/fastq`. It will also clean up file names if they are in the standard Illumina bcl2fastq structure so they can be recognized by `create-sample-sheet.tcsh` (next step). It can be run multiple times to scan multiple directories, but outputs may get overwritten if the sample names are the same. Found FASTQs are printed to keep track of what is happening.

NOTE: Certain sample names could cause issues with this, so be sure to check the output. 

### Manual Method 

Within the corresponding `<project_directory>/inputs/fastq` or `<project_directory>/inputs/bam` directory, subdirectories should be created with the name of each sample to be included in the analysis. The following naming scheme is preferable:

\<Cell_line\>-\<ChIP\>-\<treatment\>-\<SampleID\>

Each subdirectory should contain all fastq or bam files to be used for that sample through the analysis pipeline. Symlinks can be used if the files are not contained in the same location as the project analysis directory, and are preferable to save storage space. 



## 2. Create project sample sheet

A sample sheet must be created for the analysis project. Run the follow command to do so:

```bash
<project_directory>/inputs$ ./code/create-sample-sheet.tcsh <genome> <fragment-size>
```

Where `<genome>` is `hg19`, `hg38`, etc.. The `<fragment-size>` entry is optional and should be a numeric argument such as `300`, representing the library size of the sequencing sample. After creation of the sample sheet, output in `inputs/sample-sheet.tsv`, a manual review process is required to match the correct control or input samples with experimental samples, verify proper grouping names, files, and other entries. If not entered prior, `<fragment-size>` should be filled in for each sample. This process can be completed within Microsoft Excel, but saving the file in Excel should be avoided due to the introduction of formatting errors by Excel. It is advisable to instead copy the finalized sheet from Excel and paste directly into a terminal text editor such as `vi` or `nano` for saving.

## 3. Pipeline execution

Run the pipeline with:

```bash
<project_directory>$ ./code.main/pipeline-execute <project_ID> <your_email@address.edu>
```

## 4. Compile Report

A report template has been supplied in the `report` directory. It will automatically scan the pipeline output, and generate a PDF which includes sample sheets and figures found in the report. 

First, the project info text file needs to be updated with the correct parameters for the project. 


```bash
$ cat project_info.txt
PROJECT-DIR: /project/dir

PROJECT-ID: Project_ID_ChIP-Seq

PROJECT-ID-SHORT: ChIP-Seq_ID

REPORT-AUTHOR: Stephen M. Kelly

REPORT-AUTHOR-EMAIL: stephen.kelly@nyumc.org

PI-NAME: Dr. X

$ nano project_info.txt
# edit the file

$ cat project_info.txt
PROJECT-DIR: /ifs/home/kellys04/projects/SmitLab_Mike_ChIPSeq_2017-01-05

PROJECT-ID: SmitLab_Mike_ChIPSeq_2017-01-05

PROJECT-ID-SHORT: Mike_ChIPSeq

REPORT-AUTHOR: Stephen M. Kelly

REPORT-AUTHOR-EMAIL: stephen.kelly@nyumc.org

PI-NAME: Dr. Tsirigos
```

```bash
<project_directory>/report$ ./compile_report_wrapper.sh chipseq_report.Rnw
```

## Notes

Errors encountered during pipeline execution can be viewed with:

```bash
<project_directory>$ code.main/pipeline-errors
```

Analysis results can be removed with:

```bash
<project_directory>$ code/clean-all
```

