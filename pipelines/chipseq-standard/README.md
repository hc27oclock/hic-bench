# ChIP-Seq Standard Pipeline

User manual for the pipeline can be found here:

https://github.com/NYU-BFX/hic-bench_documentation/blob/master/HiC-Bench_manual.pdf.

ChIP-Seq pipeline usage follows the same structure and format as the HiC pipleline instructions listed there.

## 0. Create a new analysis

If a new analysis project has not already been created, do so with the following command:

```
~/hic-bench/code/code.main/pipeline-new-analysis chipseq-standard /path/to/<project_directory>
```

## 1. Set input files

### Option 1: manual

Within the corresponding `<project_directory>/inputs/fastq` or `<project_directory>/inputs/bam` directory, subdirectories should be created with the name of each sample to be included in the analysis. The following naming scheme is preferable:

\<Cell_line\>-\<ChIP\>-\<treatment\>-\<SampleID\>

Each subdirectory should contain all fastq or bam files to be used for that sample through the analysis pipeline. Symlinks can be used if the files are not contained in the same location as the project analysis directory, and are preferable to save storage space. 

*QQ: Include sample directory file structures*
*QQ: Include sample code for creating symlinks and dirs?*
*QQ: treatment of special characters, etc., in file/dir names*

### Option 2: automatic (may fail with certain file names)

From the project directory directory, run:
```
./code/setup-sample-files.sh <fastq_source_dir>
```
This will scan the given directory and create symlinks to any found FASTQs in `./inputs/fastq`. It will also clean up file names if they are in the standard Illumina bcl2fastq structure so they can be recognized by `create-sample-sheet.tcsh` (next step). It can be run multiple times to scan multiple directories, but outputs may get overwritten if the sample names are the same. Found FASTQs are printed to keep track of what is happening.

## 2. Create project sample sheet

A sample sheet must be created for the analysis project. Run the follow command to do so:

```
<project_directory>/inputs$ ./code/create-sample-sheet.tcsh <genome> <fragment-size>
```

Where `<genome>` is `hg19`, `hg38`, etc.. The `<fragment-size>` entry is optional and should be a numeric argument such as `300`, representing the library size of the sequencing sample. After creation of the sample sheet, output in `inputs/sample-sheet.tsv`, a manual review process is required to match the correct control or input samples with experimental samples, verify proper grouping names, files, and other entries. If not entered prior, `<fragment-size>` should be filled in for each sample. This process can be completed within Microsoft Excel, but saving the file in Excel should be avoided due to the introduction of formatting errors by Excel. It is advisable to instead copy the finalized sheet from Excel and paste directly into a terminal text editor such as `vi` or `nano` for saving.

## 3. Pipeline execution

Run the pipeline with:

```
<project_directory>$ ./code.main/pipeline-execute PROJECT-NAME E-MAIL
```

## Notes

Errors encountered during pipeline execution can be viewed with:

```
<project_directory>$ code.main/pipeline-errors
```

Analysis results can be removed with:

```
<project_directory>$ code/clean-all
```


*QQ: binaries*
*QQ: Creating a new pipeline step*
