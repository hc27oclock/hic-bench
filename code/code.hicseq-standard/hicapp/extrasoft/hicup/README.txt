HiCUP (Hi-C User Pipeline) - Scripts for analysing Hi-C sequence data
---------------------------------------------------------------------

HiCUP is a bioinformatics pipeline for processing Hi-C data. The pipeline 
receives FASTQ data which is then mapped against a reference genome and 
filtered to remove frequently encountered experimental artefacts. The 
pipeline produces paired read files in SAM/BAM format, each read pair 
corresponding to a putative Hi-C di-tag.

HiCUP comprises several Perl scripts:
i)   hicup_digester
ii)  hicup (does not analyse data but regulates data flow through the pipeline)  
iii) hicup_truncater  
iv)  hicup_mapper  
v)  hicup_filter
vi) hicup_deduplicator

along with an HTML file hicup_report.html.

HiCUP requires a working version of Perl and Bowtie installed on your machine, 
which should be running a Unix-based operating system and have gzip installed.

For further details go to:
http://bowtiebio.sourceforge.net/index.shtml

HiCUP is written in Perl and executed from the command line. To install 
HiCUP, copy the hicup_v0.X.Y.tar.gz file into a HiCUP installation folder 
and extract all files by typing:
tar xzf hicup_v0.X.Y.tar.gz

If you have any comments about HiCUP we would like to hear them. You
can either enter them in our bug tracking system at:

http://www.bioinformatics.babraham.ac.uk/bugzilla

..or send them directly to steven.wingett@babraham.ac.uk
