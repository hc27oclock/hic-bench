##
# calculate copy number for each bin region
##

## set par
args <- commandArgs(TRUE)
segfile<-args[1] # "ts543_ck_cn.seg"
binfile<-args[2] # "hg19.genome.1m"
outfile<-args[3] # "ts543_ck_cn.1m"

## load library
library(CNTools)

## function to get the copy number on genomic regions
getcn<-function(segfile, binfile, outfile){
  cn <- read.delim(segfile, as.is = TRUE)
  seg <- CNSeg(cn)
  aa<-read.delim(binfile, header=F)
  colnames(aa)<-c("chrom", "start", "end", "geneid", "genename")
  aa$chrom<-gsub("chr", "", aa$chrom)
  rsseg <- getRS(seg, by = "gene", imput = FALSE, XY = FALSE, what = "mean", geneMap = aa)  
  out<-rsseg@rs[,c(4,6)]
  write.table(out, file=outfile, quote=F, sep="\t", row.names = F, col.names=F)
}

## run copy number call on fragment and bin window
getcn(segfile, binfile, outfile)