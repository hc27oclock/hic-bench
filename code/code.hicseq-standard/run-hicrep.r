#!/usr/bin/env Rscript
usage = "\
Rscript run-hicrep.r [OPTIONS] MATRIX-1 MATRIX-2
"

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(hicrep))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(magrittr))

option_list <- list(
  make_option(c("-b","--bin"),default="40000", help="resolution of matrices, means the bin size"),
  make_option(c("-s","--smoothing"),default="0", help="smoothing parameter decides the neighborhood size of smoothing"),
  make_option(c("-r","--minreads"),default="inf", help="the size the total number one wants to adjust to, use inf for no adjustment, use auto for size equal to the smaller library"),
  make_option(c("-m","--maxdist"),default="./", help="contacts beyond this distance will not be considered")
)

arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) != 2) {
  write("Error: number of input matrix must be 2! Use --help to see help information", stderr()); quit(save='no')
}

bin_size=as.integer(opt$bin)
h=as.integer(opt$smoothing)
max_dist=as.integer(opt$maxdist)
# for hicrep version 1.0.0

#HiCR1=read.table(inputs[1], header = F, sep = "\t", skip = 1) %>% separate(V1, c("n1", "n2", "n3"), sep=":|-", remove=FALSE) %>% select(-V1) %>% set_colnames(paste0("V",seq(1:(dim(.)[2])))) %>% mutate(V2=as.integer(V2), V3=as.integer(V3))
#HiCR2=read.table(inputs[2], header = F, sep = "\t", skip = 1) %>% separate(V1, c("n1", "n2", "n3"), sep=":|-", remove=FALSE) %>% select(-V1) %>% set_colnames(paste0("V",seq(1:(dim(.)[2])))) %>% mutate(V2=as.integer(V2), V3=as.integer(V3))

# for hicrep version 1.0.1

HiCR1=as.matrix(read.table(inputs[1], header=T, check.names=F))
HiCR2=as.matrix(read.table(inputs[2], header=T, check.names=F))

# replace NAs with zeros
HiCR1[is.na(HiCR1)] = 0
HiCR2[is.na(HiCR2)] = 0

if(opt$minreads=="inf"){
  min_reads=as.integer("2000000")
  DS_HiCR1=HiCR1
  DS_HiCR2=HiCR2
}else{
  if(opt$minreads=="auto"){
    min_reads=min(sum(HiCR1[,-c(1:3)]),sum(HiCR2[,-c(1:3)]))
  }else{
    min_reads=as.integer(opt$minreads)
  }
  DS_HiCR1 <- depth.adj(HiCR1, min_reads)
  DS_HiCR2 <- depth.adj(HiCR1, min_reads)
}


SCC.out <- get.scc(DS_HiCR1, DS_HiCR2, bin_size, h, 0,  max_dist)
cat(as.numeric(SCC.out[[3]]))
