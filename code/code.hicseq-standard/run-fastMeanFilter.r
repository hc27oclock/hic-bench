#!/usr/bin/env Rscript
usage = "\
Rscript run-fastMeanFilter.r [OPTIONS] MATRIX-1
"
sessionInfo()
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-p","--prep"),default="none", help="matrix preprocessing (none or log2)"),
  make_option(c("-s","--smoothing"),default="0,1,2,3,4,5", help="smoothing parameter decides the neighborhood size of smoothing, assign multiple using comma to separate them"),
  make_option(c("-o","--output"),default="./RData.RData", help="name of the output .RData file")
)

arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) != 1) {
  write("Error: number of input matrix must be 1! Use --help to see help information", stderr()); quit(save='no')
}

suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(hicrep))
sourceCpp(file = "code/fastMeanFilter.cpp")

lambdas=as.numeric(unlist(strsplit(opt$s,",")))
x = as.matrix(read.table(inputs[1], check.names=FALSE))
y = x
prep = opt$p
if (prep=='log2') { write("Performing log2 transformation on input matrix...", stderr()); y = log2(y+1) }
y[is.na(y)] = 0

opt$'max-lambda' = max(lambdas)
opt$'min-lambda' = min(lambdas)
opt$preprocess = 'none'

solObj=array(dim=c(length(lambdas),dim(y)))
gammas = 0
ignored_cols = c()
ignored_rows = c()

for(h in seq(1:dim(solObj)[1])){
	solObj[h,,] = fastMeanFilter(y, lambdas[h])
  if (prep=='log2') { write("Perform reverse-log2 transformation on output matrix...", stderr()); solObj[h,,] = 2^solObj[h,,] - 1 }
}

save.image(opt$o)
