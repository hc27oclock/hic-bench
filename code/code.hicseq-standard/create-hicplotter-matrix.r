#!/bin/Rscript

# Rscript create-hicplotter-matrix.r OUT-MATRIX INPUT-MATRIX [NORMALIZE=0]

# Get the arguments
args <- commandArgs(trailingOnly = TRUE)

out_mat <- args[1]
input_mat <- read.table(sprintf("%s", args[2]), header=TRUE, row.names=1, check.names=FALSE)

if ((length(args)>=3)&&(args[3]==1)) {   # normalize
  input_mat = as.matrix(input_mat)
  input_mat = input_mat/mean(input_mat,na.rm=TRUE)
}

# Write to output
write.table(input_mat, file=sprintf("%s", out_mat), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")



