#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript




################################
##### MAIN   ###################
################################


args <- commandArgs(trailingOnly=T)
if (length(args)!=3) {
  cat("USAGE: create-random-vector.r N-SELECT N-TOTAL SEED\n")
  quit(save="no")
}

n_select = as.integer(args[1])
n_total = as.integer(args[2])
set.seed(as.integer(args[3]))

n_select = min(n_select,n_total)
y = sample(c(rep(1,n_select),rep(0,n_total-n_select)))

write.table(y,quote=F,col.names=F,row.names=F)



