#!/usr/bin/env Rscript
usage = "\
Rscript optimal_lambda.r MATRIX
"

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
outdir=args[1]
inputs=args[2]
min_improvement=as.numeric(args[3])

compare=read.table(inputs,header=T)
for (k in compare$LAMBDA %>% sort %>% unique %>% -1){
  if(k!=0){
    p=wilcox.test((1.0+min_improvement)*compare[compare$LAMBDA==k,]$VALUE,compare[compare$LAMBDA==k+1,]$VALUE,alternative = 'less')$p.value
    cat(paste(unique(compare$SAMPLE.1),unique(compare$SAMPLE.2),k,p,sep="\t"))
    cat("\n")
  }
}


outplot = ggplot(compare,aes(x=as.factor(LAMBDA),y=VALUE,fill=as.factor(LAMBDA))) + geom_boxplot() + xlab("Lambda") + ylab("Reproducibility") + labs(fill='Lambda') + scale_fill_manual(values=c("#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5")) 

ggsave(paste(outdir,"boxplots.pdf",sep='/'),plot=outplot)


