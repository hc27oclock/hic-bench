#!/usr/bin/env Rscript
usage = "\
Rscript [this script] Input OUTDIR
"
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('tidyr'))
suppressPackageStartupMessages(library('readr'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('ggpubr'))
expand.grid.unique <- function(x, y, include.equals=FALSE){
  x <- unique(x);  y <- unique(y)
  g <- function(i){
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}


args <- commandArgs(trailingOnly = TRUE)
dir.create(args[2], showWarnings = FALSE)
all=read_delim(args[1],delim = '\t',col_names = F) %>% mutate(X4=gsub("0","b",X4)) %>% mutate(X4=gsub("b1","I",X4)) %>% mutate(X4=gsub("b2","II",X4)) %>% mutate(X4=gsub("b3","III",X4)) %>% mutate(X4=gsub("b4","IV",X4)) %>% mutate(X4=gsub("b5","V",X4)) %>% mutate(X5=as.numeric(X5)) %>% as.data.frame()
colnames(all) <- c("Method","kappa","sample", "Boundary Strength Category","Normalized Aggregate CTCF Level")
for (myMethod in unique(all$Method)){
  data= all %>% filter(Method==myMethod)
  ggplot(data,aes(x = kappa, y = `Normalized Aggregate CTCF Level`, fill=`Boundary Strength Category`, color=`Boundary Strength Category`)) +
    stat_smooth() +
    theme(legend.position="top", legend.direction="horizontal",plot.title = element_text(hjust = 0.5)) +
    ggtitle(myMethod) +
    scale_x_continuous(breaks=sort(unique(data$kappa)))
  ggsave(paste0(args[2],"/",myMethod,".pdf"), plot = last_plot(), width = 6, height = 8)
}

for (mykappa in unique(all$kappa)){
  data= all %>% filter(kappa==mykappa)
  u=data %>% select(Method) %>% unique
  grp=u %>% separate(Method,c("V1","V2"),"-") %>% select(V2) %>% unique()
  for (mymethod in grp[,1]){
    comp=expand.grid.unique(u[,1],u[,1]) %>% as.data.frame() %>% filter(grepl(mymethod,V1)) %>% filter(grepl(mymethod,V2))
    if(dim(comp)[1]==1){
    myplot=data %>% filter(Method==comp[,1]|Method==comp[,2]) %>% separate(Method,c("Method"),"\\.") %>%
      ggplot(aes(x=`Boundary Strength Category`, y=`Normalized Aggregate CTCF Level`)) +
      geom_boxplot(aes(fill=Method),width = 0.5) +
      theme_bw() +
      theme(axis.text=element_text(size=8),axis.title=element_text(size=8),strip.text = element_text(size = 8),legend.position="top", legend.direction="vertical",legend.text=element_text(size=8),legend.title=element_blank()) +
      labs(y=ifelse(args[2]!="compare-SE","Normalized Aggregated CTCF Level","Fraction of Super-Enhancer"),x="Boundary Strength Category (I: Weakest, V: Strongest)") +
      stat_compare_means(aes(group = Method), label = "p.format",method = "wilcox.test", paired = T,size=2) +
      guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
      scale_fill_manual(values=c("#E6550D",ifelse(args[2]!="compare-SE","#31A354","#756BB1")))
    ggsave(paste0(args[2],"/",comp[,1],"_vs_",comp[,2],".kappa",mykappa,".pdf"), plot = myplot, width = 5.5, height = 2.5)
save.image(file = paste0(args[2],"/",comp[,1],"_vs_",comp[,2],".kappa",mykappa,".RData"))
    }
  }
}

