#!/usr/bin/env Rscript
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('tidyr'))
suppressPackageStartupMessages(library('stringr'))
suppressPackageStartupMessages(library('lubridate'))
suppressPackageStartupMessages(library('readr'))
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('cowplot'))
suppressPackageStartupMessages(library('pheatmap'))
suppressPackageStartupMessages(library('tools'))
suppressPackageStartupMessages(library('ggpubr'))

args <- commandArgs(trailingOnly = TRUE)

for (file in strsplit(args[1],",") %>% unlist()){
ctcf=read_delim(file,delim = '\t',col_names = F)
colnames(ctcf) <- c("sample","Boundary Strength Category","Normalized Aggregate CTCF Level")
#ctcf = ctcf %>% mutate(`Normalized Aggregate CTCF Level`=as.numeric(`Normalized Aggregate CTCF Level`)) %>% mutate(`Boundary Strength Category`=gsub("b","",`Boundary Strength Category`))
ctcf = ctcf %>% mutate(`Normalized Aggregate CTCF Level`=as.numeric(`Normalized Aggregate CTCF Level`)) %>% mutate(`Boundary Strength Category`=gsub("b1","I",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b2","II",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b3","III",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b4","IV",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b5","V",`Boundary Strength Category`))
figure4e=ctcf %>% 
  ggplot(aes(x=`Boundary Strength Category`, y=`Normalized Aggregate CTCF Level`, fill=`Boundary Strength Category`)) + 
  scale_fill_manual(values=c("#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5")) +
  geom_boxplot(width = 0.5) + 
  geom_jitter(size=0.5) + theme_bw() +
  theme(axis.text=element_text(size=8),axis.title=element_text(size=8),legend.position="none",plot.title = element_text(size=12, hjust = 0),strip.text = element_text(size = 8)) +
  labs(y="Normalized Aggregate\nCTCF Level",x="Boundary Strength Category\n(I: Weakest, V: Strongest)") +
  stat_compare_means(comparisons =list( c("V", "IV"), c("V", "III"), c("V", "II"), c("V", "I")),method = "wilcox.test", paired = T) 
ggsave(paste0(file_path_sans_ext(file),".pdf"), plot = last_plot(), width = 2.4439*2, height = 3.3067)
}

se=read_delim(args[2],delim = '\t',col_names = F)
colnames(se) <- c("sample","Boundary Strength Category","Fraction of Super-enhancers")
#se = se %>% mutate(`Fraction of Super-enhancers`=as.numeric(`Fraction of Super-enhancers`)) %>% mutate(`Boundary Strength Category`=gsub("0","",`Boundary Strength Category`))
se = se %>% mutate(`Fraction of Super-enhancers`=as.numeric(`Fraction of Super-enhancers`)) %>% mutate(`Boundary Strength Category`=gsub("b1","I",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b2","II",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b3","III",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b4","IV",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b5","V",`Boundary Strength Category`))
figure4f=se %>% 
  ggplot(aes(x=`Boundary Strength Category`,y=`Fraction of Super-enhancers`,fill=`Boundary Strength Category`)) + 
  scale_fill_manual(values=c("#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5")) +
  geom_boxplot(width = 0.5) + 
  geom_jitter(size=0.5) + theme_bw() +
  theme(axis.text=element_text(size=8),axis.title=element_text(size=8),legend.position="none",plot.title = element_text(size=12, hjust = 0),strip.text = element_text(size = 8)) +
  labs(y="Fraction\nof Super-enhancers",x="Boundary Strength Category\n(I: Weakest, V: Strongest)") +
  stat_compare_means(comparisons =list( c("V", "IV"), c("V", "III"), c("V", "II"), c("V", "I")),method = "wilcox.test", paired = T)
ggsave(paste0(file_path_sans_ext(args[2]),".pdf"), plot = last_plot(), width = 2.4439*2, height = 3.3067)

df=read_delim(args[3],delim = '\t',col_names = T)
df = df %>% as.data.frame()
rownames(df)=df[,1]
df=df[,-1]
kmean_result=kmeans(df,10)
pheatmap(df[order(kmean_result$cluster),],cluster_rows = F,show_rownames = F,fontsize=5.81, width = 4.2935*2,height = 3.4705*2, treeheight_col=15,filename=paste0(file_path_sans_ext(args[3]),".pdf"))

