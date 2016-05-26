#!/bin/Rscript

# Check for required packages and install if anything is missing
list.of.packages <- c("ggplot2", "plyr", "data.table", "reshape2", "gridExtra", "scales")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

# Required packages
library(plyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(scales)

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)

# Usage
if (length(args)!=2){print("USAGE: plot-domains-stats2.R OUTPUT-DIR INPUT-DIR"); quit(save="no")}

# Get the arguments
output <- args[1]
input  <- args[2]

dirs  <- system(sprintf("ls -d %s/*", input), intern=TRUE)
samples <- c()

#This is a list where all the information from
#all samples, lambdas, chroms will be saved
l1 <- list()

# Get into each one of the directories
for (i in 1:length(dirs)) {
	sample <- strsplit(basename(dirs[i]),"[-]")[[1]][1]
	files <- system(sprintf("ls -1 %s/domains*", dirs[i]), intern=TRUE)
	l2 <- list()
	for (j in 1:length(files)) {
		f <- read.table(sprintf("%s", files[j]), header=FALSE, stringsAsFactors=FALSE)
		domain_sizes <- f[,3] - f[,2]
		domain_number <- rep(length(domain_sizes), length(domain_sizes))
		f1 <- data.frame(cbind(f[,1],domain_sizes))
		f1[,2] <- as.numeric(as.character(domain_sizes))
		f1[,3] <- as.numeric(as.character(domain_number))
		kappa_temp <- strsplit(basename(files[j]), "[.]")[[1]][2]
		kappa <- strsplit(kappa_temp, "[=]")[[1]][2]
		# Get the sample column for the dataframe
		sample_names <- rep(sample, nrow(f))
		# Get the kappa column for the dataframe
		kappa_names <- rep(kappa, nrow(f))
		df <- data.frame(cbind(sample_names,kappa_names,f1))
		df[,1] <- as.factor(df[,1])
		df[,2] <- as.factor(df[,2])
		df[,3] <- as.factor(df[,3])
		df[,4] <- as.numeric(as.character(df[,4]))
		colnames(df) <- c("Group","Kappa","Chrom","Size","Number")
		l2[[j]] <- df
	}
	l1[[i]] <- rbindlist(l2)
}

# Convert the list to data frame 
domain_size.df <- rbindlist(l1)

# Summarize stats per chromosome
final.df <- ddply(domain_size.df, c("Group", "Kappa", "Chrom"), summarise, mean=mean(Size), min=min(Size), max=max(Size), rom=(max(Size)-min(Size))/mean(Size), N=length(Size), sd = sd(Size), se = sd/sqrt(N), seom=se/mean)

# Summarize stats per group
final.df1 <- ddply(domain_size.df, c("Group", "Kappa"), summarise, mean=mean(Size), min=min(Size), max=max(Size), rom=(max(Size)-min(Size))/mean(Size), N = length(Size), sd = sd(Size), se = sd/sqrt(N), seom=se/mean)

# Summarize domain numbers per group
domain_size.df2 <- subset(domain_size.df, domain_size.df$Group!="CD34")
final.df2 <- ddply(domain_size.df2, c("Kappa"), summarise, mean_number=mean(Number), min_number=min(Number), max_number=max(Number), rom_number=(max(Number)-min(Number))/mean(Number), sd_number = sd(Number), se_number = sd_number/sqrt(length(Number)), seom_number=se_number/mean_number) 

# Domain stats per chromosome
write.table(final.df, sprintf("%s/domain_stats_per_chr.tsv",output), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

# Domain stats per group
write.table(final.df1, sprintf("%s/domain_stats_per_group.tsv",output), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

# Domain number stats per group
write.table(final.df2, sprintf("%s/domain_number_stats_per_lambda.tsv",output), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
