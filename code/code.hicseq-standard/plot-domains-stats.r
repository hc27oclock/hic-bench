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
if (length(args)!=3){print("USAGE: plot-domains-stats.R OUTPUT-DIR INPUT-DIR KAPPA"); quit(save="no")}

# Get the arguments
output <- args[1]
input  <- args[2]
kappa  <- args[3]

file <- paste(input,'*',sprintf("domains.%s.bed", kappa),sep="/")

dirs  <- system(sprintf("ls -d %s/*", input), intern=TRUE)
files <- system(sprintf("ls -1 %s", file), intern=TRUE)

# Get the file number
file_no <- length(unique(files))

# Create a list to save all data frames
# the list will save all the dataframes with the sizes
l <- list()
# this vector will save all the numbers
n <- c()
# this will save all the names
names <- c()

for (i in 1:file_no) {
	file <- read.table(sprintf("%s", files[i]), header=FALSE)
	names[i] <- basename(dirs[i])
	sizes <- data.frame(file[,3]-file[,2])
	colnames(sizes) <- names[i]
	number <- nrow(sizes)
	l[[i]] <- sizes
	n[i] <- number
}

#Get the sizes of the domains
sizes_df <- rbindlist(l, use.names=TRUE, fill=TRUE)
#Convert to long format
sizes_long <- melt(sizes_df)
#Get the column names
colnames(sizes_long) <- c("sample","size")

#Get the domain numbers
domain_numbers <- data.frame(cbind(names,n), stringsAsFactors=FALSE)
colnames(domain_numbers) <- c("sample","number")
domain_numbers$sample <- as.factor(domain_numbers$sample)
domain_numbers$number <- as.numeric(as.character(domain_numbers$number))

#Output the domain stats
outfile <- paste(output,sprintf("domains-stats.%s.pdf", kappa),sep="/")

# Create the .pdf with two 
# plots on each page
pdf(sprintf("%s", outfile))
p1 <- ggplot(sizes_long, aes(x=factor(sample), y=size)) + geom_boxplot(aes(fill=factor(sample))) + scale_y_log10(limits = c(1e4,1e8)) + xlab("Sample") + ylab("TAD size (bp)") + theme(axis.text.y=element_text(size=7), axis.text.x = element_text(size=7, angle = 45, hjust = 1)) + ggtitle("Distribution of TAD sizes")
p2 <- ggplot(domain_numbers, aes(x=factor(sample), y=number)) + geom_bar(stat="identity") + xlab("Sample") + ylab("TAD number") + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=7, angle = 45, hjust = 1)) + scale_y_continuous(limit = c(0, 5000)) + ggtitle("TAD numbers") 
grid.arrange(p1, p2, nrow=2, heights = c(0.55, 0.45))
dev.off()
