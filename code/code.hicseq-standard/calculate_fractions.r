#!/bin/Rscript

# Load the required libraries
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)

#Rscript code/calculate_percents.R nodups "ICGC_duplications_uniqid.bed" nodels "ICGC_deletions_uniqid.bed"
args <- commandArgs(trailingOnly=TRUE)

# Check if the length of arguments is correct
if(length(args)!=2){print("USAGE: calculate_fractions.r DUPLICATION-FILE DELETION-FILE"); quit(save="no")}

# get the inputs
no_enh <- 54788
dup_file <- as.character(args[1])
no_se  <- 2283
del_file <- as.character(args[2])

# Now read in the boxplot data
codel_enh_data <- read.table("co-deletions-enhancer.boxplot.tsv",header=FALSE,sep="\t")
colnames(codel_enh_data) <- c("Sample","I","II","III","IV","V")
codel_se_data  <- read.table("co-deletions-SE.boxplot.tsv",header=FALSE,sep="\t")
colnames(codel_se_data) <- c("Sample","I","II","III","IV","V")
codup_enh_data <- read.table("co-duplications-enhancer.boxplot.tsv",header=FALSE,sep="\t")
colnames(codup_enh_data) <- c("Sample","I","II","III","IV","V")
codup_se_data  <- read.table("co-duplications-SE.boxplot.tsv",header=FALSE,sep="\t")
colnames(codup_se_data) <- c("Sample","I","II","III","IV","V")

# Convert to long format based on boundary strength
codel_enh_long <- melt(codel_enh_data)
codel_se_long  <- melt(codel_se_data)
codup_enh_long <- melt(codup_enh_data)
codup_se_long  <- melt(codup_se_data)

# Get the columns right
colnames(codel_enh_long) <- c("sample","strength","number")
colnames(codel_se_long) <- c("sample","strength","number")
colnames(codup_enh_long) <- c("sample","strength","number")
colnames(codup_se_long) <- c("sample","strength","number")

# Get the fractions
codel_enh_long$fraction <- codel_enh_long$number/no_enh
codel_se_long$fraction  <- codel_se_long$number/no_se
codup_enh_long$fraction <- codup_enh_long$number/no_enh
codup_se_long$fraction  <- codup_se_long$number/no_se

# Calculate the size of duplications and deletions
dups_file <- read.table(sprintf("%s",dup_file), header=F, sep="\t")
dels_file <- read.table(sprintf("%s",del_file), header=F, sep="\t")

dups_size_gb <- sum(as.numeric(as.character(dups_file[,3] - dups_file[,2])))/10^9
dels_size_gb <- sum(as.numeric(as.character(dels_file[,3] - dels_file[,2])))/10^9

print(dups_size_gb)
print(dels_size_gb)

# Get the final percentages (per Gb)
codel_enh_long$percent <- (codel_enh_long$fraction/dels_size_gb)
codel_se_long$percent  <- (codel_se_long$fraction/dels_size_gb)
codup_enh_long$percent <- (codup_enh_long$fraction/dups_size_gb)
codup_se_long$percent <- (codup_se_long$fraction/dups_size_gb)

# Tag enhancer
codel_enh_long$category <- c("enhancer")
codup_enh_long$category <- c("enhancer") 

# Do the same for superenhancer
codel_se_long$category <- c("SE")
codup_se_long$category <- c("SE")

# Put the data together
deletions <- rbind(codel_enh_long,codel_se_long)
#deletions$category <- as.factor(category)
#deletions$strength <- as.factor(strength)
dups      <- rbind(codup_enh_long,codup_se_long)
#dups$category <- as.factor(category)
#dups$strength <- as.factor(strength)

head(deletions)
head(dups)

# Plot and do the stats
# Get the palette inversed
rev_jco <- rev(pal_jco("default")(2))

# Get the boundary data for deletions
del <- read.table("closest-outside-deletions.boxplot.tsv", header=FALSE, sep="\t")
dup <- read.table("closest-outside-duplications.boxplot.tsv", header=FALSE, sep="\t")

# Get the row counts for deletions and duplications respectively
del_row_sum <- apply(del[,2:6],1,sum)
dup_row_sum <- apply(dup[,2:6],1,sum)

# Get the fractions for deletions and duplications
del_fraction <- del[,2:6]/del_row_sum/dels_size_gb
dup_fraction <- dup[,2:6]/dup_row_sum/dups_size_gb

# Give the names of the samples
del_fraction$sample <- del[,1]
dup_fraction$sample <- dup[,1]

# Put the column names
colnames(del_fraction)[1:5] <- c("I","II","III","IV","V")
colnames(dup_fraction)[1:5] <- c("I","II","III","IV","V")

# Get the melting done
del_frac_long <- melt(del_fraction)
dup_frac_long <- melt(dup_fraction)

# Get the colnames
colnames(del_frac_long) <- c("sample","strength","percent")
colnames(dup_frac_long) <- c("sample","strength","percent")

# Boundary categories close to deletions
pdf("boundaries_close_to_deletions.pdf", useDingbats=FALSE)
my_comparisons <- list( c("I", "V"), c("II", "V"), c("III", "V"), c("IV","V") )
ggboxplot(del_frac_long, x = "strength", y = "percent", fill="strength", palette="Blues", add="jitter") + stat_compare_means(comparisons = my_comparisons) + xlab("Boundary Strength Category\n(I:weakest, V:strongest)") + ylab("Fraction of boundaries per Gb") 
dev.off()

# Boundary categories close to dups
pdf("boundaries_close_to_duplications.pdf", useDingbats=FALSE)
my_comparisons <- list( c("I", "V"), c("II", "V"), c("III", "V"), c("IV","V") )
ggboxplot(dup_frac_long, x = "strength", y = "percent", fill="strength", palette="Reds", add="jitter") + stat_compare_means(comparisons = my_comparisons) + xlab("Boundary Strength Category\n(I:weakest, V:strongest)") + ylab("Fraction of boundaries per Gb") 
dev.off()

# Get thw category levels right
deletions$category <- factor(deletions$category, levels=c("enhancer","SE"))

# deletions
pdf("deletions_paired_true.pdf", useDingbats=FALSE)
p <- ggboxplot(deletions, x = "strength", y = "percent", color = "category", palette = rev_jco, add = "jitter") + xlab("Boundary Strength Category\n(I:weakest, V:strongest)") + ylab("Fraction of enhancer elements per Gb")
p + stat_compare_means(aes(group = category), label = "p.format", method="wilcox.test", paired=TRUE)
print(p)
dev.off()

pdf("deletions_paired_false.pdf", useDingbats=FALSE)
p <- ggboxplot(deletions, x = "strength", y = "percent", color = "category", palette = rev_jco, add = "jitter") + xlab("Boundary Strength Category\n(I:weakest, V:strongest)") + ylab("Fraction of enhancer elements per Gb")
p + stat_compare_means(aes(group = category), label = "p.format", method="wilcox.test", paired=FALSE)
print(p)
dev.off()

# Deletions (version with star)
pdf("deletions_paired_true_star.pdf", useDingbats=FALSE)
p <- ggboxplot(deletions, x = "strength", y = "percent", color = "category", palette = rev_jco, add = "jitter") + xlab("Boundary Strength Category\n(I:weakest, V:strongest)") + ylab("Fraction of enhancer elements per Gb")
p + stat_compare_means(aes(group = category), label = "p.signif", method="wilcox.test", paired=TRUE)
print(p)
dev.off()

pdf("deletions_paired_false_star.pdf", useDingbats=FALSE)
p <- ggboxplot(deletions, x = "strength", y = "percent", color = "category", palette = rev_jco, add = "jitter") + xlab("Boundary Strength Category\n(I:weakest, V:strongest)") + ylab("Fraction of enhancer elements per Gb")
p + stat_compare_means(aes(group = category), label = "p.signif", method="wilcox.test", paired=FALSE)
print(p)
dev.off()

# Get the levels of the category research right
dups$category <- factor(dups$category, levels=c("enhancer","SE"))

# duplications
pdf("duplications_paired_true.pdf", useDingbats=FALSE)
p <- ggboxplot(dups, x = "strength", y = "percent", color = "category", palette = rev_jco, add = "jitter") + xlab("Boundary Strength Category\n(I:weakest, V:strongest)") + ylab("Fraction of enhancer elements per Gb")
p + stat_compare_means(aes(group = category), label = "p.format", method="wilcox.test", paired=TRUE)
print(p)
#ggplot(dups, aes(x=strength, y=percent, fill=category)) + geom_boxplot(outlier.shape=NA) + xlab("Boundary strength") + ylab("Fraction of elements (%)") + theme_bw() + geom_point(position=position_jitterdodge(dodge.width=0.5), size=0.5) + scale_fill_manual(values=c("blue","red"))
dev.off()

pdf("duplications_paired_false.pdf", useDingbats=FALSE)
p <- ggboxplot(dups, x = "strength", y = "percent", color = "category", palette = rev_jco, add = "jitter") + xlab("Boundary Strength Category\n(I:weakest, V:strongest)") + ylab("Fraction of enhancer elements per Gb")
p + stat_compare_means(aes(group = category), label = "p.format", method="wilcox.test", paired=FALSE)
print(p)
#ggplot(dups, aes(x=strength, y=percent, fill=category)) + geom_boxplot(outlier.shape=NA) + xlab("Boundary strength") + ylab("Fraction of elements (%)") + theme_bw() + geom_point(position=position_jitterdodge(dodge.width=0.5), size=0.5) + scale_fill_manual(values=c("blue","red"))
dev.off()

# Now the version with the stars

# duplications
pdf("duplications_paired_true_star.pdf", useDingbats=FALSE)
p <- ggboxplot(dups, x = "strength", y = "percent", color = "category", palette = rev_jco, add = "jitter") + xlab("Boundary Strength Category\n(I:weakest, V:strongest)") + ylab("Fraction of enhancer elements per Gb")
p + stat_compare_means(aes(group = category), label = "p.signif", method="wilcox.test", paired=TRUE)
print(p)
dev.off()

pdf("duplications_paired_false_star.pdf", useDingbats=FALSE)
p <- ggboxplot(dups, x = "strength", y = "percent", color = "category", palette = rev_jco, add = "jitter") + xlab("Boundary Strength Category\n(I:weakest, V:strongest)") + ylab("Fraction of enhancer elements per Gb")
p + stat_compare_means(aes(group = category), label = "p.signif", method="wilcox.test", paired=FALSE)
print(p)
dev.off()
