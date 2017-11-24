#R

args <- commandArgs()

tad_activity_file <- args[3]
gene_bins_file <- args[4]
exprs_file_case_1 <- args[5]
exprs_file_case_2 <- args[6]
case_left <- args[7]
case_right <- args[8]
bin.size <-  as.numeric(args[9])
min.TAD.size <- as.numeric(args[10])
out_prefix <- args[11]

tad_activity <- read.csv(file=tad_activity_file,header=TRUE,sep="\t")
# min-TAD size requirement
tad_activity <- tad_activity[(tad_activity$j - tad_activity$i) >= (min.TAD.size/bin.size),]
tad_activity <- tad_activity[tad_activity$logFC != Inf & tad_activity$logFC != -Inf ,]
gene_bins <- read.csv(file=gene_bins_file, sep="\t", header=TRUE)

# necessary??
gene_bins <- gene_bins[gene_bins$type=="protein_coding" | gene_bins$type=="processed_transcript",]

tad_activity_up <- (tad_activity[tad_activity$FDR <= 0.1 & tad_activity$logFC >= 0.75,])
tad_activity_down <- (tad_activity[tad_activity$FDR <= 0.1 & tad_activity$logFC <= -0.75,])
tad_activity_unchanged <- (tad_activity[tad_activity$FDR >= 0.1 | abs(tad_activity$logFC) < 0.75 ,])

write.table(file=paste(out_prefix,"_active-TADs.bed",sep=""),tad_activity_up[,c(1,3,5)],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(file=paste(out_prefix,"_inactive-TADs.bed",sep=""),tad_activity_down[,c(1,3,5)],sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(file=paste(out_prefix,"_unchanged-TADs.bed",sep=""),tad_activity_unchanged,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)

write.table(file=paste(out_prefix,"_active-TADs.tsv",sep=""),tad_activity_up,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(file=paste(out_prefix,"_inactive-TADs.tsv",sep=""),tad_activity_down,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)


png(file=paste(out_prefix, "_volcano_TADs.png", sep=""),width=2048, height=2048, pointsize=110)
#max_x <- max(abs(tad_activity$logFC),na.rm=TRUE)
par(mar=c(4.2,4,0.4,0.4))
plot(tad_activity$logFC, -log(tad_activity$FDR), xlab="log(fold-change)",ylab="-log(FDR)",pch=20,col="grey",
	cex=1.5,bty="n")
#XX
#max_x <- 3.5
#max_y <- 175
#plot(tad_activity$logFC, -log(tad_activity$FDR), xlab="log(fold-change)",ylab="-log(FDR)",pch=20,col="grey",
#	cex=1.5,xlim=c(((-1)*max_x),max_x), ylim=c(0,max_y),bty="n")

#font=2,font.lab=2,
box(lwd=8)
axis(side = 1, lwd = 8)
axis(side = 2, lwd = 8)
points(tad_activity_up$logFC, -log(tad_activity_up$FDR),col="darkred",pch=20,cex=1.5)
points(tad_activity_down$logFC, -log(tad_activity_down$FDR),col="blue",pch=20,cex=1.5)
abline(h=c(2.302585), col="red", lwd=8, lty=2)
abline(v=-0.75,col="red",lwd=8,lty=2)
abline(v=0.75,col="red",lwd=8,lty=2)
dev.off()

png(file=paste(out_prefix, "_volcano_TADs_legend.png", sep=""),width=2048, height=2048, pointsize=110)
plot.new()
legend("center",legend=c("Higher activity in T-ALL", "Less activity in T-ALL"),fill=c("darkred","blue"),lwd=8)
dev.off()

#tad_extension <- 0
tad_extension <- 120000

#gene_bin_hits_exprs_up <- data.frame()
gene_hits <- vector()
for (i in 1:nrow(tad_activity)) {
	gene_bin_hits <- gene_bins[as.character(gene_bins$chr) == as.character(tad_activity[i,]$chr) & gene_bins$start >= (tad_activity[i,]$TAD_start - tad_extension) & 
		gene_bins$end <= (tad_activity[i,]$TAD_end + tad_extension),]
	if (nrow(gene_bin_hits)==0) {
		gene_hits <- rbind(gene_hits,"NA")
	} else {
		gene_hits <- rbind(gene_hits,paste(as.character(gene_bin_hits$symbol),collapse=","))
	}
#		gene_hits_up <- c(gene_hits_up, as.character(gene_bin_hits$ID))
#		gene_bin_hits_exprs <- rbind(gene_bin_hits_exprs, exprs[rownames(exprs) %in% gene_bin_hits$ID,])
}
tad_activity <- cbind(tad_activity,gene_hits)
write.table(file=paste(out_prefix,"_annotated.tsv",sep=""),tad_activity,sep="\t",quote=FALSE,row.names=FALSE)

if (exprs_file_case_1 == "FALSE") {
	stop(0)
}

exprs_case_1 <- as.data.frame(readRDS(exprs_file_case_1))
exprs_case_2 <- as.data.frame(readRDS(exprs_file_case_2))

gene_bin_hits_exprs_up <- data.frame()
gene_hits_up <- vector()
for (i in 1:nrow(tad_activity_up)) {
		gene_bin_hits <- gene_bins[as.character(gene_bins$chr) == as.character(tad_activity_up[i,]$chr) & gene_bins$start >= (tad_activity_up[i,]$TAD_start - tad_extension) & 
			gene_bins$end <= (tad_activity_up[i,]$TAD_end + tad_extension),]
		gene_hits_up <- c(gene_hits_up, as.character(gene_bin_hits$ID))
#		gene_bin_hits_exprs_up <- rbind(gene_bin_hits_exprs_up, exprs[rownames(exprs) %in% gene_bin_hits$ID,])
		
}

gene_bin_hits_exprs_down <- data.frame()
gene_hits_down <- vector()
for (i in 1:nrow(tad_activity_down)) {
		gene_bin_hits <- gene_bins[as.character(gene_bins$chr) == as.character(tad_activity_down[i,]$chr) & gene_bins$start >= (tad_activity_down[i,]$TAD_start - tad_extension) & 
			gene_bins$end <= (tad_activity_down[i,]$TAD_end + tad_extension),]
		gene_hits_down <- c(gene_hits_down, as.character(gene_bin_hits$ID))
#		gene_bin_hits_exprs_down <- rbind(gene_bin_hits_exprs_down, as.data.frame(exprs[rownames(exprs) %in% gene_bin_hits$ID,]))
}

gene_hits_up <- na.omit(gene_hits_up)
gene_hits_down <- na.omit(gene_hits_down)

#ERG
#tad_activity[tad_activity$chr=="chr21" & tad_activity$TAD_start < 40000001 & tad_activity$TAD_end > 40040000,]

write.table(file=paste(out_prefix,"_genes-in-active-TADs.tsv",sep=""),(gene_hits_up),sep="\t",quote=FALSE,row.names=FALSE)
write.table(file=paste(out_prefix,"_genes-in-inactive-TADs.tsv",sep=""),(gene_hits_down),sep="\t",quote=FALSE,row.names=FALSE)

#rpkm_cd4_up <- data.frame(rpkm_cd4_up=apply(gene_bin_hits_exprs_up[,4:16],1,mean), name=rownames(gene_bin_hits_exprs_up))
#rpkm_cd4_down <- data.frame(rpkm_cd4_down=apply(gene_bin_hits_exprs_down[,4:16],1,mean), name=rownames(gene_bin_hits_exprs_down))

rpkm_left_up <- exprs_case_1[rownames(exprs_case_1) %in% gene_hits_up,]
rpkm_left_down <- exprs_case_1[exprs_case_1$name %in% gene_hits_down,]
rpkm_right_up <- exprs_case_2[exprs_case_2$name %in% gene_hits_up,]
rpkm_right_down <- exprs_case_2[exprs_case_2$name %in% gene_hits_down,]

#rpkm_cd4_up <- data.frame(rpkm_cd4_up=apply(gene_bin_hits_exprs_up$rpkm_mean_cd4,1,mean), name=rownames(gene_bin_hits_exprs_up))
#rpkm_cd4_down <- data.frame(rpkm_cd4_down=apply(gene_bin_hits_exprs_down$rpkm_mean_cd4,1,mean), name=rownames(gene_bin_hits_exprs_down))
#rpkm_tall_up <- data.frame(rpkm_tall_up=gene_bin_hits_exprs_up[,case_left], name=rownames(gene_bin_hits_exprs_up))
#rpkm_tall_down <- data.frame(rpkm_tall_down=gene_bin_hits_exprs_down[,case_left], name=rownames(gene_bin_hits_exprs_down))

logFC_up <- log2(rpkm_right_up$rpkm / rpkm_left_up$rpkm)
logFC_down <- log2(rpkm_right_down$rpkm / rpkm_left_down$rpkm)

rpkm_right_up$logFC  <- log2(rpkm_right_up$rpkm / rpkm_left_up$rpkm)
rpkm_right_down$logFC <- log2(rpkm_right_down$rpkm / rpkm_left_down$rpkm)

write.table(file=paste(out_prefix,"_genes-in-active-TADs-logFC.tsv",sep=""),cbind(as.character(rpkm_right_up$name),rpkm_right_up$logFC),sep="\t",quote=FALSE,row.names=FALSE)
write.table(file=paste(out_prefix,"_genes-in-inactive-TADs-logFC.tsv",sep=""),cbind(as.character(rpkm_right_down$name),rpkm_right_down$logFC),sep="\t",quote=FALSE,row.names=FALSE)

logFC_up_cleaned <- logFC_up[logFC_up != Inf & logFC_up != -Inf]
logFC_up_cleaned_test <- t.test(logFC_up_cleaned, na.action="na.omit", alternative="greater")
logFC_down_cleaned <- logFC_down[logFC_down != Inf & logFC_down != -Inf]
logFC_down_cleaned_test <- t.test(logFC_down_cleaned, na.action="na.omit", alternative="less")
test_results <- rbind(cbind(logFC_up_cleaned_test$estimate, logFC_up_cleaned_test$p.value),cbind(logFC_down_cleaned_test$estimate, logFC_down_cleaned_test$p.value))
rownames(test_results) <- c("Higher TAD activity", "Less TAD activity")
colnames(test_results) <- c("Mean differences from 0", "p-value")
write.table(file=paste(out_prefix,"_logFC-differences_pvalues.tsv",sep=""),test_results,quote=FALSE, row.names=TRUE, sep="\t")


png(file=paste(out_prefix, "_boxplot_RNAexprs_logFC.png",sep=""),width=1200,height=2048,pointsize=110)
box_names <- c(case_right,case_left)
#par(mfrow=c(1,1),mar=c(8,4,4,4))
#y_max <- max(c(rpkm_left_up,gene_bin_hits_exprs_up[,tall_exprs_name],rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name]),
#	na.rm=TRUE)
#y_min <- min(c(rpkm_left_up,gene_bin_hits_exprs_up[,tall_exprs_name],rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name]),
#	na.rm=TRUE)
y_min <- -11
y_max <- 12
par(mar=c(0.1,4,0.1,0.1))
boxplot(logFC_down, logFC_up, beside=TRUE,  notch=TRUE,col=c("blue","darkred"),las=2,
	ylab="RNA logFC",na.action="na.omit", ylim=c(y_min,y_max),byt="n",pch=20,font=2,lwd=8,outline=FALSE)
#boxplot(logFC_down, logFC_up, beside=TRUE, notch=TRUE,col=c("blue","darkred"),na.action="na.omit", add=TRUE,las=2, ylim=c(y_min,y_max),bty="n",pch=20,lwd=8)
abline(h=0, lwd=10,col="black")
stripchart(list(logFC_down=logFC_down,logFC_up=logFC_up),
        col=c("blue","darkred"),add=TRUE,pch=20,cex=1.4,vertical=TRUE,method="jitter")
boxplot(logFC_down, logFC_up, beside=TRUE, notch=TRUE,col=c("blue","darkred"),na.action="na.omit", add=TRUE,las=2, ylim=c(y_min,y_max),bty="n",pch=20,lwd=8,outline=FALSE)
box(lwd=8)
axis(side = 2, lwd = 8,labels=FALSE)
#legend("topright", legend=c("Less TAD activity","Higher TAD activity"), fill=c("blue","darkred"))
#boxplot(rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name],beside=TRUE, main="Less TAD activity in T-ALL",notch=TRUE,col="lightgray",las=2,
#	names=box_names,ylab="logFPKM",ylim=c(y_min,y_max))
dev.off()



