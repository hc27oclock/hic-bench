#R

# requires rds file on logFC between any two conditions.
calculateTADDifferences <- function(conn_table_sample_1_norm, conn_table_sample_2_norm, out.prefix, tads.x, tads.y, bin.size, chr, n=NA) {
	print("Calculate differences between TADs dervied from two hi-c experiments.")

	result_pvalue <- vector()
	result_mean_x <- vector()
	result_mean_y <- vector()
	bin_range_i <- vector()
	bin_range_j <- vector()
	list_i <- vector()
	list_j <- vector()

	for (k in 1:nrow(tads.x)) {

		if (tads.x[k,1] != chr) {
			next
		}

#		tad.start <- floor(mean(tads.x[i,2], tads.y[i,2]))
#		tad.end <- ceiling(mean(tads.x[i,3], tads.y[i,3]))
		tad.start <- min(tads.x[k,2], tads.y[k,2])
		tad.end <- max(tads.x[k,3], tads.y[k,3])

		i <- tad.start / bin.size
		j <- tad.end / bin.size
		print(paste("i=",i," and j=",j,sep=""))

#		print(paste("start=",tad.start," and end=",tad.end,sep=""))

		sub_matrix_sample_1 <- conn_table_sample_1_norm[i:j, i:j]
		sub_matrix_sample_2 <- conn_table_sample_2_norm[i:j, i:j]


#		if (all(is.na(sub_matrix_sample_1[,apply(sub_matrix_sample_1, 2, var, na.rm=TRUE) != 0]),na.rm=TRUE) | all(is.na(sub_matrix_sample_2[,apply(sub_matrix_sample_2, 2, var, na.rm=TRUE) != 0]),na.rm=TRUE)) {
#			result_pvalue <- rbind(result_pvalue, 1)
#			result_mean_x <- rbind(result_mean_x, 0)
#			result_mean_y <- rbind(result_mean_y, 0)
#			bin_range_i <- rbind(bin_range_i, tad.start)
#			bin_range_j <- rbind(bin_range_j, tad.end)
#			list_i <- rbind(list_i, i)
#			list_j <- rbind(list_j, j)
#			next
#		}

#		tryCatch(t.test(unlist(sub_matrix_sample_1),unlist(sub_matrix_sample_2),paired=TRUE), error=function(x) print(x) ))

#		obj<-try(t.test(unlist(sub_matrix_sample_1), unlist(sub_matrix_sample_2),paired=TRUE), silent=TRUE)

#		result <- tryCatch(t.test(unlist(sub_matrix_sample_1),unlist(sub_matrix_sample_2), paired=TRUE, alternative="two.sided", na.action="na.omit"), error=function(e) NA)
		result <- try(t.test(as.vector(unlist(sub_matrix_sample_1)),as.vector(unlist(sub_matrix_sample_2)), paired=TRUE, alternative="two.sided", na.action="na.omit"),silent=TRUE)

		if (is(result, "try-error")) { 
			result_pvalue <- rbind(result_pvalue, 1)
			result_mean_x <- rbind(result_mean_x, 0)
			result_mean_y <- rbind(result_mean_y, 0)
			bin_range_i <- rbind(bin_range_i, tad.start)
			bin_range_j <- rbind(bin_range_j, tad.end)
			list_i <- rbind(list_i, i)
			list_j <- rbind(list_j, j)
			next
		}
		#NA else obj$p.value

		if (is.na(result$p.value)) {
			result_pvalue <- rbind(result_pvalue, 1)
			result_mean_x <- rbind(result_mean_x, 0)
			result_mean_y <- rbind(result_mean_y, 0)
			bin_range_i <- rbind(bin_range_i, tad.start)
			bin_range_j <- rbind(bin_range_j, tad.end)
			list_i <- rbind(list_i, i)
			list_j <- rbind(list_j, j)
		} else {
			result_pvalue <- rbind(result_pvalue, result$p.value)
			result_mean_x <- rbind(result_mean_x, mean(unlist(sub_matrix_sample_1),na.rm=TRUE))
			result_mean_y <- rbind(result_mean_y, mean(unlist(sub_matrix_sample_2),na.rm=TRUE))
			bin_range_i <- rbind(bin_range_i, tad.start)
			bin_range_j <- rbind(bin_range_j, tad.end)
			list_i <- rbind(list_i, i)
			list_j <- rbind(list_j, j)
		}
#		if ((j-i) > (max.range / bin.size)) {
#			break
#		}
	}

	print(paste("n=",n, " and length pval=",length(result_pvalue),sep=""))

	if (length(result_pvalue) == 0) {
		print("No results created, skipping saving results step...")
		return(NA)
	}

	if (is.na(n)) {
		n <- length(result_pvalue)
	}
	result_fdr <- p.adjust(result_pvalue, method="fdr", n=as.numeric(n))
	print("fesfes")

	logFC <- rep(NA,length(result_mean_x))
	for (z in 1:length(result_mean_x)) {
		if (result_mean_x[z] > 0 & result_mean_y[z] > 0) {
			logFC[z] <- (log2(result_mean_x[z] / result_mean_y[z]))
		}
	}

	results <- data.frame(row.names=list_i,rep(chr, length(list_i)),cbind(list_i,bin_range_i, list_j, bin_range_j,result_mean_x,result_mean_y,(result_mean_x-result_mean_y),
		logFC,result_pvalue,result_fdr))
	colnames(results) <- c("chr","i","TAD_start","j","TAD_end","sample_1_mean","sample_2_mean","mean_diff","logFC","pval","FDR")

	print("Saving results...")
	write.table(file=paste(out.prefix,"_results-table.tsv",sep=""),results,quote=FALSE, row.names=FALSE, sep="\t")
#	saveRDS(file=paste(out.prefix,"_results-table.rds",sep=""),results)
	print("Done.")
	return(results)
}


zscore <- function(x) {
	return((x - mean(x,na.rm=TRUE)) / sd(x,na.rm=TRUE))
}

norm_score3 <- function(x) { 
	return(((x / sum(x,na.rm=TRUE)) * 1000))
}

norm_score <- function(x) {
	results <- vector()
	x_wo_zero <- x
	x_wo_zero[x_wo_zero==0] <- NA
	x_median <- median(x_wo_zero,na.rm=TRUE)
	x_mad <- mad(x_wo_zero,na.rm=TRUE,constant=1.486)
	for (y in x) {
		if (is.na(y)) {
			results <- c(results,NA)
		} else if (y>0) {
			results <- c(results,(y-x_median)/x_mad)
		} else {
			results <- c(results,(y-m_median)/m_mad)
		}
	}
	return(results)
}


norm_score2 <- function(x) {
	return(((x - median(x,na.rm=TRUE)) / mad(x,na.rm=TRUE)))
}

#library(gplots)
#library(RColorBrewer)
#library(Hmisc)
#library(edgeR)

#source("differential_tad_activity_functions.r")

args <- commandArgs()

conn_table_file_sample_1 <- args[3]
conn_table_file_sample_2 <- args[4]
tads_file_sample_1 <- args[5]
tads_file_sample_2 <- args[6]
chr <- args[7]
n <- as.numeric(args[8])
is.normalize <- args[9]
bin.size <- as.numeric(args[10])
centrotelo_file <- args[11]
max_tad_size <- as.numeric(args[12])
out.prefix <- args[13]

#is.normalize <- FALSE

print(paste("Run R script on files ", conn_table_file_sample_1, ", ", conn_table_file_sample_2, ", ", out.prefix,", on chromosome ", chr, " and normalization=",is.normalize,sep=""))

#bin.n <- 10
#bin.m <- 10
#max.range <- 2000000
#bin.size <- 40000

conn_table_sample_1 <- as.matrix(read.csv(file=conn_table_file_sample_1, sep="\t", header=TRUE, row.names=1))
conn_table_sample_2 <- as.matrix(read.csv(file=conn_table_file_sample_2, sep="\t", header=TRUE, row.names=1))

# INCLUDING DIAGONAL VALUES
#for (i in 2:(ncol(conn_table_sample_1))) {
#	for (j in 1:(i-1)) {
#		conn_table_sample_1[i,j] <- NA
#		conn_table_sample_2[i,j] <- NA
#	}
#}

# EXCLUDING DIAGONAL - I would always do so, no matter what...

# Set non-mappable regions to NA

for (i in 2:(ncol(conn_table_sample_1)-1)) {
	for (j in 1:(i+1)) {
		conn_table_sample_1[i,j] <- NA
		conn_table_sample_2[i,j] <- NA
	}
}

centrotelo <- read.table(centrotelo_file, header = FALSE, sep = "\t")
centrotelo <- centrotelo[centrotelo$V1 == chr,]

for (i in 1:nrow(centrotelo)) {
        centrotelo_entry <- centrotelo[i,]
        centrotelo_entry_start_bin <- as.numeric(floor(centrotelo_entry[2] / bin.size)) + 1
        centrotelo_entry_end_bin <- as.numeric(ceiling(centrotelo_entry[3] / bin.size)) + 1
 	if (i == 3){
		centrotelo_entry_end_bin <- ncol(conn_table_sample_1)
	}
	conn_table_sample_1[centrotelo_entry_start_bin:centrotelo_entry_end_bin,] <- NA
        conn_table_sample_1[centrotelo_entry_start_bin:centrotelo_entry_end_bin,] <- NA
        
	conn_table_sample_2[,centrotelo_entry_start_bin:centrotelo_entry_end_bin] <- NA
        conn_table_sample_2[,centrotelo_entry_start_bin:centrotelo_entry_end_bin] <- NA
}
 
conn_table_sample_1[rowSums(conn_table_sample_1,na.rm=TRUE) == 0,] <- NA
conn_table_sample_1[rowSums(conn_table_sample_1,na.rm=TRUE) == 0,] <- NA

conn_table_sample_2[colSums(conn_table_sample_2,na.rm=TRUE) == 0,] <- NA
conn_table_sample_2[colSums(conn_table_sample_2,na.rm=TRUE) == 0,] <- NA

if (is.normalize) {
#	conn_table_sample_1_norm <- norm_score(conn_table_sample_1)
#	conn_table_sample_2_norm <- norm_score(conn_table_sample_2)

#	m_median_1 <- median(as.matrix(conn_table_sample_1),na.rm=TRUE)
#	m_mad_1 <- mad(as.matrix(conn_table_sample_1),na.rm=TRUE)
#	m_median_2 <- median(as.matrix(conn_table_sample_2),na.rm=TRUE)
#	m_mad_2 <- mad(as.matrix(conn_table_sample_2),na.rm=TRUE)

	# Initialize normalized matrices
	conn_table_sample_1_norm <- conn_table_sample_1
	conn_table_sample_2_norm <- conn_table_sample_2

	# Apply normalization to the diagonals of the matrices
	for (distance in 2:(ncol(conn_table_sample_1)-1)){
		vector_values_1 <- vector()
		vector_values_2 <- vector()

		for (i in 1:(ncol(conn_table_sample_1)-1-distance)) {
			vector_values_1 <- c(vector_values_1,conn_table_sample_1[i,i+distance])
			vector_values_2 <- c(vector_values_2,conn_table_sample_2[i,i+distance])
		}
#		m_median <- m_median_1
#		m_mad <- m_mad_1
		vector_values_norm_1 <- norm_score3(vector_values_1)
#		m_median <- m_median_2
#		m_mad <- m_mad_2
		vector_values_norm_2 <- norm_score3(vector_values_2)

		# Copy the normalized values to the normalized matrix
		for (i in 1:(ncol(conn_table_sample_1)-1-distance)){
			conn_table_sample_1_norm[i,i+distance] <- vector_values_norm_1[i]
			conn_table_sample_2_norm[i,i+distance] <- vector_values_norm_2[i]
		}		
	}

} else {
	conn_table_sample_1_norm <- conn_table_sample_1
	conn_table_sample_2_norm <- conn_table_sample_2
}

#saveRDS(file=paste(out.prefix,"_sample_1_norm.rds",sep=""), conn_table_sample_1_norm)
#saveRDS(file=paste(out.prefix,"_sample_2_norm.rds",sep=""), conn_table_sample_2_norm)

#conn_table_sample_1_norm <- readRDS(file=conn_table_file_sample_1)
#conn_table_sample_2_norm <- readRDS(file=conn_table_file_sample_2)

#conn_table_logFC <- log2(conn_table_sample_1_norm/conn_table_sample_2_norm)

#saveRDS(file=paste(out.prefix,"_logFC.rds",sep=""), conn_table_logFC)

tads.sample_1 <- read.csv(file=tads_file_sample_1, sep="\t", header=FALSE)
tads.sample_2 <- read.csv(file=tads_file_sample_2, sep="\t", header=FALSE)

calculateTADDifferences(conn_table_sample_1_norm=conn_table_sample_1_norm,conn_table_sample_2_norm=conn_table_sample_2_norm, 
	out.prefix=out.prefix, bin.size=bin.size, tads.x=tads.sample_1, tads.y=tads.sample_2, chr=chr,n=n)

#stop(0)
























