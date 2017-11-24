#Produce pie charts summarising the hicup_deduplicator results
#Launched by hicup_deduplicator
args <- commandArgs(TRUE)
outdir <- args[1]
file <- args[2]

data <- read.delim(file, header=FALSE, skip=1) 

for (i in 1:nrow(data)) {
  line <- data[i,]
  fileInSummary <- line[,1]
  total <- line[,2]
  uniques <- line[,3]

  percUniques <- round( (100 * uniques / total), 2 )
  
  outputfilename <- paste(outdir, fileInSummary, sep="")
  outputfilename <- paste(outputfilename,"deduplicator_uniques_barchart.svg", sep = ".")
  svg(file=outputfilename)
  
  graphTitle <- paste( fileInSummary, "\nDitag De-duplication results","\nUnique di-tags: ", percUniques, "%", sep = "")

  bpData <- c(total, uniques) 
  bp <- barplot(bpData,
          names.arg = c("Before de-duplication", "After de-duplication"),
          main = graphTitle,
          col=rainbow(2),
  )

   text( bp, bpData, bpData, cex=1, pos=1) 
 
  garbage <- dev.off()
  
} 



