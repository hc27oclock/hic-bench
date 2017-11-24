#Produce pie charts summarising the hicup_filter results
#Launched by hicup_truncater
args <- commandArgs(TRUE)
outdir <- args[1]
file <- args[2]

data <- read.delim(file, header=FALSE, skip=1) 

for (i in 1:nrow(data)) {
  line <- data[i,]
  fileInSummary <- line[,1]
  total <-line[,2]
 
  valid<-line[,3]
  noLigation<-line[,8]
  reLigation<-line[,9]
  selfLigation<-line[,10]
  internalRe2<-line[,11]
  unclassified<-line[,12]
  wrongSize<-line[,13]
  unmapped<-line[,14]
  
  percValid <- round( (100 * valid / total), 2)
  
  outputfilename <- paste(outdir, fileInSummary, sep="")
  outputfilename <- paste(fileInSummary, "filter_piechart.svg", sep = ".")
  outputfilename <- paste(outdir, outputfilename, sep = "")
  svg(file=outputfilename)
  
  pieTitle <- paste( "Filter results\n", file, "\nValid ditags: ", valid,
                     "\nPercent valid: ", percValid, "%",  sep = "")
  
  pcData <- c(valid, noLigation, reLigation, selfLigation, internalRe2, unclassified, wrongSize, unmapped)
  percLabels <- round(pcData / total * 100, 1)
  percLabels<- paste(percLabels, "%", sep="")
  
  par(oma=c(4,0,0,0))
  par(mar=c(0,0,2,0))
  
	pie  (pcData, 
    	labels=percLabels,
      main=pieTitle,
      cex.main = 1,
      col = rainbow(7),
      clockwise=TRUE
  )

	par(oma=c(0,0,0,0))
	par(mar=c(0, 0, 0, 0))


	legend("bottom", c("Valid", "No ligation", "Re-ligation", "Self-ligation", "Internal Re2", 
                    "Unclassified", "Wrong size", "Unmapped"),
                   ncol=2, cex=0.8, fill=rainbow(7))

	dev.off()
  
} 



