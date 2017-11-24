#Produce line graph summarising the di-tag size distribution
#Launched by hicup_filter
args <- commandArgs(TRUE)

outdir <- args[1]
file <- args[2]

data <- read.delim(file, header=FALSE, skip=1) 


outputfilename <- basename(file)
outputfilename <- paste(outputfilename, "svg", sep = ".")
outputfilename <- paste(outdir, outputfilename, sep = "")
svg(file=outputfilename)

data <- read.delim(file, header=FALSE, skip=0) 

plot (data, type="l", xlab="Di-tag length (bps)",
     ylab="Frequency", main="Di-tag frequency Vs Length"
    )

garbage <- dev.off()
