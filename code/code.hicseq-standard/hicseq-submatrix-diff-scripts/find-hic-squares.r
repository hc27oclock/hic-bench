# install packages
write('Loading libraries...',stderr());
for (p in c("optparse","reshape")) 
  if (!suppressPackageStartupMessages(require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE))) {
    install.packages(p,repos="http://cran.rstudio.com/") 
    library(p,character.only=TRUE,quietly=TRUE,verbose=FALSE)
  }

# setup command-line options
option_list <- list(
  make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
  make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"]."),
  make_option(c("-b","--bin-size"), default=0, help="Bin size in nucleotides for RPKB calculation [default \"%default\"]."),
  make_option(c("-p","--pseudo-count"), default=5, help="Hi-C pseudo-count for fold-change calculation [default \"%default\"]."),
  make_option(c("-s","--square-size"), default=5, help="Size of square to be examined (in bins) [default \"%default\"]."),
  make_option(c("-d","--max-dist"), default=500, help="Maximum distance off diagonal (in bins) [default \"%default\"]."),
  make_option(c("-c","--logfc-cutoff"), default=0.5, help="Absolute log2 fold-change cutoff [default \"%default\"]."),
  make_option(c("-q","--qval-cutoff"), default=0.10, help="false discovery rate cutoff [default \"%default\"].")
)
usage = 'find-hic-squares.r [OPTIONS] SAMPLE-1-DIR SAMPLE-2-DIR'
  
# get command line options (if help option encountered print help and exit)
args <- commandArgs(trailingOnly=T)
arguments <- parse_args(args=args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt <- arguments$options
files <- arguments$args
if (length(files)!=2) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

# input samples
xdir = files[1]   # 'matrices/CUTLL1-DMSO_repeat'
ydir = files[2]   # 'matrices/CUTLL1-JQ1_repeat'

# params
bin_size = opt$"bin-size"         # bin size
p = opt$"pseudo-count"            # pseudo-count
s = opt$"square-size"             # size of square in bins
maxd = opt$"max-dist"             # maximum distance off diagonal
lfc_cutoff = opt$"logfc-cutoff"   # log FC cutoff
qval_cutoff = opt$"qval-cutoff"   # FDR cutoff
outdir = opt$"output-dir"         # output directory
options(scipen=999)

# Create output directory
if (outdir=="") { write('Error: please specify output directory!',stderr()); quit(save='no'); }
if (file.exists(outdir)==FALSE) { dir.create(outdir) } else { write('Warning: output directory already exists!',stderr()) }

# find disrupted sub-TADs
matrices = list.files(xdir,pattern='matrix.*.tsv')
#matrices = 'matrix.chr21.tsv'
out = c()
for (mat in matrices) {
  # write message
  if (opt$verbose) write(paste("Processing ",mat,"...",sep=''),stderr())
  chr = unlist(strsplit(mat,split='[.]'))[2]
      
  # input files
  fx = paste(xdir,mat,sep='/')
  fy = paste(ydir,mat,sep='/')

  # read matrices
  x = as.matrix(read.table(fx,header=T))
  y = as.matrix(read.table(fy,header=T))

  # normalize by sequencing depth
  y = y*sum(x)/sum(y)

  # fold-changes
  lfc = log2((y+p)/(x+p))

  # search for squares with consistent changes
  if (opt$verbose) write("Searching for consistent changes...",stderr())
  X = matrix(0,nrow(x)-s,ncol(x)-s)               # sample #1 values
  Y = matrix(0,nrow(x)-s,ncol(x)-s)               # sample #2 values
  S = matrix(0,nrow(x)-s,ncol(x)-s)               # fold-change scores
  Q = matrix(1,nrow(x)-s,ncol(x)-s)               # q-values
  n = nrow(S)
  for (i in 1:n) {
    I = i:(i+s-1)
    for (j in seq(i,min(n,(i+maxd)),by=1)) {
      J = j:(j+s-1)
      X[i,j] = mean(x[I,J])
      Y[i,j] = mean(y[I,J])
      S[i,j] = mean(lfc[I,J])
      Q[i,j] = t.test(as.vector(x[I,J]),as.vector(y[I,J]),paired=TRUE)$p.value
#      n1 = sum(as.vector(lfc[I,J])>0.1)
#      n2 = sum(as.vector(lfc[I,J])< -0.1)
#      Q[i,j] = binom.test(c(n1,n2),0.5)$p.value
    }
  }

  # convert p-values to FDR
  if (opt$verbose) write("Adjusting p-values...",stderr())
  D = col(Q)-row(Q)
  for (d in seq(0,min(n,maxd),by=1)) Q[D==d] = p.adjust(Q[D==d],method="BH")
  
  # find sub-TAD changes
  if (opt$verbose) write("Generating table...",stderr())
  X1 = melt(X)
  Y1 = melt(Y)
  S1 = melt(S)
  Q1 = melt(Q)
  k = (abs(S1$value)>=lfc_cutoff)&(abs(Q1$value)<=qval_cutoff)
  if (sum(k)>0) {
    SQ = cbind(S1[k,],qval=Q1[k,3],xval=X1[k,3],yval=Y1[k,3])
    D = (SQ[,2]-SQ[,1])*bin_size
    fi = function(i) paste(chr,'+',i*bin_size+1,(i+s)*bin_size)
    fj = function(j) paste(chr,'+',j*bin_size+1,(j+s)*bin_size)
    SQ = cbind(SQ[,3:ncol(SQ)], dist=D, locus1=fi(SQ[,1]), locus2=fj(SQ[,2]))
    colnames(SQ)[1] = 'log2fc'
    out = rbind(out,SQ)
  }
}

# save output bed file
if (opt$verbose) write("Writing output...",stderr())
out$qval = formatC(out$qval,format="e",digits=3,drop0trailing=FALSE)
write.table(out,quote=F,row.names=F,col.names=T,sep='\t',file=paste(outdir,'out.tsv',sep='/'))

# done
if (opt$verbose) write("Done.",stderr())
quit(save='no')


