####
# caICB correction main function
####

#### set pars
args <- commandArgs(TRUE)
input<-args[1]
input2<-args[2]
outiqr<-as.numeric(args[3]) #remove ourliers with larger than upper_quantile(x) + outiqr*IQR(x)


#### local functions
easy.chrom.order<-function(v, prefix="chr"){
  ####
  # order chrom vector
  ####
  ori<-as.character(unique(v))
  med<-ori
  med<-gsub("X", "10000", med)
  med<-gsub("Y", "10001", med)
  med<-gsub("M", "10002", med)
  odr<-order(as.numeric(gsub(prefix, "", med)))
  factor(v, levels=ori[odr])
}
easy.tab2matix<-function(x, y, z){
  ####
  # Reshape three column data frame to matrix
  # x: rownames; y: colnames; z: counts in matrix
  ####
  out <- matrix(nrow=nlevels(x), ncol=nlevels(y), dimnames=list(levels(x), levels(y)))
  out[cbind(x, y)] <- z
  out
}
easy.format.table<-function(data, tmpfile=tempfile(tmpdir="."), ...){
  ##
  # format data.frame or matrix by output and input
  ##
  write.table(data, file=tmpfile, sep="\t", quote=F)
  out<-read.delim(tmpfile, row.names=1, ...)
  unlink(tmpfile)
  out
}


#### read stats data
stats<-read.delim(input2, header=F)
colnames(stats)<-c("pos", "icb")
stats$chr<-gsub("_.*", "", stats$pos)
stats.sub<-stats


#### read data
dat<-read.delim(input, header=F, stringsAsFactors=F, colClasses=c("character", "character", "integer", "integer", "numeric"))
colnames(dat)<-c("rd1", "rd2", "raw", "dist", "icb")


#### format raw data
ds<-dat
ds$chr<-gsub("_.*", "", ds$rd1)
#ds<-ds[grep("chrY|chrM", ds$chr, invert=T),]


#### ICB correction
ds$norm<-ds$raw/ds$icb


#### correct for distance within each chromosome
## average by bin for ICB normalized reads
g<-split(ds$norm, paste(ds$dist, ds$chr, sep="|||"), drop=T)
g1<-split(ds$icb, paste(ds$dist, ds$chr, sep="|||"), drop=T)
h<-t(sapply(seq_along(g), function(i) {
	x<-g[[i]]
	x1<-g1[[i]]
	outcut<-c(quantile(x, 0.25)-outiqr*IQR(x), quantile(x, 0.75)+outiqr*IQR(x))
	slt <- (x>=outcut[1] & x<=outcut[2])
	y<-x[slt]
	y1<-x1[slt]
	c(mean(y), var(y), length(y), mean(y1))})
)
rownames(h)<-names(g)
dagg2<-data.frame(do.call(rbind, strsplit(rownames(h), "\\|\\|\\|")), h)
dagg2<-easy.format.table(dagg2)
colnames(dagg2)<-c("dist", "chr", "mean", "var", "len", "avebias")
daggnz2<-dagg2[complete.cases(dagg2) & dagg2$mean>0 & dagg2$len>=10,]
## correction
chrom.fac<-easy.chrom.order(daggnz2$chr)
c2dmat<-easy.tab2matix(chrom.fac, factor(daggnz2$dist), daggnz2$mean)
c2cbias<-matrix(NA, nrow=nrow(c2dmat), ncol=nrow(c2dmat))
rownames(c2cbias)<-levels(chrom.fac)
colnames(c2cbias)<-levels(chrom.fac)
for (i in 1:nrow(c2dmat)){
	for (j in i:nrow(c2dmat)){
		fit<-lm(c2dmat[i,] ~ c2dmat[j,] - 1)
		c2cbias[i, j]<-1/as.numeric(fit$coefficients)
		c2cbias[j, i]<-as.numeric(fit$coefficients)
	}
}
c2cbias.1<-c2cbias/c2cbias[,1]
cbias<-apply(c2cbias.1, 2, median, na.rm=T)
cbias[is.na(cbias)]<-0


#### update icb
stats$cbias<-sqrt(cbias[stats$chr])
stats$aicb<-stats$icb*stats$cbias
stats$aicb<-stats$aicb/median(stats$aicb, na.rm=T)
rownames(stats)<-stats$pos


#### write out
write.table(data.frame(sqrt(cbias)), file=paste0(gsub(".icb", "", input2), ".cbias"), sep="\t", quote=F, col.names=F)
write.table(stats[,c("pos", "aicb")], file=paste0(gsub(".icb", "", input2), ".caicb"), sep="\t", quote=F, row.names=F, col.names=F)








