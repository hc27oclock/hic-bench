#!/usr/bin/Rscript

write("Installing R packages...", stderr())
for (p in c("flsa", "genlasso", "ggplot2", "optparse", "pastecs", "plotrix", "reshape2", "zoo", "plyr", "data.table", "gridExtra", "scales", "RColorBrewer", "grid", "stringr")) 
  if (!suppressPackageStartupMessages(require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE))) {
    install.packages(p, repos='http://cran.us.r-project.org') 
  }



