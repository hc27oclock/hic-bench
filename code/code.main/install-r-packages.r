#!/usr/bin/env Rscript

write("Installing R packages...", stderr())

packages = c("flsa", "genlasso", "optparse", "ggplot2", "pastecs", "plotrix",
	"reshape2", "zoo", "plyr", "data.table", "gridExtra", "scales", "grid",
	"RColorBrewer", "stringr", "chron", "corrplot")

for (p in packages)
	if (!suppressPackageStartupMessages(require(p, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))) {
		install.packages(p, repos="http://cran.us.r-project.org")
	}

