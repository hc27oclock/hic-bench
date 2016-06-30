#!/bin/bash
# set -x

##
## USAGE: compile_report_wrapper.sh document.Rnw
## This script will compile a .Rnw document with knitr in R, then compile the results .tex file with pdflatex
##

# check script args
if [ $# != 1 ] # if not enough args provided
then
  echo -e "Only one argument is supported for this script at this time" # update this in this in the future, use like '$@' in a for loop
  grep '^##' $0 # print out lines from the script that start with '##' and exit
  exit
fi

# get script arg
tmp_file="$1"

# if 'module' commands work, then load correct R version
if module unload r; then
  echo -e "Loading R module"
  module switch r/3.3.0
  echo -e "Rscript version is:"
  $(Rscript --version)
else
  echo "'module' not available, can't guarantee that R v3.2.3 will be loaded" # update this in the future, call this version directly somehow
  echo -e "Rscript version is:"
  $(Rscript --version)
fi

# compile the document with knitr
# you must have knitr installed already !!
# in the future update to make sure that the arg passed was a file that ends in .Rnw
Rscript --slave --no-save --no-restore - "$1" <<EOF
  ## R code
  cat("\nR loaded\n")
  # get script args
  args <- commandArgs(TRUE)
  cat("Script args are:\n")
  args
  
  # load knitr
  library("knitr")
  
  # compile document
  cat("Compiling document with knitr\n")
  knit(args[1])
EOF


# btw here are some other ways to pass heredoc to R
# R --slave <<EOF
# Rscript --slave --vanilla - <<EOF


# PDF compilation with pdflatex
# check if the .tex file exists; 
# # sometimes knitr can fail but still produce a tex file
if [ -f ${tmp_file%%.Rnw}.tex ]; then
  # check if pdflatex is installed
  if command -v pdflatex &>/dev/null; then
      # gdate "$@"
      echo -e "Compiling the TeX file; file is:\n${tmp_file%%.Rnw}.tex\n\n"
      # compile with pdflatex
      pdflatex "${tmp_file%%.Rnw}.tex" && pdflatex "${tmp_file%%.Rnw}.tex" && pdflatex "${tmp_file%%.Rnw}.tex"
  else
      # date "$@"
      echo -e "pdflatex not installed, unable to compile file!\n"
  fi
else
  echo -e "TeX file not found! File should be:\n${tmp_file%%.Rnw}.tex\n\n"
fi
