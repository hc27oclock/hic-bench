## ChIP-Seq Automatic Report

This report will automatically parse the `run` file of each pipeline step included in the `pipeline/index.txt` file for information to include in the report.

The report will check for the existence of listed sheets and figures, and include them in a LaTeX beamer style report.

### Report Compilation

#### Step 1.
In the project directory file `project_info.txt`, fill in the required fields with relevant file paths and information (character strings) to be used in the report. This step is required for the report to know where to look for content, and must be completed manually (for the moment). 

For manual report compilation, follow these steps:

#### Step 2a.
Alternatively, run the compile script:

#### Step 2b.
Compile the report RNW file with knitr:

```
$ # in terminal, start R
$ R
> # if knitr is not installed, install it:
> # install.packages("knitr")
> # load knitr
> library("knitr")
> # 'knit' the file
> knit("chipseq_report.Rnw")
...
output file: chipseq_report.tex

[1] "chipseq_report.tex"
> # quit R and return to the terminal
> quit()
```

Compile the resulting TEX file with `pdflatex` 2 or 3 times (yes this is necessary!)

```
$ pdflatex chipseq_report.tex && pdflatex chipseq_report.tex && pdflatex chipseq_report.tex
```

