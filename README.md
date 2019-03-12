# Common and distinct transcriptional signatures of mammalian embryonic lethality

Notes to recreate the analysis in Collins et al. 2019

Each individual step is detailed below. Once the input files are created,
you can run all the R scripts in order by typing `make`

## Setup

To get the packages for R see the [Required R packages](#packages)

Create output and plots directories and ordered sample informtaion file

```
Rscript setup.R data/counts/samples-gt-gender-stage-somites.txt
```

<h3 id="packages">Required R packages</h3>

To install the required packages, if not already installed, run packages_install.R

This installs the packages in a .R/lib directory. If using this set the R_LIBS_USER
variable to ensure the right packages are loaded

```
export R_LIBS_USER=.R/lib
```

The required packages are
```
optparse
reshape2
plyr
ggplot2
scales
cowplot
viridis
RColorBrewer
ggdendro
ggrepel
svglite
ontologyIndex
ontologyPlot
SummarizedExperiment
DESeq2
Rgraphviz
biomaRt
topgo
devtools
richysix/biovisr (GitHub)
richysix/miscr (GitHub)
```
and their dependencies
