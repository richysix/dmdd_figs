# Common and distinct transcriptional signatures of mammalian embryonic lethality

Notes to recreate the analysis in Collins et al. 2019

Each individual step is detailed below.

## Setup

To get the packages for R see the [Required R packages](#packages)

Create output and plots directories and ordered sample informtaion file

```
Rscript setup.R data/counts/samples-gt-gender-stage-somites.txt \
data/Mm_GRCm38_e88_baseline.rds
```

## Suppl Fig 1

The data for the heatmap in Supplementary Figure 1 are provided in
data/PRJEB4513-E8.25/. The heatmap can be created by running the fig1.R script.

```
Rscript fig1.R \
data/PRJEB4513-E8.25/4567_somites-counts.tsv \
data/PRJEB4513-E8.25/4567_somites-samples.tsv \
data/PRJEB4513-E8.25/tissues-counts.tsv \
data/PRJEB4513-E8.25/tissues-samples.tsv
```

Download counts for all 73 lines
```
# The quickest way is to download the compressed archive
# You'll need tar to extract the file
# Alternatively the files can be download from doi.org/10.6084/m9.figshare.6819611
curl -LO https://ndownloader.figshare.com/files/14151989
mv 14151989 data/collins_2019_counts_all.tgz
cd data
tar -xzvf collins_2019_counts_all.tgz
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
