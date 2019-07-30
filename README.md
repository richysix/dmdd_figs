# Common and distinct transcriptional signatures of mammalian embryonic lethality

Notes to recreate the analysis in Collins et al. 2019

These notes start from the processed count files. These were produced by mapping
the reads to GRCm38 using TopHat2 and counted against Ensembl version 88 annotation
using htseq-count.

Each individual step is detailed below.

## Setup

To get the packages for R see the [Required R packages](#packages) section

Download count files for all 73 lines
```
# The quickest way is to download the compressed archive
# You'll need tar to extract the file
# Alternatively the files can be download from doi.org/10.6084/m9.figshare.6819611
curl -L --output data/collins_2019_counts_all.tgz https://ndownloader.figshare.com/files/14151989
cd data
tar -xzvf collins_2019_counts_all.tgz
cd ..
```

This creates a data/counts directory with the counts for each line and a samples
file detailing the genotype, sex and stage of each sample.

Create output and plots directories and ordered sample information file

```
Rscript setup.R data/counts/samples-gt-gender-stage-somites.txt \
data/Mm_GRCm38_e88_baseline.rds
```

## Fig. 1
### Fig. 1d

The data for the Novel genes heatmap in Fig. 1d are included in the data/ folder
(data/346-from-all-novel-blacklist-header.tsv).
The heatmap was generated using this datafile and the baseline samples files
`output/samples-Mm_GRCm38_e88_baseline.txt` using the
[geneExpr](https://richysix.shinyapps.io/geneexpr/) Shiny App.
The `Max Scaled` Transformation option was used and the heatmap was clustered by
Genes. Five of the genes have zero variance in their counts and so will produce
the error below when clustering. This removes these genes from the heatmap.
```
Clustering error
Some of the genes that you are trying to cluster have zero variance across the
selected samples and have been removed:
XLOC_012085, XLOC_014947, XLOC_018956, XLOC_020760, XLOC_044545
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

## Fig. 2 and Fig. 3 

The commands to create the files needed by the fig2 script are in fig2.md


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
