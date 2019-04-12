# Create package directory
dir.create(path = file.path('.R', 'lib'), recursive = TRUE)
install.packages('optparse', lib = '.R/lib')
install.packages('reshape2', lib = '.R/lib')
install.packages('plyr', lib = '.R/lib')
install.packages('ggplot2', lib = '.R/lib')
install.packages('scales', lib = '.R/lib')
install.packages('cowplot', lib = '.R/lib')
install.packages('viridis', lib = '.R/lib')
install.packages('RColorBrewer', lib = '.R/lib')
install.packages('ggdendro', lib = '.R/lib')
install.packages('ggrepel', lib = '.R/lib')
install.packages('svglite', lib = '.R/lib')
install.packages("ontologyIndex", lib = '.R/lib')
install.packages("ontologyPlot", lib = '.R/lib')
install.packages("seriation", lib = '.R/lib')

# install Bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite("SummarizedExperiment", lib = '.R/lib', lib.loc = '.R/lib')
biocLite("DESeq2", lib = '.R/lib', lib.loc = '.R/lib')
biocLite("Rgraphviz", lib = '.R/lib', lib.loc = '.R/lib')
biocLite("biomaRt", lib = '.R/lib', lib.loc = '.R/lib')
biocLite("topgo", lib = '.R/lib', lib.loc = '.R/lib')

# install packages from GitHub
install.packages("devtools", lib = '.R/lib')
library('devtools')
install_github('richysix/biovisr', lib='.R/lib')
install_github('richysix/miscr', lib='.R/lib')

