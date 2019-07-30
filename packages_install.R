# Create package directory
dir.create(path = file.path('.R', 'lib'), recursive = TRUE)
# set mirror
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org/" 
       options(repos=r)
})

packages <- c('optparse', 'reshape2', 'plyr', 'ggplot2', 'scales', 'cowplot',
    'viridis', 'RColorBrewer', 'ggdendro', 'ggrepel', 'svglite',
    'ontologyIndex', 'ontologyPlot', 'dendextend', 'seriation')
for ( package in packages ) {
    if(!require(package, character.only = TRUE)) {
        install.packages(package, lib = '.R/lib')
        require(package, character.only = TRUE, upgrade = FALSE)
    }
}

package <- 'cowplot'
if(!require(package, character.only = TRUE)) {
    install_version("cowplot", "0.9.4", lib = '.R/lib', upgrade = FALSE)
}

# install Bioconductor packages
source("https://bioconductor.org/biocLite.R")
bioc_packages <- c("SummarizedExperiment", "DESeq2", "Rgraphviz", "biomaRt",
                    "topGO")
for ( package in bioc_packages ) {
    if(!require(package, character.only = TRUE)) {
        biocLite(package, lib = '.R/lib', lib.loc = '.R/lib', suppressUpdates = TRUE)
        require(package, character.only = TRUE)
    }
}

# install packages from GitHub
package <- 'devtools'
if(!require(package, character.only = TRUE)) {
    install.packages("devtools", lib = '.R/lib')
    require(package, character.only = TRUE)
}
install_github('richysix/biovisr', lib='.R/lib', upgrade = FALSE)
install_github('richysix/miscr', lib='.R/lib', upgrade = FALSE)

