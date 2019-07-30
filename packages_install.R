# Create package directory
package_dir <- file.path('.R', 'lib')
if (!file.exists(package_dir)) {
    dir.create(path = package_dir, recursive = TRUE)
}

# set mirror
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org/" 
       options(repos=r)
})

packages <- c('optparse', 'reshape2', 'plyr', 'ggplot2', 'scales', 'cowplot',
    'viridis', 'RColorBrewer', 'ggdendro', 'ggrepel', 'svglite',
    'ontologyIndex', 'ontologyPlot', 'dendextend', 'seriation', 'devtools')
not_installed_packages <- c()
index <- 1
for ( package in packages ) {
    if(!require(package, character.only = TRUE)) {
        not_installed_packages[index] <- package
        index = index + 1
    }
}
install.packages(not_installed_packages, lib = '.R/lib')
for ( package in not_installed_packages ) {
    require(package, character.only = TRUE)
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
install_github('richysix/biovisr', lib='.R/lib', upgrade = FALSE)
install_github('richysix/miscr', lib='.R/lib', upgrade = FALSE)

