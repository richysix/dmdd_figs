#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option("--output_base", type="character", default='mrna_abnormal-PC3-jaccard',
              help="Basename for output files [default %default]" ),
  make_option("--cluster_methods", type="character", default='ward.D2',
              help="a comma-separated list of cluster methods to use [default %default]" ),
  make_option("--output_data_file", type="character", default='ward.D2',
              help="Save an rda file of the clustering object and the heatmap. Only the last of any clustering methods will be saved [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'cluster_by_overlap.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 1
)

#cmd_line_args <- list(
#  options = list(output_base = 'mrna_abnormal-PC3-jaccard',
#                 cluster_methods = 'ward.D2',
#                 verbose = FALSE ),
#  args = c('output/mrna_abnormal-hom_vs_het_wt-sig_genes.out')
#)

################################################################################
# THIS SCRIPT ASSUMES THAT THE NAMES OF THE THINGS BEING OVERLAPPED
# ARE IN COLUMNS 1 AND 2 AND THE AMOUNT OF OVERLAP IS IN COLUMN 6
################################################################################

output_base <- cmd_line_args$options[['output_base']]
verbose <- cmd_line_args$options[['verbose']]
cluster_methods <- unlist(strsplit(cmd_line_args$options[['cluster_methods']],
                          ",", fixed = TRUE))
output_data_file <- cmd_line_args$options[['output_data_file']]

if( verbose ){
  cat( "Output files basename:", output_base, "\n", sep=" " )
  cat( "Cluster methods:", paste(cluster_methods, sep = '', collapse =', '), "\n", sep=" " )
}

packages <- c('ggplot2', 'reshape2', 'viridis', 'seriation')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# input data
input_data <- read.table(file = cmd_line_args$args[1], sep = "\t",
                         header = FALSE)

# create distance matrix
genes <- unique( c(levels(input_data[[1]]), levels(input_data[[2]])) )
jaccard_overlaps <- numeric(length = length(genes)^2)
index <- 1
for (i in seq_len(length(genes))) {
  for (j in seq_len(length(genes))) {
    if (i == j) {
      overlap <- 1
    } else {
      overlap <- input_data[ input_data[[1]] == genes[i] & input_data[[2]] == genes[j], 6 ]
    }
    if (length(overlap) == 0) {
      overlap <- input_data[ input_data[[2]] == genes[i] & input_data[[1]] == genes[j], 6 ]
    }
    if (length(overlap) == 0) {
      stop(printf('Could not find overlap for %s and %s', genes[i], genes[j]))
    }
    jaccard_overlaps[index] <- overlap
    index <- index + 1
  }
}

jaccard_overlaps <- matrix(jaccard_overlaps, nrow = length(genes))
dimnames(jaccard_overlaps) <- list(genes, genes)

jaccard_dist <- as.dist(1 - jaccard_overlaps)

matrix_heatmap <- function(molten_df, x = 'Var1', y = 'Var2', fill = 'value') {
  heatmap_plot <- ggplot(data = molten_df) +
    geom_raster( aes_(x = as.name(x), y = as.name(y),
                       fill = as.name(fill) ) ) +
    scale_x_discrete(position = 'top') +
    scale_fill_viridis(name = 'Jaccard Overlap', na.value = 'white') +
    labs(x = 'Gene KO', y = 'Gene KO') +
    theme_void() +
    theme( axis.text.x = element_text(colour = 'black', size = 12, angle = 45,
                                      hjust = 0, vjust = 0.5, debug = FALSE),
          axis.text.y = element_text(colour = 'black', size = 12,
                                     angle = 0, debug = FALSE),
          legend.title = element_text(colour = 'black', size = 14,
                                     angle = 0, debug = FALSE),
          legend.text = element_text(colour = 'black', size = 12,
                                     angle = 0, debug = FALSE),
          legend.position = c(0.75,0.75),
          plot.margin = unit(c(5.5, 11, 5.5, 5.5), 'pt'))
  return(heatmap_plot)
}

# cluster
for (cluster_method in cluster_methods) {
  
  jaccard_clust <- hclust(jaccard_dist, method = cluster_method)
  jaccard_clust <- reorder(jaccard_clust, jaccard_dist, method = "OLO")
  
  pdf(file = file.path(paste0(output_base, '-', cluster_method, '-tree.pdf') ) )
  plot(jaccard_clust)
  dev.off()
  
  # reorder matrix and plot
  jaccard_overlaps_clust <- jaccard_overlaps[ jaccard_clust$order, jaccard_clust$order ]
  # want to only plot half of the matrix
  # reverse order to match tree for plot
  genes <- rev(rownames(jaccard_overlaps_clust))

  n <- length(genes) - 1
  num_entries <- (n + 1)*n/2
  gene1 <- character(length = num_entries)
  gene2 <- character(length = num_entries)
  overlap <- numeric(length = num_entries)
  idx <- 1
  # Gene1 is for y-axis and Gene2 is for x-axis
  for (i in seq_len(length(genes))) {
    for (j in seq_len(length(genes))) {
      if (i >= j) {
        next
      } else {
        gene1[idx] <- genes[i]
        gene2[idx] <- genes[j]
        overlap[idx] <- jaccard_overlaps_clust[ genes[i], genes[j] ]
        idx <- idx + 1
      }
    }
  }
  # need to add in the first gene against itself so that all the genes appear on the x axis
  jaccard_overlaps_tophalf_m <- data.frame(
    Gene1 = factor(c(gene1, genes[1]), levels = genes),
    Gene2 = factor(c(gene2, genes[1]), levels = rev(genes)),
    Overlap = c(overlap, NA)
  )
  # plot heatmap
  clust_heatmap_plot <- matrix_heatmap(jaccard_overlaps_tophalf_m,
                                       y = 'Gene1', x = 'Gene2',
                                       fill = 'Overlap')
  
  pdf(file = file.path(paste0(output_base, '-', cluster_method, '-heatmap.pdf') ) )
  print(clust_heatmap_plot)
  dev.off()
}

# save hclust object and heatmap plot
if (!is.null(output_data_file)) {
  save(list = c('jaccard_clust', 'clust_heatmap_plot'),
       file = output_data_file)
}
