#!/usr/bin/env Rscript
.libPaths('/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd/.R/lib')

library('optparse')
option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default: FALSE]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'fig1.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 1
)

#cmd_line_args <- list(
#  options = list(directory = '/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd',
#                 verbose = FALSE ),
#  args = c('/lustre/scratch117/maz/team31/projects/mouse_DMDD/samples-minus-outliers.txt',
#           '/lustre/scratch117/maz/team31/projects/mouse_DMDD/KO_expr.tsv')
#)

if ( cmd_line_args$options[['verbose']] ){
  cat( "Working directory:", cmd_line_args$options[['directory']], "\n", sep=" " )
}

if (cmd_line_args$options[['directory']] == 'cwd') {
  dir <- getwd()
} else {
  dir <- cmd_line_args$options[['directory']]
}

# set up directories
for ( new_dir in c('plots', 'output') ) {
  dir_path <- file.path(dir, new_dir)
  if( !dir.exists(dir_path) ){
    dir.create(dir_path)
  }
}
plots_dir <- file.path(dir, 'plots')

packages <- c('ggplot2', 'plyr', 'viridis', 'reshape2')
#for( package in packages ){
#  library(package, character.only = TRUE)
#}
for( package in packages ){
  suppressPackageStartupMessages(
    suppressWarnings( library(package, character.only = TRUE) )
  )
}

# load sample info
sample_file <- cmd_line_args$args[1]
sample_info <- read.table(file = sample_file, sep = "\t", header = TRUE, row.names = 1 )
sample_info$gene <- gsub('_[a-z0-9]+', '', row.names(sample_info))

# label with Theiler stage
stage_boundaries <- c(4, 7, 12, 20, 29, 34)
stage_labels <- c('TS12b', 'TS13', 'TS14', 'TS15', 'Other')

# assign stage numbers with their TS
# this issues a warning about introducing NAs by coercion. suppress it.
stage_numbers <-
  suppressWarnings(as.integer( gsub('somites', '', sample_info$stage) ))

sample_info$Theiler_stage <- cut(stage_numbers, breaks = stage_boundaries, labels = stage_labels)
# 4 embryos are of Unknown stage. Put these in Other with the one TS16 (30s) embryo
sample_info$Theiler_stage[is.na(sample_info$Theiler_stage)] <- 'Other'

# sort KOs by delay
# count up numbers of embryos in each stage
stage_count <- ddply(sample_info, .(gene), summarise,
                     Severe = sum(Theiler_stage == 'TS12b'),
                     Moderate = sum(Theiler_stage == 'TS13'),
                     Slight = sum(Theiler_stage == 'TS14'),
                     None = sum(Theiler_stage == 'TS15')
                     )
# order by Severe, Moderate, Slight then None
stage_count <- 
  stage_count[order(stage_count$Severe, stage_count$Moderate, 
                    stage_count$Slight, stage_count$None, decreasing = TRUE), ]
# reshape for plotting
stage_count.m <- melt(stage_count, id.vars = c('gene'),
                      variable.name = 'delay_type',
                      value.name = 'count')
stage_count.m$gene <- factor(stage_count.m$gene,
                             levels = rev(stage_count$gene) )

# plot
embryo_stage_plot <- ggplot(data = stage_count.m) + 
  geom_raster( aes(x = delay_type, y = gene, fill = count) ) + 
  scale_fill_viridis(direction = -1) +
  scale_x_discrete(position = 'top') + 
  theme_void() + theme(axis.text.x = element_text(colour = 'black', angle = 90, hjust = 0, debug = FALSE),
                       legend.position = 'top' )

pdf(file = file.path(plots_dir, 'embryo_stage_colour.pdf'),
    width = 2, height = 5 )
print(embryo_stage_plot)
dev.off()

# make zeros appear as white
# convert zeros to NA
stage_count_na.m <- stage_count.m
stage_count_na.m$count[ stage_count_na.m$count == 0 ] <- NA

embryo_stage_zero_white_plot <- ggplot(data = stage_count_na.m) + 
  geom_raster( aes(x = delay_type, y = gene, fill = count) ) + 
  scale_fill_viridis(direction = -1, na.value = 'white') +
  scale_x_discrete(position = 'top') + 
  theme_void() + theme(axis.text.x =
                       element_text(size = 10, colour = 'black', angle = 90,
                                    hjust = 0, debug = FALSE),
                       legend.position = 'top',
                       legend.title = element_text(size = 10))

pdf(file = file.path(plots_dir, 'embryo_stage_colour_zero_white.pdf'),
    width = 2, height = 5 )
print(embryo_stage_zero_white_plot)
dev.off()

postscript(file = file.path(plots_dir, 'embryo_stage_colour_zero_white.eps'),
    width = 2, height = 5 )
print(embryo_stage_zero_white_plot)
dev.off()

# also plot number of embryos as size of box
# split by delay_type
stage_count.m_by_type <- split(stage_count.m, stage_count.m$delay_type)

# figure out max width of column and calculate widths and positions of boxes
calculate_boxes <- function(counts_df, x_offset = 0, x_offset_norm = 0, normalised_width = 10){
  n_genes <- nrow(counts_df)
  # max box width
  max_width <- max(counts_df$count)
  ymax <- n_genes:1
  ymin <- (n_genes - 1):0
  xmax <- rep(x_offset + max_width, n_genes)
  xmax_norm <- rep(x_offset_norm + normalised_width, n_genes)
  base_xmax <- x_offset + max_width
  xmin <- vector('integer', length = n_genes)
  boxes_list <- vector('list', length = n_genes)
  for (i in seq_len(n_genes)) {
    if (counts_df$count[i] == 0) {
      xmin[i] = base_xmax
    } else {
      xmin[i] = base_xmax - counts_df$count[i]
    }
  }
  xmin_norm <- xmax_norm - ( (xmax - xmin) / max_width * normalised_width )
  # remove one with no data
  to_keep <- xmin != xmax
  return(
    data.frame(
      gene = counts_df$gene[to_keep],
      xmin = xmin[to_keep],
      xmax = xmax[to_keep],
      ymin = ymin[to_keep],
      ymax = ymax[to_keep],
      xmin_norm = xmin_norm[to_keep],
      xmax_norm = xmax_norm[to_keep]
    )
  )
}

stage_count.for_tiles_list <- vector('list', length = length(stage_count.m_by_type))
offset <- 0
offset_norm <- 0
normalised_width <- 10
stage_separators <- data.frame(
  class = c('Left', 'Severe', 'Moderate', 'Slight', 'Right'),
  raw = rep(0, 5),
  normalised = rep(0, 5)
)
for ( i in seq_len(length(stage_count.m_by_type))) {
  stage_count.for_tiles_list[[i]] <-
    calculate_boxes(stage_count.m_by_type[[i]], offset, offset_norm, normalised_width)
  offset <- offset + max(stage_count.m_by_type[[i]]$count)
  offset_norm <- offset_norm + normalised_width
  stage_separators$raw[i+1] <- offset
  stage_separators$normalised[i+1] <- offset_norm
}
stage_count.for_tiles <- do.call(rbind, stage_count.for_tiles_list)

embryo_stage_size_plot <- ggplot(data = stage_count.for_tiles) + 
  geom_rect( aes( xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'steelblue3') +
  geom_vline(data = stage_separators, aes(xintercept = raw)) + 
  theme_void()

pdf(file = file.path(plots_dir, 'embryo_stage_size.pdf'),
    width = 2, height = 5 )
print(embryo_stage_size_plot)
dev.off()

# plot equalised column widths
embryo_stage_size_norm_plot <- ggplot(data = stage_count.for_tiles) + 
  geom_rect(aes( xmin = xmin_norm, xmax = xmax_norm, ymin = ymin, ymax = ymax),
            fill = 'steelblue3') +
  geom_vline(data = stage_separators, aes(xintercept = normalised)) + 
  theme_void()

pdf(file = file.path(plots_dir, 'embryo_stage_size_norm.pdf'),
    width = 2, height = 5 )
print(embryo_stage_size_norm_plot)
dev.off()

