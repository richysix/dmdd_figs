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
  wd <- getwd()
} else {
  wd <- cmd_line_args$options[['directory']]
}

# set up directories
for ( new_dir in c('plots', 'output') ) {
  dir_path <- file.path(wd, new_dir)
  if( !dir.exists(dir_path) ){
    dir.create(dir_path)
  }
}
plots_dir <- file.path(wd, 'plots')

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
# set levels of gt
sample_info$condition <- factor(sample_info$condition,
                                levels = c('hom', 'het', 'wt'))

# label with Theiler stage
stage_boundaries <- c(4, 7, 12, 19, 29, 34)
stage_labels <- c('TS12b', 'TS13', 'TS14', 'TS15', 'Other')

# assign stage numbers with their TS
# this issues a warning about introducing NAs by coercion. suppress it.
sample_info$stage_as_number <-
  suppressWarnings(as.integer( gsub('somites', '', sample_info$stage) ))

sample_info$Theiler_stage <-
  cut(sample_info$stage_as_number,
      breaks = stage_boundaries, labels = stage_labels)
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
stage_count_by_gt <-
  ddply(sample_info, .(gene, condition), .drop = FALSE, summarise,
    Severe = sum(Theiler_stage == 'TS12b'),
    Moderate = sum(Theiler_stage == 'TS13'),
    Slight = sum(Theiler_stage == 'TS14'),
    None = sum(Theiler_stage == 'TS15')
  )
# order genes in the same order as stage_count
stage_count_by_gt <- do.call(rbind, lapply(stage_count$gene,
       function(x){ stage_count_by_gt[ stage_count_by_gt$gene == x, ] } )
)

# melt and split by delay_type
stage_count_by_gt.m <- melt(stage_count_by_gt, id.vars = c('gene', 'condition'),
                      variable.name = 'delay_type',
                      value.name = 'count')
stage_count_by_gt.m$gene <- factor(stage_count_by_gt.m$gene,
                                   levels = rev(stage_count$gene) )

stage_count_by_gt.m_by_type <- split(stage_count_by_gt.m, stage_count_by_gt.m$delay_type)
# also split stage_count by delay_type
stage_count_by_type <- split(stage_count.m, stage_count.m$delay_type)

# figure out max width of column and calculate widths and positions of boxes
calculate_boxes <- function(delay_type, stage_count_by_gt.m_by_type,
                            stage_count_by_type, x_offset = 0){
  counts_df <- stage_count_by_gt.m_by_type[[delay_type]]
  counts_df_agg <- stage_count_by_type[[delay_type]]
  n_genes <- nrow(counts_df_agg)
  # max box width
  max_width <- max(counts_df_agg$count)
  ymax <- rep(n_genes:1, each = 3)
  ymin <- rep((n_genes - 1):0, each = 3)
  base_xmax <- x_offset + max_width
  
  xmin <- vector('integer', length = n_genes * 3)
  xmax <- vector('integer', length = n_genes * 3)
  gts <- unique(counts_df$condition)
  for (gene_i in seq_len(n_genes)) {
    current_xmax <- base_xmax
    for (gt_i in seq_len(length(gts)) ) {
      i <- (gene_i - 1) * 3 + gt_i
      xmax[i] <- current_xmax
      xmin[i] <- current_xmax - counts_df$count[i]
      current_xmax <- xmin[i]
    }
  }
  
  # add calculated values to data frame
  counts_df$xmin <- xmin
  counts_df$xmax <- xmax
  counts_df$ymin <- ymin
  counts_df$ymax <- ymax
  
  return(counts_df)
}

stage_count.for_tiles_list <- vector('list', length = length(stage_count_by_gt.m_by_type))
offset <- 0
#offset_norm <- 0
#normalised_width <- 10
stage_separators <- data.frame(
  class = c('Left', 'Severe', 'Moderate', 'Slight', 'Right'),
  raw = rep(0, 5)
  #normalised = rep(0, 5)
)
for ( i in seq_len(length(stage_count_by_gt.m_by_type))) {
  delay_type <- names(stage_count_by_gt.m_by_type)[i]
  stage_count.for_tiles_list[[i]] <-
    calculate_boxes(delay_type, stage_count_by_gt.m_by_type,
                    stage_count_by_type, offset)
  offset <- offset + max(stage_count_by_type[[i]]$count)
  stage_separators$raw[i+1] <- offset
}
stage_count.for_tiles <- do.call(rbind, stage_count.for_tiles_list)


embryo_stage_size_colour_plot <- ggplot(data = stage_count.for_tiles) + 
  geom_rect( aes( xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = condition)) +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0.5,73.5,1),
                     labels = rev(stage_count$gene) ) +
  scale_fill_manual(values = c('firebrick2', 'green', 'steelblue3'),
                    guide = guide_legend(reverse = TRUE)) +
  geom_vline(data = stage_separators, aes(xintercept = raw)) + 
  theme_void() + theme( axis.text.y = element_text(size = 8, colour = 'black', angle = 0, debug = FALSE),
                        legend.position = 'top',
                        legend.title = element_text(size = 9),
                        legend.text = element_text(size = 7),
                        legend.key.size = unit(1, 'lines') )

pdf(file = file.path(plots_dir, 'embryo_stage_size_colour.pdf'),
    width = 4, height = 8 )
print(embryo_stage_size_colour_plot)
dev.off()

# plot x position as stage
# order levels of gene_gt by delay and condition
sample_info$gene_gt <-
  factor(
    paste(sample_info$gene, sample_info$condition, sep="-"),
    levels = rev( paste(rep(stage_count$gene, each = 3),
                        c('wt', 'het', 'hom'), sep="-") )
  )

# create stage boundaries
ts_boundaries <- data.frame(
  Stage = c(3.5, 7.5, 12.5, 19.5, 29.5),
  Label = c('TS12a', 'TS12b', 'TS13', 'TS14', 'TS15')
)

embryo_stage_by_gene_by_gt_plot <-ggplot(data = sample_info) +
  geom_tile(aes(x = stage_as_number, y = gene_gt, fill = condition)) +
  scale_fill_manual(values = c('firebrick2', 'green', 'steelblue3'),
                    guide = guide_legend(reverse = TRUE)) +
  geom_vline(data = ts_boundaries, aes(xintercept = Stage)) + 
  theme_void() + theme(legend.position = 'top',
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6),
                     legend.key.size = unit(0.7, 'lines') )

pdf(file = file.path(plots_dir, 'embryo_stage_by_gene_by_gt.pdf'),
    width = 2, height = 5 )
print(embryo_stage_by_gene_by_gt_plot)
dev.off()


# expression of the knocked out gene in homs and hets
ko_expr_file <- cmd_line_args$args[2]
ko_expr <- read.table(ko_expr_file, sep = "\t", header = FALSE )
names(ko_expr) <- c('gene_id', 'symbol', 'comparison', 'log2fc')
ko_expr$gt <- factor( gsub('_vs_.*', '', ko_expr$comparison),
                     levels = c('hom', 'het') )
ko_expr$symbol <- factor(ko_expr$symbol,
                          levels = rev(stage_count$gene))

# 
embryo_ko_expr_plot <- ggplot(data = ko_expr) + 
  geom_tile(aes(x = gt, y = gene_id, fill = log2fc )) +
  scale_x_discrete(position = 'top') +
  scale_fill_gradient2(na.value = 'white') +
  theme_void() + theme(axis.text.x =
                       element_text(size = 10, colour = 'black', angle = 90,
                                    hjust = 0, debug = FALSE),
                       legend.position = 'top',
                       legend.title = element_text(size = 10))

pdf(file = file.path(plots_dir, 'embryo_ko_expr_plot.pdf'),
    width = 2, height = 5)
print(embryo_ko_expr_plot)
dev.off()


save.image(file = file.path(wd, 'fig1.RData'))
