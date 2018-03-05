#!/usr/bin/env Rscript
.libPaths('./.R/lib')

library('optparse')
option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option("--debug", action="store_true", default=FALSE,
              help="Add debugging output [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default: FALSE]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'fig2.R',
    usage = "Usage: %prog [options] expt_samples_file KO_expression_file mouse_baseline zfish_baseline sig_genes" ),
  positional_arguments = 6
)

#cmd_line_args <- list(
#  options = list(directory = '/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd',
#                 verbose = FALSE ),
#  args = c('/lustre/scratch117/maz/team31/projects/mouse_DMDD/samples-minus-outliers.txt',
#           '/lustre/scratch117/maz/team31/projects/mouse_DMDD/ko_expr/ko_expr.tsv',
#           'data/Mm_GRCm38_e88_baseline.rda',
#           '/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd/data/sig_gene_counts.tsv',
#           'output/human-mim-edited.tsv',
#           'data/Dr_GRCz10_e90_baseline.rda'
#          )
#)

# make options simpler
directory <- cmd_line_args$options[['directory']]
debug <- cmd_line_args$options[['debug']]
verbose <- cmd_line_args$options[['verbose']]

if ( verbose ){
  cat( "Working directory:", directory, "\n", sep=" " )
}

if (directory == 'cwd') {
  wd <- getwd()
} else {
  wd <- directory
}

plots_dir <- file.path(wd, 'plots')

packages <- c('ggplot2', 'plyr', 'viridis', 'RColorBrewer', 'reshape2',
              'SummarizedExperiment', 'cowplot', 'grid')
for( package in packages ){
  library(package, character.only = TRUE)
}
#for( package in packages ){
#  suppressPackageStartupMessages(
#    suppressWarnings( library(package, character.only = TRUE) )
#  )
#}

if (debug) {
  cat('LOAD SAMPLE INFO...\n')
}
# load sample info
sample_file <- cmd_line_args$args[1]
sample_info <- read.table(file = sample_file, sep = "\t", header = TRUE, row.names = 1 )
sample_info$gene <- gsub('_[a-z0-9]+', '', row.names(sample_info))
# remove Cenpl
sample_info <- sample_info[ sample_info$gene != 'Cenpl', ]
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
# set levels
stage_count$gene <- factor(stage_count$gene, levels = stage_count$gene)

# get names of genes that are delayed
delayed <-
  stage_count$Severe + stage_count$Moderate + stage_count$Slight
delayed_genes <- stage_count$gene[ delayed > 0 ]
# write genes to file
write.table(delayed_genes, file = file.path(wd, 'output', 'KOs_delayed.txt'),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# output delayed order
delay <- factor(rep('None', nrow(stage_count)),
                levels = c('Severe', 'Moderate', 'Slight', 'None'))
for ( delay_categoy in c('Slight', 'Moderate', 'Severe')) {
  delay[ stage_count[[delay_categoy]] > 0 ] <- delay_categoy
}
write.table(data.frame(stage_count$gene, delay = delay),
            file = file.path(wd, 'output', 'KOs_ordered_by_delay.txt'),
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

# reshape for plotting
stage_count_m <- melt(stage_count, id.vars = c('gene'),
                      variable.name = 'delay_type',
                      value.name = 'count')
stage_count_m$gene <- factor(stage_count_m$gene,
                             levels = rev(stage_count$gene) )

# function to get legend of a ggplot object
get_gg_legend <- function(ggplot_obj){ 
  tmp <- ggplot_gtable(ggplot_build(ggplot_obj)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

# plot
embryo_stage_plot <- ggplot(data = stage_count_m) + 
  geom_raster( aes(x = delay_type, y = gene, fill = count) ) + 
  scale_fill_viridis(direction = -1) +
  scale_x_discrete(position = 'top') + 
  theme_void() + theme(axis.text.x = element_text(colour = 'black', angle = 90, hjust = 0, debug = FALSE),
                       legend.position = 'top' )

if (debug) {
  cat('EMBRYO STAGE PLOT...\n')
}

pdf(file = file.path(plots_dir, 'embryo_stage_colour.pdf'),
    width = 2, height = 5 )
print(embryo_stage_plot)
dev.off()

# make zeros appear as white
# convert zeros to NA
stage_count_na_m <- stage_count_m
stage_count_na_m$count[ stage_count_na_m$count == 0 ] <- NA

embryo_stage_zero_white_plot <- ggplot(data = stage_count_na_m) + 
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
stage_count_by_type <- split(stage_count_m, stage_count_m$delay_type)

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

# plot with stage as bar length and gt as colour
## TO DO: Add annotations at top for Delay category, somite numbers,
##        Theiler Stage and dpc
embryo_stage_size_colour_plot <- ggplot(data = stage_count.for_tiles) + 
  geom_rect( aes( xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = condition)) +
  geom_vline(data = stage_separators, aes(xintercept = raw)) +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0.5,length(stage_count$gene) - 0.5,1),
                     labels = rev(stage_count$gene) ) +
  scale_fill_manual(values = c('firebrick2', 'steelblue3', 'green'),
                    guide = 'none') +  
  theme_void() + theme( axis.text.y = element_text(size = 8, colour = 'black',
                                                   angle = 0, debug = FALSE) )

# get legend for plotting separately
embryo_stage_size_colour_plot_plus_legend <-
  embryo_stage_size_colour_plot +
    scale_fill_manual(values = c('firebrick2', 'steelblue3', 'green'),
                      guide = guide_legend(reverse = TRUE, title = "Genotype")) +
    theme( axis.text.y = element_text(size = 8, colour = 'black',
                                                   angle = 0, debug = FALSE),
      legend.position = 'top',
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7)
    )
gt_legend <- get_gg_legend(embryo_stage_size_colour_plot_plus_legend)

# plot with legend
pdf(file = file.path(plots_dir, 'embryo_stage_size_colour.pdf'),
    width = 4, height = 8 )
print(embryo_stage_size_colour_plot_plus_legend)
dev.off()
# and plot legend separately
postscript(file = file.path(plots_dir, 'embryo_stage_size_colour.legend.eps'),
           width = 6, height = 4)
grid.draw(gt_legend)
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

embryo_stage_by_gene_by_gt_plot <- ggplot(data = sample_info) +
  geom_tile(aes(x = stage_as_number, y = gene_gt, fill = condition), alpha = 0.2) +
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

# heatmap
embryo_ko_expr_plot <- ggplot(data = ko_expr) + 
  geom_tile(aes(x = gt, y = symbol, fill = log2fc )) +
  scale_x_discrete(position = 'top') +
  scale_fill_gradient2(na.value = 'grey 80',
    guide = 'none') +
  theme_void() + theme(axis.text.x =
                       element_text(size = 10, colour = 'black', angle = 90,
                                    hjust = 0, debug = FALSE))

# get legend for plotting separately
embryo_ko_expr_plot_plus_legend <-
  embryo_ko_expr_plot +
    scale_fill_gradient2(
      na.value = 'grey 80',
      guide = guide_colourbar(title = expression(paste(log[2], "[Fold Change]", sep = '') ) ) ) +
    theme( axis.text.x = element_text(size = 10, colour = 'black', angle = 90,
                                    hjust = 0, debug = FALSE),
      axis.text.y = element_text(size = 10, colour = 'black',
                                  angle = 0, debug = FALSE),
      legend.position = 'top',
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7)
    )
ko_legend <- get_gg_legend(embryo_ko_expr_plot_plus_legend)

# plot with legend
pdf(file = file.path(plots_dir, 'embryo_ko_expr_plot.pdf'),
    width = 2, height = 5)
print(embryo_ko_expr_plot_plus_legend)
dev.off()
# and plot legend separately
postscript(file = file.path(plots_dir, 'embryo_ko_expr_plot.legend.eps'),
           width = 6, height = 4)
grid.draw(ko_legend) 
dev.off()

if (debug) {
  cat('BASELINE EXPRESSION PLOT...\n')
}

## BASELINE
# heatmap for expression of the knocked out genes in baseline
# load baseline data
load(cmd_line_args$args[3])

# get gene_ids in the same order as the heatmap
id_for <- as.character(unique(ko_expr$gene_id))
names(id_for) <- as.character(unique(ko_expr$symbol))
id_for <- id_for[ stage_count$gene ]

# get counts for those genes
baseline_counts <- assay(Mm_GRCm38_e88_baseline)[id_for, ]
names(dimnames(baseline_counts)) <- c('Gene', 'Stage')
# get means for each gene for each stage
baseline_counts_by_stage <-
  split.data.frame(t(baseline_counts),
                    colData(Mm_GRCm38_e88_baseline)$condition )
mean_counts_by_stage <- lapply(baseline_counts_by_stage, colMeans)
mean_counts <- do.call(cbind, mean_counts_by_stage)
names(dimnames(mean_counts)) <- c('Gene', 'Stage')
# mean center and scale, melt
mean_counts_scaled <- t(scale(t(mean_counts)))
mean_counts_scaled.m <- melt(mean_counts_scaled)

# heatmap
# diverging colour palette
max_fill <- max(abs(mean_counts_scaled.m$value))
mouse_baseline_heatmap <-
  ggplot(data = mean_counts_scaled.m) +
    geom_tile( aes(x = Stage, y = Gene, fill = value)) +
    scale_fill_distiller(limits = c(-max_fill, max_fill), type= 'div', palette = "RdBu") +
    scale_x_discrete(position = 'top') + 
    theme_void() +
    theme( axis.text = element_text(colour = 'black', angle = 90,
                                    hjust = 0, vjust = 1),
          axis.text.y = element_blank(),
          legend.position = 'top')

pdf(file = file.path(plots_dir, 'mouse_baseline_heatmap.pdf'))
print(mouse_baseline_heatmap)
dev.off()

# do average per Theiler stage
baseline_counts_by_theiler_stage <-
  split.data.frame(t(baseline_counts),
                    colData(Mm_GRCm38_e88_baseline)$Theiler_stage )
mean_counts_by_theiler_stage <- lapply(baseline_counts_by_theiler_stage, colMeans)
mean_counts_ts <- do.call(cbind, mean_counts_by_theiler_stage)
names(dimnames(mean_counts_ts)) <- c('Gene', 'Stage')
# mean center and scale, melt
mean_counts_ts_scaled <- t(scale(t(mean_counts_ts)))
names(dimnames(mean_counts_ts_scaled)) <- c('Gene', 'Stage')
mean_counts_ts_scaled.m <- melt(mean_counts_ts_scaled)

max_fill <- max(abs(mean_counts_ts_scaled.m$value))
mouse_baseline_ts_heatmap <-
  ggplot(data = mean_counts_ts_scaled.m) +
    geom_tile( aes(x = Stage, y = Gene, fill = value)) +
    scale_fill_distiller(limits = c(-max_fill, max_fill), type= 'div',
      palette = "RdBu", guide = 'none') +
    scale_x_discrete(position = 'top') + 
    theme_void() +
    theme( axis.text = element_text(colour = 'black', size = 10, angle = 90,
                                    hjust = 0, vjust = 1),
          axis.text.y = element_blank(),
 )

# get legend for plotting separately
mouse_baseline_ts_heatmap_plus_legend <-
  mouse_baseline_ts_heatmap +
    scale_fill_distiller(
      limits = c(-max_fill, max_fill), type= 'div',
      palette = "RdBu",
      guide =
        guide_colourbar(title = "Normalised Counts\n(Mean Centred and Scaled)")
    ) +
    theme( axis.text.x = element_text(colour = 'black', size = 10, angle = 90),
      axis.text.y = element_text(colour = 'black', size = 10,
                                 angle = 0, debug = FALSE),
      legend.position = 'top',
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7)
    )

mouse_baseline_legend <- get_gg_legend(mouse_baseline_ts_heatmap_plus_legend)

# plot with legend
pdf(file = file.path(plots_dir, 'mouse_baseline_theiler_stage_heatmap.pdf'),
    width = 2, height = 8)
print(mouse_baseline_ts_heatmap_plus_legend)
dev.off()

# and plot legend separately
postscript(file = file.path(plots_dir, 'mouse_baseline_theiler_stage_heatmap.legend.eps'),
           width = 6, height = 4)
grid.draw(mouse_baseline_legend) 
dev.off()


# log10 heatmap
mean_counts_ts.m <- melt(mean_counts_ts)
mean_counts_ts.m$log10 <- log10(mean_counts_ts.m$value + 1)
mouse_baseline_ts_log10_heatmap <-
  ggplot(data = mean_counts_ts.m) +
    geom_tile( aes(x = Stage, y = Gene, fill = log10)) +
    scale_fill_viridis() +
    scale_x_discrete(position = 'top') + 
    theme_void() +
    theme( axis.text = element_text(colour = 'black', angle = 90,
                                    hjust = 0, vjust = 1),
          axis.text.y = element_blank(),
          legend.position = 'top')

pdf(file = file.path(plots_dir, 'mouse_baseline_theiler_stage_log10_heatmap.pdf'),
    width = 2, height = 8)
print(mouse_baseline_ts_log10_heatmap)
dev.off()

if (debug) {
  cat('SIG GENES PLOT...\n')
}

# numbers of significant genes
sig_genes_file <- cmd_line_args$args[4]
sig_genes <- read.delim(sig_genes_file)
# subset to ko_response
sig_genes <- sig_genes[ sig_genes$Set == 'ko_response', ]
# order genes
sig_genes$Gene <- factor(sig_genes$Gene,
                         levels = rev(stage_count$gene))

# create new column with combination of type and comparison
sig_genes$Category <-
  factor( paste(sig_genes$Comparison, sig_genes$Type, sep = '-' ),
          levels = c('hom_vs_het_wt-unfiltered', 'hom_vs_het_wt-filtered',
                      'het_vs_wt-unfiltered', 'het_vs_wt-filtered'))

sig_genes_heatmap <- ggplot(data = sig_genes) +
  geom_tile( aes(x = Category, y = Gene, fill = Count) ) +
  scale_fill_viridis(direction = -1, na.value = 'grey 90') +
  scale_x_discrete(position = 'top') + 
  theme_void() +
  theme(axis.text.x = element_text(size = 10, colour = 'black', angle = 90,
                                    hjust = 0, vjust = 0.5, debug = FALSE),
        legend.position = 'top',
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))

# plot as log10 count as well
sig_genes$log10_count <- log10(sig_genes$Count + 1)
sig_genes_log10_heatmap <- ggplot(data = sig_genes) +
  geom_tile( aes(x = Category, y = Gene, fill = log10_count) ) +
  scale_fill_viridis(direction = -1, na.value = 'grey 90', guide = 'none') +
  scale_x_discrete(position = 'top',
                   labels = c('hom vs sibs', "post filter",
                              'het vs wt', "post filter") ) + 
  theme_void() +
  theme(axis.text.x = element_text(size = 10, colour = 'black', angle = 90,
                                    hjust = 0, debug = FALSE))

# get legend for plotting separately
sig_genes_log10_heatmap_plus_legend <-
  sig_genes_log10_heatmap +
    scale_fill_viridis(direction = -1, na.value = 'grey 90',
    guide = guide_colourbar(title =
              expression(paste(log[10], "[Sig genes count]", sep = '') ) ) ) +
    theme( axis.text.x = element_text(size = 10, colour = 'black', angle = 90,
                                    hjust = 0, debug = FALSE),
      axis.text.y = element_text(size = 10, colour = 'black',
                                 angle = 0, debug = FALSE),
      legend.position = 'top',
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7)
    )
sig_genes_log10_legend <- get_gg_legend(sig_genes_log10_heatmap_plus_legend)

# plot with legend
pdf(file = file.path(plots_dir, 'sig_genes_heatmap.pdf'),
    width = 4, height = 8)
print(sig_genes_log10_heatmap_plus_legend)
print(sig_genes_heatmap)
print(sig_genes_heatmap + scale_fill_viridis(direction = 1, na.value = 'grey 85'))
dev.off()

# and plot legend separately
postscript(file = file.path(plots_dir, 'sig_genes_heatmap.legend.eps'),
           width = 6, height = 4)
grid.draw(sig_genes_log10_legend) 
dev.off()

# make composite figure
save_plot(file.path(plots_dir, "Figure2.eps"),
          plot_grid(embryo_stage_size_colour_plot, embryo_ko_expr_plot,
          sig_genes_log10_heatmap, mouse_baseline_ts_heatmap, 
          ncol = 4, rel_widths = c(6,1,2,3), align = 'h' ),
          ncol = 4, device = 'eps',
          base_height = 9, base_aspect_ratio = 0.2 )

# save plot objects
save.image(file = file.path(wd, 'output', 'fig2.RData'))
#object_to_save = c('embryo_stage_size_colour_plot', 'embryo_ko_expr_plot',
#                   'mouse_baseline_ts_heatmap', 'sig_genes_heatmap')
#save(list = object_to_save, file = file.path(wd, 'output', 'fig1.RData'))
