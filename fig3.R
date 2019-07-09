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
    usage = "Usage: %prog [options] gene_info_file expt_samples_file Genes_order_file KO_expression_file mouse_baseline sig_genes OMIM_data" ),
  positional_arguments = 7
)

#cmd_line_args <- list(
#  options = list(directory = 'cwd',
#                 verbose = FALSE ),
#  args = c('/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/dmdd-genes.txt',
#           'output/sample_info.txt',
#           'output/KOs_ordered_by_delay.txt',
#           '/lustre/scratch117/maz/team31/projects/mouse_DMDD/ko_expr/ko_expr.tsv',
#           'data/Mm_GRCm38_e88_baseline.rda',
#           'data/sig_gene_counts.tsv',
#           'output/human-mim-edited.tsv'
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
# load gene info
gene_file <- cmd_line_args$args[1]
gene_info <- read.delim(gene_file)
# construct named vector of gene names for directory names
gene_name_for_dir <- as.character(gene_info$Gene.Name)
names(gene_name_for_dir) <- gene_info$Dir

# load sample info
sample_file <- cmd_line_args$args[2]
sample_info <- read.table(file = sample_file, sep = "\t", header = TRUE, row.names = 1 )

# set levels of gt
sample_info$condition <- factor(sample_info$condition,
                                levels = c('hom', 'het', 'wt'))

# calculate percentage of hom embryos that are delayed (< 20 somites)
pc_delayed_homs <- sum( sample_info$somite_number[ sample_info$condition == 'hom'] < 20, na.rm = TRUE) / nrow(sample_info[ sample_info$condition == 'hom', ]) * 100
#print(sprintf('%.1f%% of homs are delayed\n', pc_delayed_homs))
print(sprintf('num = %d total = %d %.1f%% of homs are delayed',
sum( sample_info$somite_number[ sample_info$condition == 'hom'] < 20, na.rm = TRUE),
nrow(sample_info[ sample_info$condition == 'hom', ]), pc_delayed_homs))

# calculate mean number of embryos per condition per expt
sample_info_by_expt <- split(sample_info, sample_info$gene)
median_homs_by_expt <- median(sapply(sample_info_by_expt,
                            function(df){
                              num_homs <- sum(df$condition == 'hom')
                              if (num_homs == 0) {
                                return(NA)
                              } else {
                                return(num_homs)
                              }
                            }), na.rm = TRUE)

median_sibs_by_expt <- median(sapply(sample_info_by_expt,
                            function(df){
                              num_homs <- sum(df$condition == 'hom')
                              num_sibs <- sum(df$condition != 'hom')
                              if (num_homs == 0) {
                                return(NA)
                              } else {
                                return(num_sibs)
                              }
                            }), na.rm = TRUE)
print(sprintf('The median number of homs and sibs per expt are %d and %d respectively\n', median_homs_by_expt, median_sibs_by_expt))

# set levels of Theiler stage
sample_info$Theiler_stage <-
  factor(sample_info$Theiler_stage,
          levels = c('TS12b', 'TS13', 'TS14', 'TS15', 'Other'))

if (debug) {
  cat('EMBRYO STAGE PLOT...\n')
}

# plot number of embryos as size of box
# split by delay_type
stage_count_by_gt <-
  ddply(sample_info, .(gene, condition), .drop = FALSE, summarise,
    Severe = sum(Theiler_stage == 'TS12b'),
    Moderate = sum(Theiler_stage == 'TS13'),
    Slight = sum(Theiler_stage == 'TS14'),
    None = sum(Theiler_stage == 'TS15')
  )
  
# read in KO ordering
ko_order_file <- cmd_line_args$args[3]
ko_order <- read.table(ko_order_file)
names(ko_order) <- c('Gene', 'Delay Category')

# order genes in the same order as stage_count
stage_count_by_gt <- do.call(rbind, lapply(as.character(ko_order$Gene),
       function(x){ stage_count_by_gt[ stage_count_by_gt$gene == x, ] } )
)

# melt and split by delay_type
stage_count_by_gt.m <- melt(stage_count_by_gt, id.vars = c('gene', 'condition'),
                      variable.name = 'delay_type',
                      value.name = 'count')
stage_count_by_gt.m$gene <- factor(stage_count_by_gt.m$gene,
                                   levels = rev(as.character(ko_order$Gene)) )

stage_count_by_gt.m_by_type <- split(stage_count_by_gt.m, stage_count_by_gt.m$delay_type)
# also split stage_count by delay_type
stage_count.m <- ddply(stage_count_by_gt.m, .(gene, delay_type), summarise, count = sum(count))
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
  offset <- offset + max(stage_count_by_type[[i]]$count) #max( sapply(split(stage_count_by_gt[['Severe']], stage_count_by_gt$gene), sum) )
  stage_separators$raw[i+1] <- offset
}
stage_count.for_tiles <- do.call(rbind, stage_count.for_tiles_list)

# plot with stage as bar length and gt as colour
# convert gene (directory) name to correct e88 gene name
gene_names <- gene_name_for_dir[ rev(as.character(ko_order$Gene)) ]

embryo_stage_size_colour_plot <- ggplot(data = stage_count.for_tiles) + 
  geom_rect( aes( xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = condition)) +
  geom_vline(data = stage_separators, aes(xintercept = raw)) +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0.5,length(stage_count$gene) - 0.5,1),
                     labels =  gene_names) +
  scale_fill_manual(values = c('firebrick2', 'steelblue3', 'green'),
                    guide = 'none') +  
  theme_void() + theme( axis.text.y = element_text(size = 8, colour = 'black',
                                                   face = 'italic', angle = 0, debug = FALSE) )

# function to get legend of a ggplot object
get_gg_legend <- function(ggplot_obj){ 
  tmp <- ggplot_gtable(ggplot_build(ggplot_obj)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

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
invisible(dev.off())
# and plot legend separately
postscript(file = file.path(plots_dir, 'embryo_stage_size_colour.legend.eps'),
           width = 6, height = 4)
grid.draw(gt_legend)
invisible(dev.off())

if (debug) {
  cat('KO EXPRESSION PLOT...\n')
}

# expression of the knocked out gene in homs and hets
ko_expr_file <- cmd_line_args$args[4]
ko_expr <- read.table(ko_expr_file, sep = "\t", header = FALSE )
names(ko_expr) <- c('gene_id', 'symbol', 'comparison', 'log2fc')
ko_expr$gt <- factor( gsub('_vs_.*', '', ko_expr$comparison),
                     levels = c('hom', 'het', 'wt') )
ko_expr$symbol <- factor(ko_expr$symbol,
                          levels = rev(stage_count$gene))

# add a df for wt (log2fc = 0)
num_genes <- nlevels(ko_expr$gene_id)

ko_expr <-
  rbind(ko_expr,
    data.frame(
      gene_id = unique(ko_expr$gene_id),
      symbol = unique(ko_expr$symbol),
      comparison = rep('wt', num_genes),
      log2fc = rep(0, num_genes),
      gt = rep('wt', num_genes)
    )
  )

# get count data for ko expression to add to heatmap
get_mean_counts <- function(gene_id, sample_info){
  # open counts file
  #counts_file <- file.path(wd, 'data', paste0(gene_id, '.expr.tsv'))
  counts_file <-
    file.path('/lustre/scratch117/maz/team31/projects/mouse_DMDD/ko_expr',
                paste0(gene_id, '.expr.tsv'))
  counts_data <- read.delim(counts_file, check.names = FALSE)
  # get counts columns
  counts <- counts_data[ , grepl('normalised.count', names(counts_data)) ]
  names(counts) <- gsub('.normalised.count', '', names(counts))
  counts <- melt(counts, variable.name = 'sample',
                 value.name = 'normalised count', id.vars = c())
  
  counts <- merge( sample_info, counts,
                    by.x = "row.names", by.y = c('sample') )
  # calculate mean for each genotype
  mean_count_by_gt <-
    melt( sapply( split(counts, counts$condition),
            function(counts_df){ mean(counts_df[['normalised count']]) } ),
         value.name = 'Mean.Count', id.vars = c())
  mean_count_by_gt$gt <- factor( row.names(mean_count_by_gt),
                                levels = c('hom', 'het', 'wt') )
  mean_count_by_gt$gene_id <- rep(gene_id, nrow(mean_count_by_gt))
  # add gene name in from ko_expr data frame
  mean_count_by_gt <- merge(mean_count_by_gt,
                            unique(ko_expr[ , c('gene_id', 'symbol')]) )
  
  return(mean_count_by_gt)
}

mean_count_by_gt <- do.call(rbind,
                            lapply(as.character(unique(ko_expr$gene_id)),
                                    get_mean_counts,
                                    sample_info) )
# round to 0dp
mean_count_by_gt$Mean.Count <- round(mean_count_by_gt$Mean.Count)

# subset to wt only
wt_mean_count <- mean_count_by_gt[ mean_count_by_gt$gt == 'wt', ]
# add in means for Ift140 and Oaz1
mean_count <- rbind(wt_mean_count,
                    mean_count_by_gt[ mean_count_by_gt$symbol == 'Ift140' &
                                     mean_count_by_gt$gt == 'het', ],
                    mean_count_by_gt[ mean_count_by_gt$symbol == 'Oaz1' &
                                     mean_count_by_gt$gt == 'het', ]
                   )

write.table(mean_count, file = file.path('output', 'mean_expr_wt.tsv'),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# heatmap
embryo_ko_expr_plot <- ggplot(data = ko_expr) + 
  geom_tile(aes(x = gt, y = symbol, fill = log2fc )) +
  geom_text(data = mean_count,
            aes(x = gt, y = symbol, label = Mean.Count),
            size = 2, hjust = 0.5 ) +
  scale_x_discrete(position = 'top') +
  scale_fill_gradient2(na.value = 'grey 80', guide = 'none') +
  theme_void() +
  theme(axis.text.x = element_text(size = 10, colour = 'black', angle = 90,
                                    hjust = 0, debug = FALSE) )

# get legend for plotting separately
embryo_ko_expr_plot_plus_legend <-
  embryo_ko_expr_plot +
    scale_y_discrete(labels = gene_names) +
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
    width = 4, height = 8)
print(embryo_ko_expr_plot_plus_legend)
invisible(dev.off())
# and plot legend separately
postscript(file = file.path(plots_dir, 'embryo_ko_expr_plot.legend.eps'),
           width = 6, height = 4)
grid.draw(ko_legend) 
invisible(dev.off())

if (debug) {
  cat('BASELINE EXPRESSION PLOT...\n')
}

## BASELINE
# heatmap for expression of the knocked out genes in baseline
# load baseline data
load(cmd_line_args$args[5])

# get gene_ids in the same order as the heatmap
id_for <- as.character(unique(ko_expr$gene_id))
names(id_for) <- as.character(unique(ko_expr$symbol))
gene_ids <- id_for[ as.character(stage_count$gene) ]

# get counts for those genes
baseline_subset <- Mm_GRCm38_e88_baseline[gene_ids, ]
baseline_counts <- assays(baseline_subset)$norm_counts
# get means for each gene for each stage
baseline_counts_by_stage <-
  split.data.frame(t(baseline_counts),
                    colData(baseline_subset)$condition )
mean_counts_by_stage <- lapply(baseline_counts_by_stage, colMeans)
mean_counts <- do.call(cbind, mean_counts_by_stage)
# mean center and scale, melt
mean_counts_scaled <- t(scale(t(mean_counts)))
mean_counts_scaled_df <- data.frame(
  gene_name = factor(rowData(baseline_subset)$Name,
                     levels = rev(as.character(rowData(baseline_subset)$Name))),
  mean_counts_scaled,
  check.names = FALSE
)
mean_counts_scaled_m <- melt(mean_counts_scaled_df, id.vars = c('gene_name'),
                             variable.name = 'Stage')

# heatmap
# diverging colour palette
max_fill <- max(abs(mean_counts_scaled_m$value))
mouse_baseline_heatmap <-
  ggplot(data = mean_counts_scaled_m) +
    geom_tile( aes(x = Stage, y = gene_name, fill = value)) +
    scale_fill_distiller(limits = c(-max_fill, max_fill), type= 'div', palette = "RdBu") +
    scale_x_discrete(position = 'top') + 
    theme_void() +
    theme( axis.text.x = element_text(colour = 'black', size = 10, angle = 90,
                                      debug = FALSE, hjust = 0),
          axis.text.y = element_text(colour = 'black', size = 10, angle = 0,
                                      debug = FALSE, vjust = 0),
          legend.position = 'top')

pdf(file = file.path(plots_dir, 'mouse_baseline_heatmap.pdf'))
print(mouse_baseline_heatmap)
invisible(dev.off())

# do average per Theiler stage
baseline_counts_by_theiler_stage <-
  split.data.frame(t(baseline_counts),
                    colData(baseline_subset)$Theiler_stage )
mean_counts_by_theiler_stage <- lapply(baseline_counts_by_theiler_stage, colMeans)
mean_counts_ts <- do.call(cbind, mean_counts_by_theiler_stage)
# mean center and scale, melt
mean_counts_ts_scaled <- t(scale(t(mean_counts_ts)))
mean_counts_ts_scaled_df <- data.frame(
  gene_name = factor(rowData(baseline_subset)$Name,
                     levels = rev(as.character(rowData(baseline_subset)$Name))),
  mean_counts_ts_scaled,
  check.names = FALSE
)

write.table(mean_counts_ts_scaled_df, file = file.path('output', 'mean_expr_by_TS_baseline.tsv'),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

mean_counts_ts_scaled_m <- melt(mean_counts_ts_scaled_df, id.vars = c('gene_name'),
                             variable.name = 'Stage')

max_fill <- max(abs(mean_counts_ts_scaled_m$value))
mouse_baseline_ts_heatmap <-
  ggplot(data = mean_counts_ts_scaled_m) +
    geom_tile( aes(x = Stage, y = gene_name, fill = value)) +
    scale_fill_distiller(limits = c(-max_fill, max_fill), type= 'div',
      palette = "RdBu", guide = 'none') +
    scale_x_discrete(position = 'top') + 
    theme_void() +
    theme( axis.text = element_text(colour = 'black', size = 10, angle = 90),
          axis.text.y = element_blank(),
 )

# get legend for plotting separately
mouse_baseline_ts_heatmap_plus_legend <-
  mouse_baseline_ts_heatmap +
    scale_y_discrete(labels = gene_names) +
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
    width = 5, height = 8)
print(mouse_baseline_ts_heatmap_plus_legend)
invisible(dev.off())

# and plot legend separately
postscript(file = file.path(plots_dir, 'mouse_baseline_theiler_stage_heatmap.legend.eps'),
           width = 6, height = 4)
grid.draw(mouse_baseline_legend) 
invisible(dev.off())

# log10 heatmap
mean_counts_ts_df <- data.frame(
  gene_name = factor(rowData(baseline_subset)$Name,
                     levels = rev(as.character(rowData(baseline_subset)$Name))),
  mean_counts_ts,
  check.names = FALSE
)
mean_counts_ts_m <- melt(mean_counts_ts_df, variable.name = 'Stage')
mean_counts_ts_m$log10 <- log10(mean_counts_ts_m$value + 1)
mouse_baseline_ts_log10_heatmap <-
  ggplot(data = mean_counts_ts_m) +
    geom_tile( aes(x = Stage, y = gene_name, fill = log10)) +
    scale_fill_viridis() +
    scale_x_discrete(position = 'top') + 
    theme_void() +
    theme( axis.text = element_text(colour = 'black', angle = 90,
                                    hjust = 0, vjust = 0.5),
          axis.text.y = element_blank(),
          legend.position = 'top')

pdf(file = file.path(plots_dir, 'mouse_baseline_theiler_stage_log10_heatmap.pdf'),
    width = 2, height = 8)
print(mouse_baseline_ts_log10_heatmap)
invisible(dev.off())

if (debug) {
  cat('SIG GENES PLOT...\n')
}

# numbers of significant genes
sig_genes_file <- cmd_line_args$args[6]
sig_genes <- read.delim(sig_genes_file)
# subset to ko_response and remove hom_vs_het
sig_genes <- sig_genes[ sig_genes$Set == 'ko_response' &
                        sig_genes$Comparison != 'hom_vs_het', ]
sig_genes <- droplevels(sig_genes)

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
    scale_y_discrete(labels = gene_names) +
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
invisible(dev.off())

# and plot legend separately
postscript(file = file.path(plots_dir, 'sig_genes_heatmap.legend.eps'),
           width = 6, height = 4)
grid.draw(sig_genes_log10_legend) 
invisible(dev.off())

if (debug) {
  cat('OMIM DATA TABLE PLOT...\n')
}

# plot MIM data as graphical table
mim_data_file <- cmd_line_args$args[7]
mim_data <- read.delim(mim_data_file)
mim_data$Dir <- factor(mim_data$Dir, levels = rev(stage_count$gene))

mim_reformatted <- do.call(rbind,
                    lapply(split(mim_data, mim_data$Dir),
                      function(x){
                        data.frame(
                          x = 1,
                          gene = x$Dir[1],
                          mim = paste(x$mim_short[ order(x$mim_short) ],
                                      collapse = '; ')
                        )
                      }
                    )
                  )

mim_table_plot <- ggplot( data = mim_reformatted ) + 
  geom_tile( aes(x = x, y = gene), fill = "white", colour = "white" ) +
  geom_text( aes(x = x, y = gene, label = mim), size = 2.1 ) + 
  theme_void()

pdf(file = file.path(plots_dir, 'MIM-table.pdf'))
print(mim_table_plot +
      scale_y_discrete(labels = gene_names) +
      theme(axis.text.y = element_text(colour = 'black', size = 10,
                                       angle = 0, debug = FALSE) )
     )
invisible(dev.off())

if (debug) {
  cat('COMPOSITE FIGURE...\n')
}

# make composite figure
save_plot(file.path(plots_dir, "Sample_overview.eps"),
          plot_grid(embryo_stage_size_colour_plot, embryo_ko_expr_plot,
          sig_genes_log10_heatmap, mouse_baseline_ts_heatmap, mim_table_plot,
          ncol = 5, rel_widths = c(12,3,4,6,4), align = 'h' ),
          ncol = 5, device = 'eps',
          base_height = 9, base_width = 1.54 )

# save plot objects
save.image(file = file.path(wd, 'output', 'fig3.RData'))
#object_to_save = c('embryo_stage_size_colour_plot', 'embryo_ko_expr_plot',
#                   'mouse_baseline_ts_heatmap', 'sig_genes_heatmap')
#save(list = object_to_save, file = file.path(wd, 'output', 'fig1.RData'))

# remove Rplots.pdf file
unlink('Rplots.pdf')
