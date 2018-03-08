#!/usr/bin/env Rscript
.libPaths('./.R/lib')

library('optparse')

option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'mean_expression_by_mut_gt.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 3
)

#cmd_line_args <- list(
#  options = list(directory = 'cwd',
#                 verbose = FALSE ),
#  args = c('output/all_samples_merged.counts.tsv',
#           'output/all_mutants-samples.tsv',
#           'output/KOs_ordered_by_delay.txt')
#)


if (cmd_line_args$options[['directory']] == 'cwd') {
  wd <- getwd()
} else {
  wd <- cmd_line_args$options[['directory']]
}

if( cmd_line_args$options[['verbose']] ){
  cat( "Working directory:", cmd_line_args$options[['directory']], "\n", sep=" " )
  cat( "Output files basename:", cmd_line_args$options[['basename']], "\n", sep=" " )
}

packages <- c('plyr')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

gene_expr_data_file <- cmd_line_args$args[1]
gene_expr_data <- read.delim(gene_expr_data_file, check.names = FALSE)

samples_file <- cmd_line_args$args[2]
samples <- read.delim(samples_file, col.names = c('sample_name', 'Genotype', 'Gene'))
# set levels of condition and reorder
samples$Genotype <- factor(samples$Genotype,
                            levels = c('wt', 'het', 'hom'))
# order samples by Gene and then Genotype
samples <- arrange(samples, Gene, Genotype)

# get just counts
count_data <- gene_expr_data[, grepl(" count$", names(gene_expr_data)) &
                              !grepl(" normalised count$", names(gene_expr_data))]
names(count_data) <- gsub(" count$", "", names(count_data))

# get just gene info
gene_annotation <- gene_expr_data[, !grepl(" count$", names(gene_expr_data))]

# Subset and reorder count data
count_data <- count_data[, samples$sample_name]

# split matrix by mut and gt
mut_gt <- factor( paste(samples$Gene, samples$Genotype, sep = '.'),
                 levels = unique(paste(samples$Gene, samples$Genotype, sep = '.')) )
counts_by_mut_gt <- split.data.frame(t(count_data), mut_gt)
# calc mean of each column
mean_counts_by_mut_gt_list <- lapply(counts_by_mut_gt, colMeans)

mean_counts_by_mut_gt <- do.call(cbind, mean_counts_by_mut_gt_list)
row.names(mean_counts_by_mut_gt) <- gene_expr_data[['Gene ID']]

# read in delayed info
delayed_info <- read.delim(cmd_line_args$args[3], header = FALSE)
names(delayed_info) <- c('Gene', 'Delay Category')
# set levels of Delay
delayed_info[['Delay Category']] <-
  factor(delayed_info[['Delay Category']],
          levels = rev(unique(delayed_info[['Delay Category']])) )

# create samples file
samples_merged <- merge(unique(samples[, c('Gene', 'Genotype')]), delayed_info)
samples_merged$condition <- paste(samples_merged$Gene, samples_merged$Genotype, sep = '.')
samples_sorted <- arrange(samples_merged, `Delay Category`)
rownames(samples_sorted) <- samples_sorted$condition
write.table(samples_sorted,
            file = file.path(wd, 'output', 'samples_by_mut_by_gt.txt'),
            quote = FALSE, sep = "\t",
            row.names = TRUE, col.names = NA)

# round counts to integers and join to annotation
mean_counts_by_mut_gt <- round(mean_counts_by_mut_gt)
# add count to column names
colnames(mean_counts_by_mut_gt) <- gsub('$', ' count', colnames(mean_counts_by_mut_gt))
mean_counts_by_mut_gt_merged <- merge(gene_annotation, mean_counts_by_mut_gt,
                               by.x = 'Gene ID', by.y = 'row.names')

write.table(mean_counts_by_mut_gt_merged,
            file = file.path(wd, 'output', 'mean_by_mut_gt.counts.tsv'),
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
