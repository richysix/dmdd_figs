#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'fig1.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 4
)

#cmd_line_args <- list(
#  options = list(directory = 'cwd',
#                 verbose = FALSE ),
#  args = c('data/PRJEB4513-E8.25/4567_somites-counts.tsv',
#           'data/PRJEB4513-E8.25/4567_somites-samples.tsv',
#           'data/PRJEB4513-E8.25/tissues-counts.tsv',
#           'data/PRJEB4513-E8.25/tissues-samples.tsv' )
#)

if (cmd_line_args$options[['directory']] == 'cwd') {
  wd <- getwd()
} else {
  wd <- cmd_line_args$options[['directory']]
}
plots_dir <- file.path(wd, 'plots')

if( cmd_line_args$options[['verbose']] ){
  cat( "Working directory:", cmd_line_args$options[['directory']], "\n", sep=" " )
}

packages <- c('ggplot2', 'viridis', 'reshape2', 'SummarizedExperiment',
              'DESeq2')
for( package in packages ){
    library(package, character.only = TRUE)
}

################################################################################
## baseline data
# 4, 5, 6 and 7 somites data
we_count_file <- cmd_line_args$args[1]
we_count_data <- read.delim(we_count_file, row.names = 1, check.names = FALSE)

we_samples_file <- cmd_line_args$args[2]
we_sample_info <- read.delim(we_samples_file, row.names = 1)

################################################################################
# tissue data
tissue_count_file <- cmd_line_args$args[3]
tissue_data <- read.delim(tissue_count_file, row.names = 1)
tissue_count_data <- tissue_data[, grepl('count$', colnames(tissue_data))]
colnames(tissue_count_data) <- gsub('.count$', '', colnames(tissue_count_data))

tissue_annotation_data <- cbind('Gene ID' = rownames(tissue_data),
                                tissue_data[, !grepl('count$', colnames(tissue_data))] )

tissue_sample_file <- cmd_line_args$args[4]
tissue_sample_info <- read.delim(tissue_sample_file, row.names = 1)
tissue_sample_info$condition <- factor(tissue_sample_info$condition,
                                       levels = tissue_sample_info$condition)

################################################################################
# merge counts together
merged_counts <- merge(we_count_data, tissue_count_data, by = "row.names")
rownames(merged_counts) <- merged_counts[, 'Row.names']
merged_counts <- merged_counts[ , !grepl('Row.names', colnames(merged_counts)) ]

# make merged samples df
#merged_samples <- rbind(sample_info[, c('condition', 'type')], tissue_sample_info)
merged_samples <- rbind(we_sample_info, tissue_sample_info)

# make DESeq2 object and normalise
dds <- DESeqDataSetFromMatrix(merged_counts, merged_samples, design = ~ 1)
dds <- estimateSizeFactors(dds)
mcols(dds) <- tissue_annotation_data

norm_counts <- counts(dds, normalized = TRUE)

# count (protein-coding) genes expressed at >= 10 counts
norm_counts_protein <- norm_counts[ mcols(dds)$Biotype == 'protein_coding', ]

we_gt10 <- as.data.frame(norm_counts_protein[, merged_samples$type == 'Whole Embryo'] >= 10)
we_gt10_genes <- rownames(norm_counts_protein)[ Reduce('|', we_gt10) ]

tissue_gt10 <- as.data.frame(norm_counts_protein[, merged_samples$type == 'Tissue'] >= 10)
tissue_gt10_genes <- rownames(norm_counts_protein)[ Reduce('|', tissue_gt10) ]

venn_numbers <- data.frame(
  total_genes = length(union(we_gt10_genes, tissue_gt10_genes)),
  we = length(we_gt10_genes),
  we_only = length(setdiff(we_gt10_genes, tissue_gt10_genes)),
  intersection = length(intersect(we_gt10_genes, tissue_gt10_genes)),
  tissue_only = length(setdiff(tissue_gt10_genes, we_gt10_genes)),
  tissue = length(tissue_gt10_genes)
)
print(venn_numbers)

write.table(venn_numbers,
            file = file.path('output', 'tissue_vs_WE_Venn_numbers.tsv'),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# plot heatmap of tissue only genes
tissue_only_genes <- setdiff(tissue_gt10_genes, we_gt10_genes)
tissue_only_counts <- norm_counts[tissue_only_genes, ]

# sort by row sums
row_sums <- rowSums(tissue_only_counts)
tissue_only_counts <- tissue_only_counts[names(sort(row_sums)), ]

# calculate mean by stage for whole embryo samples for heatmap
# split matrix by stage
counts_by_stage <- split.data.frame(t(tissue_only_counts), merged_samples$condition)
# sum columns
mean_counts_by_stage_list <- lapply(counts_by_stage, colMeans)
mean_counts_by_stage <- do.call(cbind, mean_counts_by_stage_list)

# melt
tissue_only_counts_m <- melt(mean_counts_by_stage)
colnames(tissue_only_counts_m) <- c('Gene', 'Sample', 'Norm_counts')
# reverse levels of Gene for plot
tissue_only_counts_m$Gene <-
  factor(tissue_only_counts_m$Gene,
         levels = rev(levels(tissue_only_counts_m$Gene)))

max_fill <- max(abs(tissue_only_counts_m$Norm_counts))

# cut into different groups
max_fill_thou <- ceiling(max_fill/1000) * 1000
tissue_only_counts_m$Expr_group <-
  cut(tissue_only_counts_m$Norm_counts,
      breaks = c(0,1,10,20,50,max_fill_thou),
      labels = c('Zero Counts', '1-9', '10-19', '20-49', '>= 50'),
      right = FALSE)

tissue_vs_WE_heatmap <-
  ggplot(data = tissue_only_counts_m) +
    geom_tile(aes(x = Sample, y = Gene, fill = Expr_group)) +
    scale_fill_viridis(discrete=TRUE, name = 'Normalised\nCounts') +
    scale_x_discrete(position = 'top') + 
    theme_void() +
    theme( axis.text = element_text(colour = 'black', angle = 45,
                                    hjust = 0, vjust = 1),
          axis.text.y = element_blank(),
          legend.position = 'right')

pdf(file = file.path(plots_dir, 'tissue_vs_WE_heatmap.pdf'),
    width = 6, height = 9, paper = "special")
print(tissue_vs_WE_heatmap)
dev.off()

# save an eps as well
postscript(file = file.path(plots_dir, 'tissue_vs_WE_heatmap.eps'),
    width = 6, height = 9, paper = "special", horizontal = FALSE)
print(tissue_vs_WE_heatmap)
dev.off()

# count up genes with >= 50 counts in tissues
highly_expressed <-
  as.data.frame(
    table(
      tissue_only_counts_m$Sample[ tissue_only_counts_m$Expr_group == '>= 50' ]
    )
  )
names(highly_expressed) <- c('Sample', '>= 50')

# find out how many genes have at least one tissue with >= 50 counts
tissue_only_gt50 <- as.data.frame(tissue_only_counts[, merged_samples$type == 'Tissue'] >= 50)
tissue_only_gt50_genes <- rownames(tissue_only_counts)[ Reduce('|', tissue_only_gt50) ]

all_tissues_gt50_genes <- rownames(tissue_only_counts)[ Reduce('&', tissue_only_gt50) ]

highly_expressed <- rbind(highly_expressed,
                          data.frame( 'Sample' = c('Any tissue', 'All tissues'),
                                     '>= 50' = c(length(tissue_only_gt50_genes), length(all_tissues_gt50_genes)),
                                     check.names = FALSE
                                    ))
write.table(highly_expressed,
            file = file.path('output', 'highly_expressed-tissues.tsv'),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# output normalised counts of tissue only genes
colnames(mean_counts_by_stage) <- gsub('$', ' normalised count', colnames(mean_counts_by_stage))

tissue_only_counts_with_anno <- cbind(as.data.frame(rowData(dds[rownames(mean_counts_by_stage),])),
                                      mean_counts_by_stage)

write.table(tissue_only_counts_with_anno,
            file = file.path('output', 'tissue_only-counts.tsv'),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

################################################################################
# save plot objects
save.image(file = file.path(wd, 'output', 'fig1.RData'))

################################################################################
