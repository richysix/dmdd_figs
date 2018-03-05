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
  positional_arguments = 2
)

#cmd_line_args <- list(
#  options = list(directory = '/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd',
#                 verbose = FALSE ),
#  args = c('data/Mm_GRCm38_e88_baseline.rda',
#           '/lustre/scratch117/maz/team31/projects/mouse_DMDD/PRJEB4513-E8.25/all.tsv' )
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
load(cmd_line_args$args[1])
counts <- assays(Mm_GRCm38_e88_baseline)$counts[ ,
              grepl('^[567]somites_[0-9]$', colnames(Mm_GRCm38_e88_baseline)) ]

# subset to 5, 6 and 7 somites and drop levels
sample_info <-
  colData(Mm_GRCm38_e88_baseline)[
            grepl('^[567]somites_[0-9]$', colnames(Mm_GRCm38_e88_baseline)), ]
sample_info <- droplevels(sample_info)

# split matrix by stage
counts_by_stage <- split.data.frame(t(counts), sample_info$condition)
# sum columns
sum_counts_by_stage_list <- lapply(counts_by_stage, colSums)

sum_counts_by_stage <- do.call(cbind, sum_counts_by_stage_list)
row.names(sum_counts_by_stage) <- rowData(Mm_GRCm38_e88_baseline)$Gene.ID

samples_by_stage <- data.frame(
  row.names = colnames(sum_counts_by_stage),
  condition = colnames(sum_counts_by_stage),
  Type      = rep('Whole Embryo', ncol(sum_counts_by_stage))
)


################################################################################
# tissue data
tissue_count_file <- cmd_line_args$args[2]
tissue_count_data <- read.delim(tissue_count_file)

# subset to count columns
tissue_counts <- tissue_count_data[ , grepl('count$', colnames(tissue_count_data)) &
                                      !grepl('normalised', colnames(tissue_count_data)) ]
names(tissue_counts) <- gsub('.count$', '', names(tissue_counts))
rownames(tissue_counts) <- tissue_count_data[['Gene.ID']]

tissue_samples <- data.frame(
  row.names = c('Head', 'Heart', 'Somites', 'PSC', 'CEM', 'Carcass'),
  condition = c('Head', 'Heart', 'Somites', 'PSC', 'CEM', 'Carcass'),
  Type = rep('Tissue', 6)
)

# Subset and reorder count data
tissue_counts <- tissue_counts[, row.names(tissue_samples)]

################################################################################
# merge counts together
merged_counts <- merge(sum_counts_by_stage, tissue_counts, by = "row.names")
rownames(merged_counts) <- merged_counts[, 'Row.names']
merged_counts <- merged_counts[ , !grepl('Row.names', colnames(merged_counts)) ]

# make merged samples df
merged_samples <- rbind(samples_by_stage, tissue_samples)

# make DESeq2 object and normalise
dds <- DESeqDataSetFromMatrix(merged_counts, merged_samples, design = ~ 1)
dds <- estimateSizeFactors(dds)
mcols(dds) <- rowData(Mm_GRCm38_e88_baseline)

norm_counts <- counts(dds, normalized = TRUE)

# count (protein-coding) genes expressed at >= 10 counts
norm_counts_protein <- norm_counts[ mcols(dds)$Biotype == 'protein_coding', ]

we_gt10 <- as.data.frame(norm_counts_protein[, merged_samples$Type == 'Whole Embryo'] >= 10)
we_gt10_genes <- rownames(norm_counts_protein)[ Reduce('|', we_gt10) ]

tissue_gt10 <- as.data.frame(norm_counts_protein[, merged_samples$Type == 'Tissue'] >= 10)
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

# plot heatmap of tissue only genes
tissue_only_genes <- setdiff(tissue_gt10_genes, we_gt10_genes)

tissue_only_counts <- norm_counts[tissue_only_genes, ]

# sort by row sums
row_sums <- rowSums(tissue_only_counts)
tissue_only_counts <- tissue_only_counts[names(sort(row_sums)), ]

# mean center and scale, melt
#tissue_only_counts_scaled <- t(scale(t(tissue_only_counts)))
tissue_only_counts_m <- melt(tissue_only_counts)
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
    scale_fill_viridis(discrete=TRUE, name = 'Counts') +
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
    width = 6, height = 9, paper = "special")
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

write.table(highly_expressed,
            file = file.path('output', 'highly_expressed-tissues.tsv'),
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

################################################################################
# save plot objects
save.image(file = file.path(wd, 'output', 'fig1.RData'))

################################################################################
