#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'fig7c.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 3
)

#cmd_line_args <- list(
#  options = list(verbose = FALSE),
#  args = c('output/fig7c-for-enrichment-test.tsv',
#           1639267,
#           'output/fig7c-dhx35-repeats-introns-genes-sig.tsv')
#)

if( cmd_line_args$options[['verbose']] ){
}

packages <- c('ggplot2', 'viridis', 'ggrepel')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# read in data
repeats_de_introns <- read.delim(file = cmd_line_args$args[1])
# replace blacklist with NA and convert gene_padj and gene_log2fc to numeric
repeats_de_introns$gene_padj[ repeats_de_introns$gene_padj == 'blacklist' ] <- NA
repeats_de_introns$gene_padj <- as.numeric(levels(repeats_de_introns$gene_padj))[repeats_de_introns$gene_padj]

repeats_de_introns$gene_log2fc[ repeats_de_introns$gene_log2fc == 'blacklist' ] <- NA
repeats_de_introns$gene_log2fc <- as.numeric(levels(repeats_de_introns$gene_log2fc))[repeats_de_introns$gene_log2fc]

# test for enrichment
test_enrichment <- function(i, data, total_repeats) {
    binom_res <- binom.test(data$repeats_intron_de[i], sum(data$repeats_intron_de),
                            p = data$repeats_intron[i]/total_repeats,
                            alternative = 'greater')
    return(binom_res$p.value)
}

total_repeats <- as.integer( cmd_line_args$args[2] )
#total_repeats_de <- c(all = 1293, full_length = 592)
repeats_de_introns$p.value <- sapply(1:nrow(repeats_de_introns), test_enrichment,
                                        repeats_de_introns, total_repeats )

# adjust p values with Benjamini-Hochberg
repeats_de_introns$padj <- p.adjust(repeats_de_introns$p.value, method = "BH")
# order by padj
repeats_de_introns <- repeats_de_introns[ order(repeats_de_introns$padj), ]
# subset to genes with a significant enrichment and more than 1 repeat
# also remove blacklist ones
gene_sig_enrichment <- repeats_de_introns[ repeats_de_introns$padj < 0.05 &
                                          repeats_de_introns$repeats_intron_de > 1 &
                                          !is.na(repeats_de_introns$gene_padj), ]

fold_change_data <- read.delim(file = cmd_line_args$args[3])

# get fold change of repeats in de in each gene
gene_sig_enrichment$repeat_mean_log2fc <-
    sapply(as.character(gene_sig_enrichment$gene_id),
            function(gene_id, fold_change_data){
                mean(fold_change_data[ fold_change_data$gene_id == gene_id, 'log2fc'])
            }, fold_change_data )

gene_sig_enrichment$gene_padj_na <- gene_sig_enrichment$gene_padj
gene_sig_enrichment$gene_padj_na[ gene_sig_enrichment$gene_padj_na > 0.05 ] <- NA
gene_sig_enrichment$gene_log10p <- -log10(gene_sig_enrichment$gene_padj_na)

gene_sig_enrichment$repeat_enrich_log10p <- -log10(gene_sig_enrichment$padj)

gene_sig_enrichment$gene_name <-
    sapply(as.character(gene_sig_enrichment$gene_id),
            function(gene_id, fold_change_data){
                as.character(fold_change_data[ fold_change_data$gene_id == gene_id, 'gene_name' ][1])
            }, fold_change_data )

## label genes above 2.5 log2fc
#gene_sig_enrichment$label <- paste0(gene_sig_enrichment$gene_name, ' (',
#                                     gene_sig_enrichment$repeats_intron_de, ')')
#genes_to_label <- c('ENSMUSG00000038235', 'ENSMUSG00000046574', 'ENSMUSG00000031555',
#                    'ENSMUSG00000028811', 'ENSMUSG00000020547', 'ENSMUSG00000027012',
#                    'ENSMUSG00000070544', 'ENSMUSG00000037236')
#gene_sig_enrichment$label[ !(Reduce('|', lapply(genes_to_label,
#                                                function(gene_id){
#                                                   gene_sig_enrichment$gene_id == gene_id
#                                                }
#                            ))) ] <- ''

#gene_sig_enrichment$label[ !(gene_sig_enrichment$repeat_mean_log2fc > 3) ] <- ""

# calculate pearson correlation coefficient for Dhx35_plus ones
Dhx35_plus <- gene_sig_enrichment[ gene_sig_enrichment$source == 'Dhx35_plus', ]
# subset to ones with a +ve log2fc and pval < 0.05
Dhx35_plus <- Dhx35_plus[ Dhx35_plus$gene_padj < 0.05 & Dhx35_plus$gene_log2fc > 0, ]

pearson_coeff <- cor(Dhx35_plus$gene_log2fc, Dhx35_plus$repeat_mean_log2fc)
print(sprintf('Pearson correlation coefficient (Dhx35 plus others) = %.3f', pearson_coeff))

Dhx35_only <- gene_sig_enrichment[ gene_sig_enrichment$source == 'Dhx35_only', ]
# subset to ones with a +ve log2fc and pval < 0.05
Dhx35_only <- Dhx35_only[ Dhx35_only$gene_padj < 0.05 & Dhx35_only$gene_log2fc > 0, ]

pearson_coeff_dhx35 <- cor(Dhx35_only$gene_log2fc, Dhx35_only$repeat_mean_log2fc)
print(sprintf('Pearson correlation coefficient (Dhx35 only) = %.3f', pearson_coeff_dhx35))

# plot repeat mean log2fc against gene log2fc
repeats_log2fc_plot <- ggplot(data = gene_sig_enrichment,
                              aes(x = gene_log2fc, y = repeat_mean_log2fc)) +
    geom_point(aes(shape = source, colour = gene_log10p), size = 2) +
    #geom_text_repel(aes(label = label), size = 4.2, fontface = "italic") +
    scale_color_viridis(name = expression(-log[10]*'['*pvalue*']'),
                        direction = -1, na.value = 'grey80') +
    scale_shape(name = '', labels = c('Dhx35 only', 'Dhx35 plus \nother(s)')) +
    labs(x = expression(paste(log[2], '[Fold Change] (Gene)', sep = "")),
         y = expression(paste('Mean ', log[2], '[Fold Change] (Repeat)', sep = "")) ) +
    theme_minimal() +
    theme(title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12)
          )

#repeats_log2fc_plot_enrich_p <-
#    ggplot(data = gene_sig_enrichment,
#            aes(x = gene_log2fc, y = repeat_mean_log2fc)) +
#    geom_point(aes(shape = source, colour = repeat_enrich_log10p),
#               size = 2) +
#    scale_color_viridis(direction = -1) +
#    geom_text_repel(aes(label = label)) +
#    theme_minimal()

pdf(file = file.path('plots', 'fig7c-plots.pdf'),
    width = 11, height = 8, paper = 'special')
print(repeats_log2fc_plot)
#print(repeats_log2fc_plot_enrich_p)
dev.off()

# output eps
postscript(file = file.path('plots', 'fig7c.eps'),
            width = 9, height = 5, paper = 'special',
            horizontal = FALSE)
print(repeats_log2fc_plot)
dev.off()

# output gene ids to file
write.table(gene_sig_enrichment, file = file.path('output', 'fig7c-plot_data.tsv'),
            sep = "\t",  quote = FALSE, row.names = FALSE)
