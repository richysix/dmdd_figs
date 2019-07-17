#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option("--output_file", type="character", default='plots/fig4-heatmap.svg',
              help="Name of output file [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'fig5b.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 1
)

#cmd_line_args <- list(
#  options = list(directory = 'cwd',
#                 output_file = file.path('plots', 'fig5b-heatmap.svg'),
#                 verbose = FALSE),
#  args = c('data/fig5b_log2fc.tsv')
#)

if (cmd_line_args$options[['directory']] == 'cwd') {
  wd <- getwd()
} else {
  wd <- cmd_line_args$options[['directory']]
}

if( cmd_line_args$options[['verbose']] ){
  cat( "Working directory:", cmd_line_args$options[['directory']], "\n", sep=" " )
  cat( "Output file:", cmd_line_args$options[['output_file']], "\n", sep=" " )
}

packages <- c('ggplot2', 'scales', 'plyr', 'svglite')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# load data
input_data <- read.table(file = cmd_line_args$args[1], sep = "\t",
                         header = TRUE)
# if header exists match column names by name, else match by number
x_column <- 'Line'
y_column <- 'Gene.ID'
data_column <- 'log2fc'

# set levels of Line
input_data[[x_column]] <- factor(input_data[[x_column]],
                           levels = unique(as.character(input_data[[x_column]])) )

# order genes
num_appearances <- ddply(input_data[ !is.na(input_data$log2fc), ], .(Gene.ID),
                         summarise, times = length(Line))
input_data <- merge(input_data, num_appearances)
# set levels of class
input_data$Class <- factor(input_data$Class,
                           levels = c('Downstream of Shh signalling', 'Shh signalling interactors', 'Novel'))

# sort df based on Class and times
input_data <- input_data[ order(match(input_data$Class, levels(input_data$Class)), -input_data$times), ]

# set levels of Gene.ID
input_data[[y_column]] <- factor(input_data[[y_column]],
                           levels = rev(unique(as.character(input_data[[y_column]]))) )

genes_tmp <- unique(input_data[, c('Gene.ID', 'Gene.Name')])
genes <- as.character(genes_tmp$Gene.Name[ !is.na(genes_tmp$Gene.Name) ])
names(genes) <- genes_tmp$Gene.ID[ !is.na(genes_tmp$Gene.Name) ]

# create new data column
input_data$truncated_log2fc <- input_data$log2fc
# set everything above 2 to 2 and everything below -2 to -2
input_data$truncated_log2fc[ input_data$truncated_log2fc > 2 ] <- 2
input_data$truncated_log2fc[ input_data$truncated_log2fc < -2 ] <- -2

heatmap_plot <- ggplot(data = input_data) +
  geom_raster( aes_q(x = as.name(x_column), y = as.name(y_column),
                     fill = quote(truncated_log2fc) ) ) +
  guides( fill = guide_colourbar(title = expression(log[2]*"[Fold Change]") ) ) +
  scale_x_discrete(position = 'top') +
  scale_y_discrete(labels = genes) + 
  scale_fill_gradientn(limits = c(-2,1), colours = c('blue', 'white', 'red'),
                       values = rescale(c(-2,0,1), to = c(0,1)),
                       na.value = 'grey 80') +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_text(face = 'italic') )

svglite(file = cmd_line_args$options[['output_file']])
print(heatmap_plot)
dev.off()

postscript(file = sub('svg$', 'eps', cmd_line_args$options[['output_file']]),
           horizontal = FALSE, width = 7, height = 4)
print(heatmap_plot)
dev.off()


################################################################################
# save plot objects
save.image(file = file.path(wd, 'output', 'fig5.RData'))

################################################################################
