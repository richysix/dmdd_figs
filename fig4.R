#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option("--output_file", type="character", default='plots/fig5-heatmap.pdf',
              help="Name of output file [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'fig5.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 1
)

#cmd_line_args <- list(
#  options = list(directory = 'cwd',
#                 output_file = file.path('plots', 'fig5-heatmap.svg'),
#                 verbose = FALSE),
#  args = c('output/fig5_table.tsv')
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

packages <- c('ggplot2', 'viridis', 'scales')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# load data
input_data <- read.table(file = cmd_line_args$args[1], sep = "\t",
                         header = TRUE)
# if header exists match column names by name, else match by number
x_column <- 'KO'
y_column <- 'Gene'
data_column <- 'log2fc'

# set levels of x and y as order in which they appear in input file
input_data[[x_column]] <-
  factor(input_data[[x_column]], levels = unique(input_data[[x_column]]))
input_data[[y_column]] <-
    factor(input_data[[y_column]],
      levels = rev(unique(input_data[[y_column]])))

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
  scale_fill_gradientn(limits = c(-2,1), colours = c('blue', 'white', 'red'),
                       values = rescale(c(-2,0,1), to = c(0,1)),
                       na.value = 'grey 80') +
  theme_minimal() +
  theme(axis.title = element_blank())

if(grepl('pdf$', cmd_line_args$options[['output_file']])) {
  pdf(file = cmd_line_args$options[['output_file']])
} else if (grepl('ps$', cmd_line_args$options[['output_file']])) {
  postscript(file = cmd_line_args$options[['output_file']],
             paper = "A4", horizontal = FALSE)
} else if (grepl('svg$', cmd_line_args$options[['output_file']])) {
  library(svglite)
  svglite(file = cmd_line_args$options[['output_file']])
} else { # default to pdf 
  pdf(file = cmd_line_args$options[['output_file']])
}

print(heatmap_plot)
dev.off()

################################################################################
# save plot objects
save.image(file = file.path(wd, 'output', 'fig5.RData'))

################################################################################
