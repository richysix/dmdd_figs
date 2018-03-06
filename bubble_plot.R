#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option("--output_file", type="character", default='bubble_plot.pdf',
              help="Name of output file [default %default]" ),
  make_option("--output_type", type="character", default='pdf',
              help="Image Format of output file (pdf, png, ps, svg) [default %default]" ),
  make_option("--x_var", type="character", default='x',
              help="Column to plot on the x axis [default %default]" ),
  make_option("--y_var", type="character", default='y',
              help="Column to plot on the y axis [default %default]" ),
  make_option("--size_var", type="character", default='size',
              help="Column to map to size [default %default]" ),
  make_option("--fill_var", type="character", default='fill',
              help="Column to map to fill [default %default]" ),
  make_option("--x_labels", type="character", default=NULL,
              help="Column to use as labels for the x axis [default %default]" ),
  make_option("--y_labels", type="character", default=NULL,
              help="Column to use as labels for the y axis [default %default]" ),
  make_option("--reverse_y", action="store_true", default=FALSE,
              help="Switch to reverse the ordering of the y axis [default %default]" ),
  make_option("--width", type="integer", default=NULL,
              help="Width of the plot [default %default]" ),
  make_option("--height", type="integer", default=NULL,
              help="Height of the plot [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'bubble_plot.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 1
)

#cmd_line_args <- list(
#  options = list(directory = 'cwd', output_file = 'plots/topgo-shared-terms-bubble.pdf',
#                 x_var = 'Gene', y_var = 'GO.ID', size_var = 'Fold Enrichment', fill_var = '-log10(pvalue)',
#                 x_labels = NULL, y_labels = NULL, reverse_y = TRUE,
#                 verbose = TRUE),
#  args = c('data/topgo-shared-terms.tsv')
#)

# convert options to less cumbersome variable names
x_col <- cmd_line_args$options[['x_var']]
y_col <- cmd_line_args$options[['y_var']]
size_col <- cmd_line_args$options[['size_var']]
fill_col <- cmd_line_args$options[['fill_var']]
x_lab_col <- cmd_line_args$options[['x_labels']]
y_lab_col <- cmd_line_args$options[['y_labels']]
reverse_y <- cmd_line_args$options[['reverse_y']]
output_file <- cmd_line_args$options[['output_file']]
output_type <- cmd_line_args$options[['output_type']]
width <- cmd_line_args$options[['width']]
height <- cmd_line_args$options[['height']]

if( cmd_line_args$options[['verbose']] ){
  cat( "Working directory:",
      cmd_line_args$options[['directory']], "\n", sep=" " )
}

if (cmd_line_args$options[['directory']] == 'cwd') {
  wd <- getwd()
} else {
  wd <- cmd_line_args$options[['directory']]
}

packages <- c('ggplot2', 'viridis', 'reshape2', 'devtools')
for( package in packages ){
  suppressPackageStartupMessages(
    suppressWarnings( library(package, character.only = TRUE) )
  )
}
# load biovisr. If not installed, install from github repo
package <- 'biovisr'
if(!require( package, character.only = TRUE )){
  install_github('richysix/biovisr')
  library( package, character.only = TRUE )
}

################################################################################
#' bubble_plot
#'
#' \code{bubble_plot} create 'bubble' plot
#' 
#' This function takes a data.frame and creates a bubble plot.
#' An x-y scatterplot where an attribute is represented by the size of the
#' points and optionally an another attribute is mapped to the colour.
#' 
#' @param plot_df   data.frame  data to plot
#' @param x         character   Name of the variable to plot on the x-axis
#' @param y         character   Name of the variable to plot on the y-axis
#' @param size      character   Name of the variable to use for the size of the points
#' @param fill      character   Name of the variable to use for the fill colour
#' @param x_labels  character   Optional labels to use instead of levels of x-axis variable
#' @param y_labels  character   Optional labels to use instead of levels of y-axis variable
#'
#' @return plot - ggplot2 object
#'
#bubble_plot <- function(plot_df, x = 'x', y = 'y', size = 'size',
#                        fill = 'fill', x_labels = NULL,
#                        y_labels = NULL ){
#  # check labels is the same length as the levels of the x/y column
#  if (!is.null(x_labels) ) {
#    if (length(x_labels) != nlevels(plot_df[[x]]) ) {
#      stop('Supplied labels vector for x axis is the wrong length')
#    }
#  }
#  if (!is.null(y_labels) ) {
#    if (length(y_labels) != nlevels(plot_df[[y]]) ) {
#      stop('Supplied labels vector for y axis is the wrong length')
#    }
#  }  
#  bubble_plot <- ggplot(data = plot_df) +
#    geom_point(aes_(x = as.name(x), y = as.name(y), size = as.name(size),
#                   fill = as.name(fill)), shape = 21) +
#    scale_fill_viridis(direction = -1) +
#    theme_void() +
#    theme(
#      axis.text.x = element_text(size = 10, colour = 'black', angle = 90,
#                                  hjust = 0, debug = FALSE),
#      axis.text.y = element_text(size = 8, colour = 'black', angle = 0,
#                                  hjust = 0, debug = FALSE),
#      panel.grid.major = element_line(colour = 'grey80', linetype = 'dotted'),
#      legend.position = 'top',
#      legend.title = element_text(size = 14),
#      legend.text = element_text(size = 12))
#  
#  if( !is.null(x_labels) ) {
#    bubble_plot <- bubble_plot +
#      scale_x_discrete(position = 'top', labels = x_labels)
#  } else {
#    bubble_plot <- bubble_plot +
#      scale_x_discrete(position = 'top')
#  }
#
#  if( !is.null(y_labels) ) {
#    bubble_plot <- bubble_plot +
#      scale_y_discrete(labels = y_labels)
#  }
#  
#  return(bubble_plot)
#}

# load input data
input_data <- read.delim(cmd_line_args$args[1], stringsAsFactors = FALSE,
                         check.names = FALSE)
# set levels of input_data
# check class of x and y data
# if character, convert to factors in the order that they appear in the file
if( class(input_data[[ x_col ]]) == 'character' ) {
  input_data[[ x_col ]] <- factor(input_data[[ x_col ]],
                                levels = unique(input_data[[ x_col ]]))
}
# reverse the order for the y axis if option --reverse_y is TRUE
if( class(input_data[[ y_col ]]) == 'character' ) {
  if( reverse_y ) {
    y_levels <- rev(unique(input_data[[ y_col ]]))
  } else {
    y_levels <- unique(input_data[[ y_col ]])    
  }
  input_data[[ y_col ]] <- factor(input_data[[ y_col ]],
                                levels = y_levels)
}

# check labels
if ( is.null(x_lab_col) ) {
  x_labels <- NULL 
} else {
  x_labels <- unique(input_data[[ x_lab_col ]])
}

# Need to make sure the y's and y labels are the same length
if ( is.null(y_lab_col) ) {
  y_labels <- NULL 
} else {
  y_lab_df <- unique(input_data[, c(y_col, y_lab_col) ])
  if( reverse_y ) {
    y_labels <- rev(y_lab_df[[ y_lab_col ]])
  } else {
    y_labels <- y_lab_df[[ y_lab_col ]]
  }
}

plot <- bubble_plot(input_data, x = x_col,
                    y = y_col,
                    size = size_col,
                    fill = fill_col,
                    x_labels = x_labels,
                    y_labels = y_labels)

# plot to file
if (output_type == "png") {
  if (is.null(width)) {
    width <- 480
  }
  if (is.null(height)) {
    height <- 480
  }
  png(file = output_file, width = width, height = height)
} else if (output_type == "ps") {
  if (is.null(width)) {
    width <- 10
  }
  if (is.null(height)) {
    height <- 11
  }
  postscript(file = output_file, width = width, height = height,
             paper = "special", horizontal = FALSE)
} else if (output_type == "svg") {
  if (is.null(width)) {
    width <- 10
  }
  if (is.null(height)) {
    height <- 11
  }
  library(svglite)
  svglite(file = output_file, width = width, height = height)
} else {
  if (is.null(width)) {
    width <- 10
  }
  if (is.null(height)) {
    height <- 11
  }
  pdf(file = output_file, width = width, height = height)
}
print(plot)
dev.off()
