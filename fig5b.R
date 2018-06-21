#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
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
#  options = list(verbose = FALSE),
#  args = c('data/repeats-location_vs_genes.tsv')
#)

if( cmd_line_args$options[['verbose']] ){
}

packages <- c('ggplot2', 'reshape2')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# read in data
repeats_de_location <- read.delim(file = cmd_line_args$args[1])
# melt with only exon, intron and intergenic columns
repeats_de_location_m <- melt(repeats_de_location,
                              id.vars = c('KO'),
                              measure.vars = c('exon', 'intron', 'intergenic'),
                              variable.name = 'Location',
                              value.name = 'count')

# reverse levels of Location and KO
repeats_de_location_m$Location <- factor(repeats_de_location_m$Location,
                                         levels = rev(as.character(levels(repeats_de_location_m$Location))))
repeats_de_location_m$KO <- factor(repeats_de_location_m$KO,
                                         levels = rev(as.character(levels(repeats_de_location_m$KO))))

# this is for plotting the DGE on the same plot
# transforming the data by 10-fold is totally arbitrary but happens to allow
# the two different sorts of data to be plotted together
repeats_de_location$DGE_transform <- repeats_de_location$DGE / 10

colour_palette <- c(
    'blue'  = "#0073B3",
    'purple' = '#CC99B3',
    'vermillion' = '#CC6600',
    'blue_green' = '#009980'
)

bars_palette <- colour_palette[c('purple', 'vermillion', 'blue_green')]
names(bars_palette) <- levels(repeats_de_location_m$Location)

location_bar_plot <- ggplot(data = repeats_de_location_m) +
    geom_col(aes(x = KO, y = count, fill = Location),
             position = "dodge") +
    geom_point(data = repeats_de_location, aes(x = KO, y = DGE_transform),
               shape = 4, size = 3, colour = colour_palette['blue']) +
    scale_y_continuous(name = "Number of repeats", position = 'right',
                       sec.axis = sec_axis(~ . * 10, name = "Number of differentially expressed genes") ) +
    scale_fill_manual(values = bars_palette,
                      guide = guide_legend(reverse = TRUE)) +
    coord_flip() +
    labs(x = 'KO Line') +
    theme_minimal(base_size = 15) +
    theme(legend.position = 'top',
          legend.direction = 'horizontal',
          title = element_text(size = 14),
          axis.text.y = element_text(face = "italic")
          )

pdf(file = file.path('plots', 'fig5b.pdf'))
print(location_bar_plot)
dev.off()

postscript(file = file.path('plots', 'fig5b.eps'),
           width = 5.5, height = 10, horizontal = FALSE, paper = 'special')
print(location_bar_plot)
dev.off()
