#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'fig5d.R',
    usage = "Usage: %prog [options] data_file heatmap_rda_file" ),
  positional_arguments = 2
)

#cmd_line_args <- list(
#  options = list(directory = 'cwd',
#                 plot_file = file.path('plots', 'fig5d.svg'),
#                 verbose = FALSE),
#  args = c('data/fig5d_data_go.tsv',
#           'output/fig5d-heatmap.rda')
#)

packages <- c('ggplot2', 'reshape2', 'cowplot')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

go_data <- read.delim(cmd_line_args$args[1], header = TRUE)
# reshape to long format
go_data_m <- melt(go_data, id.vars = c('Gene', 'Name'),
                  measure.vars = c('CNS', 'cardiac', 'microtubule'),
                  variable.name = 'GO_category')
# set levels of Gene and GO_category
go_data_m$Gene <- factor(go_data_m$Gene,
                         levels = rev(unique(as.character(go_data_m$Gene))))
go_data_m$Name <- factor(go_data_m$Name,
                         levels = rev(unique(as.character(go_data_m$Name))))
go_data_m$GO_category <- factor(go_data_m$GO_category,
                         levels = unique(as.character(go_data_m$GO_category)))

# subset to where value == 1
subset_go_data_m <- go_data_m[ go_data_m$value == 1, ]

colour_blind_palette <- 
  c( 'blue' = rgb(0,0.45,0.7),
     'yellow' = rgb(0.95, 0.9, 0.25),
     'sky_blue' = rgb(0.35, 0.7, 0.9),
     'purple' = rgb(0.8, 0.6, 0.7),
     'orange' = rgb(0.9, 0.6, 0),
     'vermillion' = rgb(0.8, 0.4, 0),
     'blue_green' = rgb(0, 0.6, 0.5),
     'black' = rgb(0, 0, 0) )

category_palette <- colour_blind_palette[c('blue', 'purple', 'orange')]
names(category_palette) <- levels(go_data_m$GO_category)

category_plot <- ggplot(data = subset_go_data_m) +
  geom_tile(aes(x = GO_category, y = Gene, fill = GO_category)) +
  scale_x_discrete(position = 'top') +
  scale_fill_manual(values = category_palette) +
  theme_void() +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 45,
                                   hjust = 0, vjust = 0.5, debug = FALSE),
        legend.text = element_text(size = 12, colour = 'black'),
        legend.title = element_text(size = 14, colour = 'black'),
        legend.position = 'bottom')

# load heatmap and plot together
load(cmd_line_args$args[2])

heatmap_plot <- heatmap_plot +
  theme(axis.text.x = element_text(size = 12, colour = 'black', angle = 45,
                                   hjust = 0, vjust = 0.5, debug = FALSE),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.title = element_blank())

# make composite figure
save_plot(file.path('plots', "fig5d.eps"),
          plot_grid(heatmap_plot, category_plot,
          nrow = 1, ncol = 2, rel_widths = c(4,1), align = 'h', axis = 'tb'),
          ncol = 2, device = 'eps',
          base_height = 4, base_aspect_ratio = 1 )


