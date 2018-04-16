#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'fig5c.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 2
)

#cmd_line_args <- list(
#  options = list(verbose = FALSE),
#  args = c('data/fig5c_repeats_de.tsv',    
#           'data/fig5c_repeats_location.tsv')
#)

if( cmd_line_args$options[['verbose']] ){
}

packages <- c('ggplot2', 'reshape2', 'cowplot', 'grid', 'miscr')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

colour_palette <- 
  c('blue'  = '#A1BAE2',
    'red'   = '#C1504D',
    'green' = '#9FC257')

# read in data and set levels of Family factor
repeats_de <- read.delim(file = cmd_line_args$args[1])
repeats_de$Family <- factor(repeats_de$Family, levels = rev(unique(repeats_de$Family)))

# melt data
repeats_de_m <- melt(repeats_de[, !grepl('padj', names(repeats_de))],
                        id.vars = c('Family', 'full_length'),
                        variable.name = 'category', value.name = 'count')

# set levels of factors
repeats_de_m$category <- factor(repeats_de_m$category, levels = c('DE', 'notDE'))

# create colour palette for DE vs notDE
category_palette <- colour_palette[c('red', 'blue')]
names(category_palette) <- levels(repeats_de_m$category)

# create df for labels
repeats_de$total <- repeats_de$notDE + repeats_de$DE

# create df for lines
lines_df <- data.frame(
    'posn' = seq(0.5, nlevels(repeats_de$Family) - 0.5, 1)
)

# get max count to set limits on full-length plot
max_y <- ceiling( max( repeats_de$total, na.rm = TRUE ) / 1000 ) * 1000
all_bar_plot <- ggplot(data = repeats_de_m[repeats_de_m$full_length == "all", ]) +
                    geom_col(aes(x = Family, y = count, fill = category)) +
                    geom_text(data = repeats_de[ repeats_de$full_length == "all", ],
                                aes(x = Family, y = total, label = DE),
                                hjust = 0, nudge_y = 100) +
                    geom_vline(data = lines_df, aes(xintercept = posn),
                               linetype = 5, colour = 'grey80') +
                    scale_y_continuous(position = 'right', limits = c(0, max_y)) +
                    scale_fill_manual(values = category_palette) +
                    coord_flip() +
                    theme_void() +
                    theme(axis.text.x = element_text(colour = 'black', size = 12,
                                                     angle = 0, debug = FALSE),
                          axis.line.x = element_line(colour = 'black'),
                          axis.ticks.x = element_line(colour = 'black'),
                          legend.position = 'none')

pdf(file = file.path('plots', 'fig5c-all_repeats_de.pdf'))
print(all_bar_plot)
invisible(dev.off())

# To make the x scales match up for the all and full-length repeats plots
# remove the y axis labels and make a separate plot of the just the y axis labels
y_labels_plot <- ggplot(data = repeats_de_m[repeats_de_m$full_length == "all", ]) +
    geom_text(aes(x = full_length, y = Family, label = NA)) +
    theme_void() +
    theme(axis.text.y = element_text(colour = 'black', size = 12,
                                        angle = 0, debug = FALSE))

# plot of full_length data
full_length_bar_plot <- ggplot(data = repeats_de_m[repeats_de_m$full_length == "full_length", ]) +
                geom_col(aes(x = Family, y = count, fill = category)) +
                geom_text(data = repeats_de[ repeats_de$full_length == "full_length", ],
                            aes(x = Family, y = total, label = DE),
                            hjust = 0, nudge_y = 100) +
                geom_vline(data = lines_df, aes(xintercept = posn),
                           linetype = 5, colour = 'grey80') +
                scale_y_continuous(position = 'right', limits = c(0, max_y)) +
                scale_fill_manual(values = category_palette) +
                coord_flip()
                
pdf(file = file.path('plots', 'fig5c-full_length-repeats_de.pdf'))
print(full_length_bar_plot)
invisible(dev.off())

# remove legend and plot separately
repeats_legend <- get_gg_legend(full_length_bar_plot)

# plot legend
postscript(file = file.path('plots', 'repeats_de-legend.eps'),
    width = 2, height = 2, paper = 'special', horizontal = FALSE)
grid.draw(repeats_legend) 
invisible(dev.off())

# add theme to plot
full_length_bar_plot <- full_length_bar_plot +
                theme_void() +
                theme(axis.text.x = element_text(colour = 'black', size = 12,
                                                angle = 0, debug = FALSE),
                      axis.line.x = element_line(colour = 'black'),
                      axis.ticks.x = element_line(colour = 'black'),
                      legend.position = 'none')


# add in adjusted pvalues for enrichment tests
all_padj_data <- repeats_de[ repeats_de$full_length == "all", c("Family", "padj")]
# format numbers
all_padj_data$padj_formatted <- sprintf('%.1e', all_padj_data$padj)
all_padj_data$padj_formatted[ all_padj_data$padj_formatted == '0.0e+00' ] <- '0'
all_padj_data$x = factor(rep('all', nrow(all_padj_data)))

all_pvalue_plot <- ggplot(data = all_padj_data) +
                    geom_tile(aes(x = x, y = Family), fill = 'grey90') +
                    geom_text(aes(x = x, y = Family, label = padj_formatted),
                              hjust = 1, nudge_x = 0.25) +
                    geom_hline(data = lines_df, aes(yintercept = posn),
                               linetype = 5, colour = 'grey80') +
                    theme_void()

pdf(file = file.path('plots', 'fig5c-all-pvalues.pdf'))
print(all_pvalue_plot)
invisible(dev.off())

# same for full_length ones
full_length_padj_data <- repeats_de[ repeats_de$full_length == "full_length", c("Family", "padj")]
full_length_padj_data$padj_formatted <- sprintf('%.1e', full_length_padj_data$padj)
full_length_padj_data$padj_formatted[ full_length_padj_data$padj_formatted == '0.0e+00' ] <- '0'
full_length_padj_data$padj_formatted[ full_length_padj_data$padj_formatted == "NA" ] <- NA
full_length_padj_data$x = factor(rep('full_length', nrow(full_length_padj_data)))

# make a factor for the background fill
full_length_padj_data$fill_colour <- rep('empty', nrow(full_length_padj_data))
full_length_padj_data$fill_colour[ !is.na(full_length_padj_data$padj) ] <- 'fill'
full_length_padj_data$fill_colour <- factor(full_length_padj_data$fill_colour,
                                            levels = c('fill', 'empty'))
fill_palette <- c('fill' = 'grey90', 'empty' = 'white')

full_length_pvalue_plot <- ggplot(data = full_length_padj_data) +
                                geom_tile(aes(x = x, y = Family, fill = fill_colour)) +
                                geom_text(aes(x = x, y = Family, label = padj_formatted),
                                          hjust = 1, nudge_x = 0.25) +
                                scale_fill_manual(values = fill_palette) +
                                theme_void() +
                                theme(legend.position = 'none')

pdf(file = file.path('plots', 'fig5c-full_length-pvalues.pdf'))
print(full_length_pvalue_plot)
invisible(dev.off())

# exon, intron, intergenic plot
# read in data and set levels of Family
repeats_location <- read.delim(file = cmd_line_args$args[2])
repeats_location$Family <- factor(repeats_location$Family,
                                  levels = rev(unique(repeats_location$Family)))
# melt
repeats_location_m <- melt(repeats_location, id.vars = 'Family',
                           variable.name = 'Location',
                           value.name = 'count')
# set levels of Location
repeats_location_m$Location <- factor(repeats_location_m$Location,
                                    levels = rev(levels(repeats_location_m$Location)))

# create colour palette for Location
location_palette <- colour_palette[c('blue', 'red', 'green')]
names(location_palette) <- levels(repeats_location_m$Location)

# create repeats location plot
repeats_location_plot <- ggplot(data = repeats_location_m,
                                aes(x = Family, y = count, fill = Location)) +
                            geom_tile(height = 0.1, position = position_dodge()) +
                            geom_vline(data = lines_df, aes(xintercept = posn),
                                       linetype = 5, colour = 'grey80') +
                            scale_fill_manual(values = location_palette,
                                        guide = guide_legend(reverse = TRUE)) +
                            scale_y_log10(position = 'right') +
                            coord_flip()

pdf(file = file.path('plots', 'fig5c-all-location.pdf'))
print(repeats_location_plot)
invisible(dev.off())

# remove legend and plot at the end
location_legend <- get_gg_legend(repeats_location_plot)

# plot legend
postscript(file = file.path('plots', 'repeats_location-legend.eps'),
    width = 2, height = 2, paper = 'special', horizontal = FALSE)
grid.draw(location_legend) 
invisible(dev.off())

# add theme to plot
repeats_location_plot <- repeats_location_plot +
                            theme_void() +
                            theme(axis.text.x = element_text(colour = 'black', size = 12,
                                                angle = 0, debug = FALSE),
                                    axis.line.x = element_line(colour = 'black'),
                                    axis.ticks.x = element_line(colour = 'black'),
                                    legend.position = 'none')

# make version of the location plot with a break in the axis
# artificially set the 2 large values to something else
# then break the axis in Illustrator
repeats_location_m$count[ repeats_location_m$Family == 'L1MdGf_I' &
                            repeats_location_m$Location == 'intergenic' ] <- 125

repeats_location_plot_broken_axis <-
    ggplot(data = repeats_location_m) +
        geom_tile(aes(x = Family, y = count, fill = Location),
                    height = 4, position = position_dodge()) +
        geom_vline(data = lines_df, aes(xintercept = posn),
                    linetype = 5, colour = 'grey80') +
        scale_fill_manual(values = location_palette,
                            guide = guide_legend(reverse = TRUE)) +
        scale_y_continuous(position = 'right') +
        coord_flip() +
        theme_void() +
        theme(axis.text.x = element_text(colour = 'black', size = 12,
                angle = 0, debug = FALSE),
                axis.line.x = element_line(colour = 'black'),
                axis.ticks.x = element_line(colour = 'black'),
                legend.position = 'none')


postscript(file = file.path('plots', 'fig5c-all-location.eps'),
           width = 2, height = 9, paper = 'special', horizontal = FALSE)
print(repeats_location_plot_broken_axis)
invisible(dev.off())


# make composite figure
save_plot(file.path('plots', "morc2a-repeats.eps"),
          plot_grid(y_labels_plot, all_bar_plot, all_pvalue_plot,
                    repeats_location_plot_broken_axis, full_length_bar_plot,
                    full_length_pvalue_plot,
                    ncol = 6, rel_widths = c(3,6,2,4,6,2), align = 'h' ),
          ncol = 6, device = 'eps',
          base_height = 9, base_width = 2)



