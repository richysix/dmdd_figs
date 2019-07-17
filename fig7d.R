#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'fig5d.R',
    usage = "Usage: %prog [options] repeats_file num_repeats num_repeats_gt_4kbp repeats_location_file" ),
  positional_arguments = 4
)

#cmd_line_args <- list(
#  options = list(verbose = FALSE),
#  args = c('output/fig5d_repeats_de.tsv',
#           3765374, 24162, 
#           'output/fig5d_repeats_location.tsv')
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
    'purple' = '#CC99B3',
    'vermillion' = '#CC6600',
    'blue_green' = '#009980')

# read in data
repeats_de_all <- read.delim(file = cmd_line_args$args[1])

# test for enrichment
test_enrichment <- function(i, data_subset, total_repeats, full_length) {
    binom_res <- binom.test(data_subset$de[i], sum(data_subset$de),
                            p = data_subset$repeats[i]/total_repeats[full_length],
                            alternative = 'greater')
    return(binom_res$p.value)
}

#total_repeats <- c(all = 3765374, full_length = 24162)
total_repeats <- c("all" = as.integer( cmd_line_args$args[2] ),
                   "full_length" = as.integer( cmd_line_args$args[3] ) )

pvals <- vector('list', length = nlevels(repeats_de_all$full_length))
list_index <- 1
for (full_length in levels(repeats_de_all$full_length)) {
    data_subset <- repeats_de_all[ repeats_de_all$full_length == full_length, ]
    pvals[[list_index]] <- sapply(1:nrow(data_subset), test_enrichment,
                          data_subset, total_repeats, full_length )
    list_index <- list_index + 1
}

# adjust p values with Benjamini-Hochberg
padj_list <- lapply(pvals, p.adjust, method = "BH")

repeats_de_all$p.value <- do.call(c, pvals)
repeats_de_all$padj <- do.call(c, padj_list)

enriched_repeats <- repeats_de_all[ repeats_de_all$padj < 0.05, ]
# set levels of Group
enriched_repeats$Group <- factor(enriched_repeats$Group,
                                 levels = c("LINE", "LTR", "DNA", "Satellite", "Unknown") )

enriched_repeats <- enriched_repeats[ order(enriched_repeats$Group, enriched_repeats$padj), ]
# set levels of Family
enriched_repeats$Family <- factor(enriched_repeats$Family,
                                  levels = rev(unique(enriched_repeats$Family)) )
# create entries for families that aren't enriched as full_length
enriched_repeats <-
    rbind( enriched_repeats,
            data.frame(
                Family = setdiff(enriched_repeats$Family[ enriched_repeats$full_length == "all" ],
                                enriched_repeats$Family[ enriched_repeats$full_length == "full_length" ]),
                Group = NA,
                full_length = "full_length",
                repeats = NA,
                de = NA,
                not_de = NA,
                p.value = NA,
                padj = NA
            )
    )

# write out repeat families that are enriched
write.table(enriched_repeats, file = file.path('output', 'repeats-enriched_families.tsv'),
            sep = "\t",  quote = FALSE, row.names = FALSE)


# melt data
repeats_de_sig_m <- melt(enriched_repeats[ , c('Group', 'Family', 'full_length', 'not_de', 'de')],
                            id.vars = c('Group', 'Family', 'full_length'),
                            variable.name = 'category', value.name = 'count')

# set levels of factors
repeats_de_sig_m$category <- factor(repeats_de_sig_m$category, levels = c('de', 'not_de'))

# create colour palette for DE vs notDE
category_palette <- colour_palette[c('red', 'blue')]
names(category_palette) <- levels(repeats_de_sig_m$category)

## create df for labels
#repeats_de$total <- repeats_de$notDE + repeats_de$DE

# create df for lines
lines_df <- data.frame(
    'posn' = seq(0.5, nlevels(repeats_de_sig_m$Family) - 0.5, 1)
)

# get max count to set limits on full-length plot
max_y <- ceiling( max( enriched_repeats$repeats, na.rm = TRUE ) / 1000 ) * 1000
# use max_y to set the limits on the full-length plot so that the bars are scaled the same
all_bar_plot <- ggplot(data = repeats_de_sig_m[ repeats_de_sig_m$full_length == 'all', ]) +
                    geom_col(aes(x = Family, y = count, fill = category)) +
                    geom_text(data = enriched_repeats[ enriched_repeats$full_length == 'all', ],
                                aes(x = Family, y = repeats, label = de),
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

pdf(file = file.path('plots', 'fig5d-all_repeats_de.pdf'))
print(all_bar_plot)
invisible(dev.off())

# To make the x scales match up for the all and full-length repeats plots
# remove the y axis labels and make a separate plot of the just the y axis labels
y_labels_plot <- ggplot(data = repeats_de_sig_m[ repeats_de_sig_m$full_length == 'all', ]) +
    geom_text(aes(x = full_length, y = Family, label = NA)) +
    theme_void() +
    theme(axis.text.y = element_text(colour = 'black', size = 12,
                                    angle = 0, hjust = 1, debug = FALSE))

# plot of full_length data
full_length_bar_plot <- ggplot(data = repeats_de_sig_m[ repeats_de_sig_m$full_length == 'full_length', ]) +
                geom_col(aes(x = Family, y = count, fill = category)) +
                geom_text(data = enriched_repeats[ enriched_repeats$full_length == 'full_length', ],
                            aes(x = Family, y = repeats, label = de),
                            hjust = 0, nudge_y = 100) +
                geom_vline(data = lines_df, aes(xintercept = posn),
                           linetype = 5, colour = 'grey80') +
                scale_y_continuous(position = 'right', limits = c(0, max_y)) +
                scale_fill_manual(values = category_palette) +
                coord_flip()
                
pdf(file = file.path('plots', 'fig5d-full_length-repeats_de.pdf'))
print(full_length_bar_plot)
invisible(dev.off())

# remove legend and plot separately
repeats_legend <- get_gg_legend(full_length_bar_plot + theme(legend.position = 'top', legend.direction = 'horizontal'))

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
all_padj_data <- enriched_repeats[ enriched_repeats$full_length == "all", c("Family", "padj")]
# format numbers
all_padj_data$padj_formatted <- sprintf('%.1e', all_padj_data$padj)
all_padj_data$padj_formatted[ all_padj_data$padj_formatted == '0.0e+00' ] <- '< 5e-324'
all_padj_data$x = factor(rep('all', nrow(all_padj_data)))

all_pvalue_plot <- ggplot(data = all_padj_data) +
                    geom_tile(aes(x = x, y = Family), fill = 'grey90') +
                    geom_text(aes(x = x, y = Family, label = padj_formatted),
                              hjust = 1, nudge_x = 0.45) +
                    geom_hline(data = lines_df, aes(yintercept = posn),
                               linetype = 5, colour = 'grey80') +
                    theme_void()

pdf(file = file.path('plots', 'fig5d-all-pvalues.pdf'))
print(all_pvalue_plot)
invisible(dev.off())

# same for full_length ones
full_length_padj_data <- enriched_repeats[ enriched_repeats$full_length == "full_length", c("Family", "padj")]
full_length_padj_data$padj[ full_length_padj_data$padj >= 0.05 ] <- NA
full_length_padj_data$padj_formatted <- sprintf('%.1e', full_length_padj_data$padj)
full_length_padj_data$padj_formatted[ full_length_padj_data$padj_formatted == '0.0e+00' ] <- '< 5e-324'
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
                                          hjust = 1, nudge_x = 0.45) +
                                scale_fill_manual(values = fill_palette) +
                                theme_void() +
                                theme(legend.position = 'none')

pdf(file = file.path('plots', 'fig5d-full_length-pvalues.pdf'))
print(full_length_pvalue_plot)
invisible(dev.off())

# exon, intron, intergenic plot
# read in data and set levels of Family
repeats_location <- read.delim(file = cmd_line_args$args[4])
# subset to significant families
rownames(repeats_location) <- repeats_location$Family
repeats_location <- repeats_location[unique( as.character(enriched_repeats$Family) ), ]
repeats_location$Family <- factor(repeats_location$Family,
                                  levels = levels(enriched_repeats$Family) )
# melt
repeats_location_m <- melt(repeats_location, id.vars = 'Family',
                           variable.name = 'Location',
                           value.name = 'count')
# set levels of Location
repeats_location_m$Location <- factor(repeats_location_m$Location,
                                    levels = rev(levels(repeats_location_m$Location)))

# create colour palette for Location
location_palette <- colour_palette[c('purple', 'vermillion', 'blue_green')]
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

pdf(file = file.path('plots', 'fig5d-all-location.pdf'))
print(repeats_location_plot)
invisible(dev.off())

# remove legend and plot at the end
location_legend <- get_gg_legend(repeats_location_plot + theme(legend.position = 'top', legend.direction = 'horizontal'))

# plot legend
postscript(file = file.path('plots', 'repeats_location-legend.eps'),
    width = 4, height = 2, paper = 'special', horizontal = FALSE)
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
# artificially set the largest value to something else
# then break the axis in Illustrator
repeats_location_m$count[ repeats_location_m$Family == 'L1MdGf_I' &
                            repeats_location_m$Location == 'intergenic' ] <- 125

repeats_location_plot_broken_axis <-
    ggplot(data = repeats_location_m) +
        geom_col(aes(x = Family, y = count, fill = Location),
                    position = position_dodge()) +
        geom_vline(data = lines_df, aes(xintercept = posn),
                    linetype = 5, colour = 'grey80') +
        scale_fill_manual(values = location_palette,
                            guide = guide_legend(reverse = TRUE)) +
        scale_y_continuous(position = 'right', breaks = c(0,50,100) ) +
        coord_flip() +
        theme_void() +
        theme(axis.text.x = element_text(colour = 'black', size = 12,
                angle = 0, debug = FALSE),
                axis.line.x = element_line(colour = 'black'),
                axis.ticks.x = element_line(colour = 'black'),
                legend.position = 'none')


postscript(file = file.path('plots', 'fig5d-all-location.eps'),
           width = 2, height = 9, paper = 'special', horizontal = FALSE)
print(repeats_location_plot_broken_axis)
invisible(dev.off())


# make composite figure
save_plot(file.path('plots', "morc2a-repeats.eps"),
          plot_grid(y_labels_plot, all_bar_plot, all_pvalue_plot,
                    repeats_location_plot_broken_axis, full_length_bar_plot,
                    full_length_pvalue_plot,
                    ncol = 6, rel_widths = c(6,9,3,6,9,3), align = 'h' ),
          ncol = 6, device = 'eps',
          base_height = 9, base_width = 1.5)
