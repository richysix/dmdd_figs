#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'SupplFig-repeats.R',
    usage = "Usage: %prog [options] input_files" ),
  positional_arguments = 6
)

wd <- getwd()

if( cmd_line_args$options[['verbose']] ){
}

packages <- c('ggplot2', 'plyr', 'reshape2', 'grid')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

#cmd_line_args <- list(
#  options = list(verbose = FALSE),
#  args = c('output/num_repeats_by_line.tsv',
#           'output/repeats_by_line_by_type.tsv',
#           'output/repeats-for_volcano_plot.tsv',
#           '/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Morc2a/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv',
#           '/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Morc2a/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.tsv',
#           'output/repeats-enriched_families.tsv')
#)

colour_blind_palette <- 
  c( 'blue' = rgb(0,0.45,0.7),
     'yellow' = rgb(0.95, 0.9, 0.25),
     'sky_blue' = rgb(0.35, 0.7, 0.9),
     'purple' = rgb(0.8, 0.6, 0.7),
     'orange' = rgb(0.9, 0.6, 0),
     'vermillion' = rgb(0.8, 0.4, 0),
     'blue_green' = rgb(0, 0.6, 0.5),
     'black' = rgb(0, 0, 0) )

## Fig S4a
# load number of repeats data
num_repeats_by_line <- read.delim(file = cmd_line_args$args[1])
num_repeats_by_line <-
  num_repeats_by_line[ order(num_repeats_by_line$Num_repeats, decreasing = TRUE), ]
# set levels of KO_Line
num_repeats_by_line$KO_line <- factor(num_repeats_by_line$KO_line,
                                      levels = num_repeats_by_line$KO_line)
# make a factor for colour of text
num_repeats_by_line$text_colour <-
  factor(ifelse(num_repeats_by_line$Num_repeats > 100 | num_repeats_by_line$KO_line == 'Hmgxb3',
         TRUE, FALSE) )

repeats_by_line_bar_chart <- ggplot(data = num_repeats_by_line) +
  geom_col(aes(x = KO_line, y = Num_repeats), fill = 'steelblue3') +
  geom_text(aes(x = KO_line, y = Num_repeats, label = Num_repeats, colour = text_colour),
            angle = 45, nudge_y = 50, nudge_x = 0.25) +
  labs(x = 'Knockout Line', y = 'Number of repeat instances in DE\n(mutant vs sibling)') +
  scale_colour_manual(values = c('TRUE' = 'firebrick3', 'FALSE' = 'black'),
                      guide = 'none') +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = 'black', size = 10, face = 'italic',
                                   angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = 'black', size = 10),
        axis.title = element_text(size = 12),
        panel.grid = element_blank())

pdf(file = file.path('plots', 'num_repeats_by_line.pdf'))
print(repeats_by_line_bar_chart)
dev.off()

postscript(file = file.path('plots', 'num_repeats_by_line.eps'),
           width = 9, height = 5,
           horizontal = FALSE)
print(repeats_by_line_bar_chart)
dev.off()

## Fig S4b
# load repeats by type data
repeats_by_line_by_type <- read.delim(file = cmd_line_args$args[2])
# subset to sense and antisense
repeats_by_line_by_type_sense <- repeats_by_line_by_type[ repeats_by_line_by_type$Orientation == "sense", ]
repeats_by_line_by_type_antisense <- repeats_by_line_by_type[ repeats_by_line_by_type$Orientation == "antisense", ]
# and merge together
repeats_by_line_by_type <- merge(repeats_by_line_by_type_sense, repeats_by_line_by_type_antisense,
                                 by = c('KO_Line', 'Type'), all = TRUE)

repeats_by_line_by_type <- repeats_by_line_by_type[, c('KO_Line', 'Type', 'Count.x', 'Count.y')]
names(repeats_by_line_by_type) <- c('KO_Line', 'Type', 'Sense', 'Antisense')

# set NAs to zero
repeats_by_line_by_type$Sense[ is.na(repeats_by_line_by_type$Sense) ] <- 0
repeats_by_line_by_type$Antisense[ is.na(repeats_by_line_by_type$Antisense) ] <- 0

# subset to types DNA, LINE, LTR, SINE
repeats_by_line_by_type_subset <- do.call(rbind,
                                   lapply(c('DNA', 'LINE', 'LTR', 'SINE'),
                                          function(type){
                                            repeats_by_line_by_type[repeats_by_line_by_type$Type == type, ]
                                          }))
# drop levels
repeats_by_line_by_type_subset <- droplevels(repeats_by_line_by_type_subset)

# df for drawing a highlight box
inset_box <- data.frame(
    xmin = 0,
    xmax = 30,
    ymin = 0,
    ymax = 140
)

# set colour palette
repeat_type_palette <- c('#4e80bd', '#c04f4c', '#9bbb58', '#8164a2')
names(repeat_type_palette) <- levels(repeats_by_line_by_type_subset$Type)

repeats_by_line_by_type_plot <-
  ggplot(data = repeats_by_line_by_type_subset) +
  geom_rect(data = inset_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = 'grey95') +
  geom_hline(yintercept = 0, colour = "grey70", size = 0.5, linetype = 1) +
  geom_vline(xintercept = 0, colour = "grey70", size = 0.5, linetype = 1) +
  geom_point(aes(x = Antisense, y = Sense, colour = Type)) +
  scale_colour_manual(values = repeat_type_palette) +
  labs(x = 'DE repeats in introns (antisense)',
       y = 'DE repeats in introns (sense)') +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

pdf(file = file.path('plots', 'num_repeats_by_type.pdf'))
print(repeats_by_line_by_type_plot)
dev.off()

# Make inset with smaller limits to expand area close to origin
repeats_by_line_by_type_inset_plot <-
  ggplot(data = repeats_by_line_by_type_subset) +
  geom_hline(yintercept = 0, colour = "grey70", size = 0.5, linetype = 1) +
  geom_vline(xintercept = 0, colour = "grey70", size = 0.5, linetype = 1) +
  geom_point(aes(x = Antisense, y = Sense, colour = Type)) +
  scale_colour_manual(values = repeat_type_palette,
                      guide = 'none') +
  xlim(c(0,30)) + ylim(0,140) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_blank(),
        plot.background = element_rect(fill = 'grey95', colour = 'white') )

pdf(file = file.path('plots', 'num_repeats_by_type_inset.pdf'))
print(repeats_by_line_by_type_inset_plot)
dev.off()

# make composite figure
pdf(file = file.path('plots', 'num_repeats_by_type_combined.pdf'))
print(repeats_by_line_by_type_plot)
vp <- viewport(width = 0.5, height = 0.5, x = 0.9, y = 0.95,
               just = c('right', 'top'))
print(repeats_by_line_by_type_inset_plot, vp = vp)
dev.off()

postscript(file = file.path('plots', 'num_repeats_by_type_combined.eps'),
           width = 6, height = 5,
           horizontal = FALSE)
print(repeats_by_line_by_type_plot)
vp <- viewport(width = 0.5, height = 0.5, x = 0.85, y = 0.95,
               just = c('right', 'top'))
print(repeats_by_line_by_type_inset_plot, vp = vp)
dev.off()

## Fig S4d
# load padj/log2fc data
repeats_padj_log2fc <- read.delim(file = cmd_line_args$args[3])

# plot volcano plot
repeats_volcano_plot <- ggplot(data = repeats_padj_log2fc) +
    geom_point(aes(x = log2fc, y = -log10(padj), colour = Dhx35)) +
    scale_colour_manual(values = c('Dhx35_only' = '#4c81bd', 'Dhx35_plus' = '#c04e4b'),
                        name = '', labels = c('Dhx35_only' = 'Dhx35 only', 'Dhx35_plus' = 'Dhx35 plus \nother(s)')) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    labs(x = expression(log[2]*'['*Fold~Change*']'),
         y = expression('-'*log[10]*'['*Adjusted~pvalue*']') ) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 10))


pdf(file = file.path('plots', 'repeats-volcano.pdf'))
print(repeats_volcano_plot)
dev.off()

postscript(file = file.path('plots', 'repeats-volcano.eps'),
    width = 8, height = 4.5,
    horizontal = FALSE)
print(repeats_volcano_plot)
dev.off()

# Fig S4e
# lengths of repeat by type
repeats_de <- read.delim(file = cmd_line_args$args[4])

repeats_de$Type <- gsub('/.*$', '', repeats_de$Class)
repeats_de$Type[ repeats_de$Type != 'DNA' &
                repeats_de$Type != 'LINE' &
                repeats_de$Type != 'LTR' &
                repeats_de$Type != 'SINE' ] <- 'Other'
repeats_de$Type <- factor(repeats_de$Type,
                          levels = c('DNA', 'LINE', 'LTR', 'SINE', 'Other'))
# create family column
repeats_de$Family <- gsub(':.*$', '', repeats_de$Name)

# calculate length
repeats_de$Length <- repeats_de$End - repeats_de$Start + 1

# order by class and then length
repeats_de <- repeats_de[ order(repeats_de$Type, repeats_de$Family, -repeats_de$Length), ]

# set levels of Name to match ordering of data frame
repeats_de$Name <- factor(repeats_de$Name,
                          levels = unique(as.character(repeats_de$Name)))

# colour palette
type_palette <- colour_blind_palette[c('vermillion', 'blue_green', 'purple', 'blue', 'orange')]
names(type_palette) <- c('DNA', 'LINE', 'LTR', 'SINE', 'Other')

# calculate boxes to highlight L1MdGf_I, MMERGLN-int and MMETn-int
round_to <- 500
highlight_boxes <-
    do.call(rbind,
            lapply(c('L1MdGf_I', 'MMERGLN-int', 'MMETn-int'),
                function(family){
                    family_subset <- repeats_de[ repeats_de$Family == family, ]
                    data.frame(
                        Family = family, repeat_start = family_subset$Name[1],
                        repeat_end = family_subset$Name[nrow(family_subset)],
                        length_start = 0,
                        length_end = ceiling(family_subset$Length[1]/round_to)*round_to
                    )
                }
            )
    )


# plot of length
repeat_length_plot <- ggplot(data = repeats_de) +
    geom_point(aes(x = Name, y = Length, colour = Type)) +
    geom_rect(data = highlight_boxes,
              aes(xmin = repeat_start, xmax = repeat_end,
                  ymin = length_start, ymax = length_end),
              fill = 'grey90') + 
    geom_point(aes(x = Name, y = Length, colour = Type)) +
    geom_hline(yintercept = 4000, linetype = 'longdash', colour = 'firebrick3') +
    labs(x = 'Repeat Element', y = 'Repeat Length (bases)') +
    scale_y_continuous(expand = c(0.01,0) ) +
    scale_x_discrete(expand = c(0,15)) +
    scale_colour_manual(values = type_palette) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.5, linetype = 1,
                                   lineend = 'butt', arrow = FALSE),
          axis.title = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10)
        )

pdf(file = file.path('plots', 'repeats-length.pdf'))
print(repeat_length_plot)
dev.off()

postscript(file = file.path('plots', 'repeats-length.eps'),
           width = 8, height = 5,
           horizontal = FALSE)
print(repeat_length_plot)
dev.off()

# Fig S4f
repeats_all <- read.delim(file = cmd_line_args$args[5])

# subset to rows where at least one embryo is > 2 normalised counts
to_keep <- Reduce('|',
                    lapply(colnames(repeats_all)[grepl('normalised.count', colnames(repeats_all))],
                        function(col_name){
                            repeats_all[[col_name]] > 2
                        }
                    )
                )

repeats_subset <- repeats_all[to_keep, ]

# make family column
repeats_subset$Family <- gsub(':.*$', '', repeats_subset$Name)
# get normalised count columns plus family and melt
repeats_subset <-
    repeats_subset[ , c( colnames(repeats_subset)[ grepl('normalised.count', colnames(repeats_subset)) ],
                        'Family')]
names(repeats_subset) <- gsub('.normalised.count', '', names(repeats_subset))
rm(repeats_all)

repeats_subset_m <- melt(repeats_subset, id.vars = c('Family'),
                         variable.name = 'Sample', value.name = 'Normalised.Count')
# make column for genotype
repeats_subset_m$genotype <- gsub('Morc2a_', '', repeats_subset_m$Sample)
repeats_subset_m$genotype <- gsub('[0-9]+$', '', repeats_subset_m$genotype)

# make condition
repeats_subset_m$condition <- repeats_subset_m$genotype
repeats_subset_m$condition[ repeats_subset_m$condition == 'het' ] <- 'sib'
repeats_subset_m$condition[ repeats_subset_m$condition == 'wt' ] <- 'sib'

# sum up counts for each family for each sample
sample_sums_by_Family <- ddply(repeats_subset_m, .(Family, Sample, condition), summarise,
                              count_sum = sum(Normalised.Count)
                             )

mean_by_Family_by_condition <- ddply(sample_sums_by_Family, .(Family, condition), summarise,
                                       mean_counts = mean(count_sum)
                                )
# cast and calc mut/sib
mean_by_Family <- dcast(mean_by_Family_by_condition, Family ~ condition,
                        value.var = 'mean_counts')
mean_by_Family$log2fc <- log2(mean_by_Family$hom/mean_by_Family$sib)

# sum up counts for each family
total_by_Family <- ddply(repeats_subset_m, .(Family), summarise,
                              total_counts = sum(Normalised.Count)
                        )

plot_data <- merge(total_by_Family, mean_by_Family)

# get list of enriched family names to plot in a different colour
enriched_families <- read.delim(file = cmd_line_args$args[6])

plot_data$enriched <- rep('not_enriched', nrow(plot_data))
for(family in as.character(unique(enriched_families$Family))) {
    plot_data$enriched[ plot_data$Family == family ] <- 'enriched'
}
plot_data$enriched <- factor(plot_data$enriched,
                             levels = c('not_enriched', 'enriched'))

# remove Inf
plot_data <- plot_data[ is.finite(plot_data$log2fc), ]
counts_vs_fc_plot <-
    ggplot(data = plot_data,
            aes(x = log2fc, y = total_counts, colour = enriched)) +
    geom_hline(yintercept = 10^(0:6), colour = 'grey80') +
    geom_vline(xintercept = 0) +
    geom_point() +
    scale_colour_manual(values = c('enriched' = '#c04e4b', 'not_enriched' = '#4c81bd'),
                        name = '', labels = c('enriched' = 'Enriched', 'not_enriched' = 'Not Enriched')) +
    scale_y_log10(limits = c(1,1e6),
                  breaks = c(0,1,10,100,1000,10000,100000,1000000),
                  minor_breaks = NULL,
                  labels = function(breaks){ sprintf('%.0f', breaks) } ) +
    scale_x_continuous(limits = c(-4,8),
                       breaks = seq(-4,8,2)) +
    labs(x = expression(log[2]*'['*Fold~Change*']'),
         y = 'Total counts per repeat family') +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10)
         )


pdf(file = file.path('plots', 'repeats-counts_vs_fc.pdf'))
print(counts_vs_fc_plot)
dev.off()

postscript(file = file.path('plots', 'repeats-counts_vs_fc.eps'),
           width = 7, height = 5,
           horizontal = FALSE)
print(counts_vs_fc_plot)
dev.off()


#save_plot(file.path('plots', "repeats-supplfig.eps"),
#          plot_grid(repeats_by_line_bar_chart, repeats_by_line_by_type_plot,
#                    ggplot(), repeats_volcano_plot,
#                    ncol = 2, nrow = 3, rel_widths = c(6,9,3,6,9,3), align = 'h' ),
#          ncol = 6, device = 'eps',
#          base_height = 9, base_width = 1.5)

