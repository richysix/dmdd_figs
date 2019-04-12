library('optparse')

option_list <- list(
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'setup.R',
    usage = "Usage: %prog [options] samples_file baseline_rds_file" ),
  positional_arguments = 2
)

packages <- c('plyr', 'SummarizedExperiment')
for( package in packages ){
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

# setup directories
wd <- getwd()
for ( new_dir in c('plots', 'output') ) {
  dir_path <- file.path(wd, new_dir)
  if( !dir.exists(dir_path) ){
    dir.create(dir_path)
  }
}

# create KOs ordered by delay file
# load sample info
sample_file <- cmd_line_args$args[1]
sample_info <- read.table(file = sample_file, sep = "\t", header = TRUE, row.names = 1 )
sample_info$gene <- gsub('_[a-z0-9]+', '', row.names(sample_info))

# label with Theiler stage
stage_boundaries <- c(4, 7, 12, 19, 29, 34)
stage_labels <- c('TS12b', 'TS13', 'TS14', 'TS15', 'Other')

# assign stage numbers with their TS
sample_info$Theiler_stage <-
  cut(sample_info$somite_number,
      breaks = stage_boundaries, labels = stage_labels)
# 4 embryos are of Unknown stage. Put these in Other with the one TS16 (30s) embryo
sample_info$Theiler_stage[is.na(sample_info$Theiler_stage)] <- 'Other'

# sort KOs by delay
# count up numbers of embryos in each stage
stage_count <- ddply(sample_info, .(gene), summarise,
                     Severe = sum(Theiler_stage == 'TS12b'),
                     Moderate = sum(Theiler_stage == 'TS13'),
                     Slight = sum(Theiler_stage == 'TS14'),
                     None = sum(Theiler_stage == 'TS15')
                     )
# order by Severe, Moderate, Slight then None
stage_count <- 
  stage_count[order(stage_count$Severe, stage_count$Moderate, 
                    stage_count$Slight, stage_count$None, decreasing = TRUE), ]
# set levels
stage_count$gene <- factor(stage_count$gene, levels = stage_count$gene)

# get names of genes that are delayed
delayed <-
  stage_count$Severe + stage_count$Moderate + stage_count$Slight
delayed_genes <- stage_count$gene[ delayed > 0 ]
# write genes to file
write.table(delayed_genes, file = file.path(wd, 'output', 'KOs_delayed.txt'),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# output delayed order
delay <- factor(rep('None', nrow(stage_count)),
                levels = c('Severe', 'Moderate', 'Slight', 'None'))
for ( delay_category in c('Slight', 'Moderate', 'Severe')) {
  delay[ stage_count[[delay_category]] > 0 ] <- delay_category
}
write.table(data.frame(stage_count$gene, delay = delay),
            file = file.path(wd, 'output', 'KOs_ordered_by_delay.txt'),
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

# output sample info with gene and Theiler Stage
write.table(sample_info,
            file = file.path(wd, 'output', 'sample_info.txt'),
            quote = FALSE, sep = "\t",
            row.names = TRUE, col.names = NA)

# load baseline rds file and make all.tsv file and samples file
Mm_baseline <- readRDS(cmd_line_args$args[2])
counts <- assays(Mm_baseline)[['counts']]
colnames(counts) <- sub("$", " count", colnames(counts))
normalised_counts <- assays(Mm_baseline)[['normalised_counts']]
colnames(normalised_counts) <- sub("$", " normalised count", colnames(normalised_counts))

all_data <- data.frame(
  rowData(Mm_baseline),
  counts,
  normalised_counts,
  check.names = FALSE
)

write.table(all_data,
            file = file.path(wd, 'output', 'Mm_GRCm38_e88_baseline.tsv'),
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

sample_data <- colData(Mm_baseline)
sample_data <- sample_data[ , !grepl("sample_id", colnames(sample_data)) ]
write.table(sample_data,
            file = file.path(wd, 'output', 'samples-Mm_GRCm38_e88_baseline.txt'),
            quote = FALSE, sep = "\t",
            row.names = TRUE, col.names = NA)

file.create('setup.done')
