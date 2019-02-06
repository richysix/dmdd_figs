.libPaths('./.R/lib')

library('optparse')

option_list <- list(
  make_option("--output_file", type="character", default='Mm_baseline_data.rda',
              help="Working directory [default %default]" ),
  make_option("--ensembl_version", type="integer", default=88,
              help="Working directory [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'baseline.R',
    usage = "Usage: %prog [options] countFile sampleFile\n",
    description = paste("Create SummarizedExperiment object for baseline data",
                   "from a count and sample file")
  ),
  positional_arguments = 2
)

#cmd_line_args <- list(
#  options = list(output_file = 'Mm_GRCm38_e90_baseline_data.rda',
#                 ensembl_version = 90,
#                 verbose = FALSE ),
#  args = c('/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/baseline-grandhet/baseline-90.tsv',
#           '/lustre/scratch117/maz/team31/projects/mouse_DMDD/mouse_dmdd_figs/baseline-samples.txt')
#)

# load packages
packages <- c('SummarizedExperiment')
for( package in packages ){
  suppressPackageStartupMessages(
    suppressWarnings( library(package, character.only = TRUE) )
  )
}

# load baseline data
baseline_counts <- read.delim(cmd_line_args$args[1], check.names = FALSE)

baseline_samples <- read.delim(cmd_line_args$args[2], row.names = 1)
# order levels of condition
baseline_samples$condition <-
    factor(baseline_samples$condition,
            levels = unique(baseline_samples$condition))

# label with Theiler stage
stage_boundaries <- c(1, 7, 12, 20, 29, 34, 39)
stage_labels <- c('TS12', 'TS13', 'TS14', 'TS15', 'TS16', 'TS17')

# assign stage numbers with their TS
# this issues a warning about introducing NAs by coercion. suppress it.
somites <- gsub('[A-Z]+_', '', baseline_samples$condition)
baseline_samples$somites <-
  as.integer( gsub('somites', '', somites) )

baseline_samples$Theiler_stage <-
  cut(baseline_samples$somites,
      breaks = stage_boundaries, labels = stage_labels)

# get normalised counts
counts <- 
    baseline_counts[, grepl(' count$', colnames(baseline_counts)) &
                      !grepl('normalised count$', colnames(baseline_counts))]
colnames(counts) <-
    gsub(' count$', '', colnames(counts))
rownames(counts) <- baseline_counts[, 'Gene ID']

norm_counts <-
    baseline_counts[, grepl('normalised count$', colnames(baseline_counts))]
colnames(norm_counts) <-
    gsub(' normalised count$', '', colnames(norm_counts))
rownames(norm_counts) <- baseline_counts[, 'Gene ID']

# row data. chr, start etc for each gene
baseline_row_data <- baseline_counts[, c('Gene ID', 'Chr', 'Start', 'End',
                                         'Strand', 'Biotype', 'Name',
                                         'Description')]

Mm_GRCm38_e88_baseline <-
    SummarizedExperiment(assays = list(counts = as.matrix(counts),
                                       norm_counts = as.matrix(norm_counts)),
                        rowData = DataFrame(baseline_row_data),
                        colData = DataFrame(baseline_samples))

# write out to rdata file
save(Mm_GRCm38_e88_baseline,
    file = file.path('data', 'Mm_GRCm38_e88_baseline.rda'),
    compress = "xz", compression_level = 9)

