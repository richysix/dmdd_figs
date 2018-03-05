.libPaths('./.R/lib')

# load packages
packages <- c('SummarizedExperiment')
for( package in packages ){
  suppressPackageStartupMessages(
    suppressWarnings( library(package, character.only = TRUE) )
  )
}

# load baseline data
baseline_count_file <-
    '/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/baseline-grandhet/all.tsv'
baseline_counts <- read.delim(baseline_count_file, check.names = FALSE)

baseline_sample_file <-
    '/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/baseline-grandhet/deseq2/samples.txt'
baseline_samples <- read.delim(baseline_sample_file, row.names = 1)
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

