#!/usr/bin/env Rscript

suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(DESeq2))

Args        <- commandArgs()
countFile   <- ifelse(is.na(Args[6]),  "deseq2/counts.txt",             Args[6])
samplesFile <- ifelse(is.na(Args[7]),  "deseq2/samples.txt",            Args[7])
outputFile   <- ifelse(is.na(Args[8]),  "deseq2/normalised-counts.tsv",  Args[8])

# Get data and samples
countData <- read.delim(   countFile, header=TRUE, row.names=1, check.names = FALSE )
samples   <- read.delim( samplesFile, header=TRUE, row.names=1 )

counts <- countData[ , grepl("count", names(countData)) ]
names(counts) <- sub(" count", "", names(counts) )

# order counts by sample
counts <- counts[ , rownames(samples) ]

# Normalise
dds <- DESeqDataSetFromMatrix(counts, samples, design = ~ 1)
dds <- estimateSizeFactors(dds)

normalised_counts <- counts(dds, normalized=TRUE)
colnames(normalised_counts) <-
    sub("$", " normalised count", colnames(normalised_counts))

norm_counts <- data.frame(
    'Gene ID' = rownames(countData),
    countData[ , !grepl("count", names(countData)) ],
    normalised_counts,
    check.names = FALSE
)

# Write out normalised counts
write.table(norm_counts, file=outputFile, col.names=TRUE, quote=FALSE, sep="\t")
