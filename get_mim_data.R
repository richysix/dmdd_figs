#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'get_mim_data.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 1
)

if (cmd_line_args$options[['directory']] == 'cwd') {
  wd <- getwd()
} else {
  wd <- cmd_line_args$options[['directory']]
}

if( cmd_line_args$options[['verbose']] ){
  cat( "Working directory:", cmd_line_args$options[['directory']], "\n", sep=" " )
  cat( "Output files basename:", cmd_line_args$options[['basename']], "\n", sep=" " )
}

#cmd_line_args <- list(
#  options = list(directory = '/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd',
#                 verbose = FALSE ),
#  args = c('data/dmdd-genes.txt')
#)

packages <- c('biomaRt')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# get OMIM data from BioMart for all KOs
# input file contains the gene ids
gene_info_file <- cmd_line_args$args[1]
gene_info <- read.delim(gene_info_file)

# remove Cenpl
gene_info <- gene_info[ gene_info$Dir != 'Cenpl', ]

# use e88 version of biomart
Mm_ensembl88 <- useMart(host='mar2017.archive.ensembl.org', 
                        biomart='ENSEMBL_MART_ENSEMBL', 
                        dataset='mmusculus_gene_ensembl')

# get human homologues
hs_homologues <- getBM( filters = c('ensembl_gene_id'), 
                        values = gene_info[['Ensembl.ID']],
                        attributes = c('ensembl_gene_id',
                                    'hsapiens_homolog_ensembl_gene',
                                    'hsapiens_homolog_orthology_type',
                                    'hsapiens_homolog_orthology_confidence'),
                        mart = Mm_ensembl88)

# most are one-to-one but a few are one-to-many because of duplicate mouse genes
# ENSMUSG00000065586 has no orthologue. Mir96 so not surprising

# remove mouse genes that don't have a homologue
hs_homologues <-
  hs_homologues[ !is.na(hs_homologues$hsapiens_homolog_orthology_confidence), ]

# get MIM data for human homologues
Hs_ensembl88 <- useMart(host='mar2017.archive.ensembl.org', 
                        biomart='ENSEMBL_MART_ENSEMBL', 
                        dataset='hsapiens_gene_ensembl')

hs_mim_data <- getBM( filters = c('ensembl_gene_id'), 
                      values = hs_homologues$hsapiens_homolog_ensembl_gene,
                      attributes = c('ensembl_gene_id',
                                    'mim_morbid_accession',
                                    'mim_morbid_description'),
                      mart = Hs_ensembl88)

# merge all 3 data frame together
mouse_to_mim <- merge(hs_homologues, hs_mim_data,
                      by.x = 'hsapiens_homolog_ensembl_gene',
                      by.y = 'ensembl_gene_id')

mouse_to_human_mim <- merge(gene_info, mouse_to_mim, all.x = TRUE,
                            by.x = 'Ensembl.ID', by.y = 'ensembl_gene_id')

# write to file
write.table(mouse_to_human_mim, quote = FALSE, sep = "\t", 
            file = file.path('output', 'human-mim.tsv'),
            row.names = FALSE, col.names = TRUE)
