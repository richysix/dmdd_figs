#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'process_emap.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 3
)

#cmd_line_args <- list(
#  options = list(directory = '/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd',
#                 verbose = FALSE ),
#  args = c('data/emap_results.all.tsv',
#           '/lustre/scratch117/maz/team31/projects/mouse_DMDD/emap/EMAPA_Nov_17_2017.obo',
#           'data/root_terms.txt')
#)

if (cmd_line_args$options[['directory']] == 'cwd') {
  wd <- getwd()
} else {
  wd <- cmd_line_args$options[['directory']]
}

if( cmd_line_args$options[['verbose']] ){
  cat( "Working directory:", cmd_line_args$options[['directory']], "\n", sep=" " )
}

plots_dir <- file.path(wd, 'plots')

packages <- c('ggplot2', 'viridis', 'reshape2', 'ontologyIndex', 'ontologyPlot',
              'plyr')
for( package in packages ){
    suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

# read in data for different sets
emap_results_file <- cmd_line_args$args[1]
emap_results <- read.delim(file = emap_results_file)
# split by results set
emap_results_list <- split(emap_results, emap_results$Set)

# read in EMAPA ontology
emap_obo_file <- cmd_line_args$args[2]
emap_ontology <- get_ontology(emap_obo_file,
                              propagate_relationships=c("is_a", "part_of"))
# make vector of ids named with term descriptions
ids_for <- emap_ontology$id
names(ids_for) <- emap_ontology$name

# collapse terms to specified levels
root_terms_file <- cmd_line_args$args[3]
root_terms <- read.table(root_terms_file, sep = "\t", header = FALSE,
                         col.names = c('Description'), stringsAsFactors = FALSE )
root_terms$Term.ID <- ids_for[root_terms$Description]
if(any(is.na(root_terms$Term.ID))) {
  # find first problem
  index <- which((is.na(root_terms$Term.ID)))
  stop("Couldn't match all terms to ids: ", paste(root_terms$Description[index], collapse = " "))
}

# go through each root term and get all descendants
# plot the graph as well
if( cmd_line_args$options[['verbose']] ){
  cat("Getting decendants for root terms.\n")
}
root_df_list <- vector('list', length = nrow(root_terms))
for(i in seq_len(nrow(root_terms))) {
  term <- root_terms$Term.ID[i]
  if( cmd_line_args$options[['verbose']] ){
    cat(root_terms$Description[i], "\n")
  }
  all_descendants <- get_descendants(emap_ontology, term)
  
  # make a df with the EMAPA ids, names and parental name/id
  num_rows <- length(all_descendants)
  root_df_list[[i]] <- data.frame(
    Term.ID = all_descendants,
    name = emap_ontology$name[all_descendants],
    parent_id = rep(term, num_rows),
    parent_name = rep(emap_ontology$name[term], num_rows)
  )
}
root_df <- do.call(rbind, root_df_list)
# save as rda file
save(root_df, file = file.path('output', 'root_terms.rda'))

# some terms appears as descendants of more than one root term
# output and make a descision on each one
emap_results <- do.call(rbind, emap_results_list)
emap_results_root_merge <- merge(emap_results, root_df)
term_counts <-
  table(unique(emap_results_root_merge[, c('Term.ID', 'parent_id')] )$Term.ID)
multiple_parents <-
  do.call(rbind,
          lapply(names(term_counts[ term_counts > 1 ]),
          function(term){
            unique(
              emap_results_root_merge[ emap_results_root_merge$Term.ID == term,
                                        c('Term.ID', 'Description',
                                          'parent_id', 'parent_name')])
          } ) )

write.table(multiple_parents,
            file = file.path('output', 'duplicated_terms.tsv'),
            quote = FALSE, sep = "\t", row.names = FALSE)

# save workspace
save.image(file = file.path(wd, 'output', 'emap.RData'))
