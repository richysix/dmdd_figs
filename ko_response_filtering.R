#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'ko_response_filtering.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 2
)

#cmd_line_args <- list(
#  options = list( verbose = FALSE ),
#  args = c('output/dir_number_to_gene.tsv',
#           'output/go_blind'
# )
#)

if( cmd_line_args$options[['verbose']] ){
}

packages <- c('GO.db', 'biovisr', 'miscr', 'Rgraphviz', 'plyr', 'ggplot2')
for (package in packages){
  library(package, character.only = TRUE)
}

key_file <- read.delim(cmd_line_args$args[1], header = FALSE)
colnames(key_file) <- c('Gene', 'Filtered', 'Dir')
# set order of levels of Filtered
key_file$Filtered <- factor(key_file$Filtered,
                            levels = c('unfiltered', 'filtered'))

base_dir <- cmd_line_args$arg[2]

# build GO graph
bp_terms <- select(GO.db, keytype = 'ONTOLOGY', keys = 'BP', columns = c('GOID', 'TERM'))
ontology_graph <- graphNEL( nodes = bp_terms$GOID,
                            edgemode='directed' )
bp_parents_list <- as.list(GOBPPARENTS)
bp_parents_list <- bp_parents_list[!is.na(bp_parents_list)]
num_edges <- sum(sapply(bp_parents_list, length)) - 1 # root node has an entry for parent called 'all'

from_nodes <- character(length = num_edges)
to_nodes <- character(length = num_edges)
i <- 1
for (from_node in bp_terms$GOID) {
#  cat(sprintf('%s\n', from_node))
  for (parent_node in bp_parents_list[[from_node]]) {
    # skip root node
    if (parent_node == 'all'){
      next
    }
    # get node index for parent
    from_nodes[i] <- from_node
    to_nodes[i] <- parent_node
    i <- i + 1
  }
}
  
# add edges to graph
ontology_graph <- addEdge(from_nodes, to_nodes, ontology_graph)

# go through each directory and get top terms
# do BP only for now
get_top_results_by_pvalue <- function(sig_file, num_results, gene,
                                      filtered_status){
  go_results <- read.delim(sig_file)
  # make sure that results are ordered by p.value
  go_results <- go_results[ order(go_results$pval), ]
  pval_threshold <- go_results$pval[num_results]
  go_results_subset <- go_results[ go_results$pval <= pval_threshold, ]
  results <- data.frame(
    Gene = gene,
    Filtered = filtered_status,
    GO.ID = go_results_subset$GO.ID,
    Term = go_results_subset$Term,
    pval = go_results_subset$pval,
    FoldEnrichment = go_results_subset$Significant / go_results_subset$Expected
  )
  return(results)
}

num_results <- 20
subset_results <-
  do.call(rbind,
          lapply(seq_len(nrow(key_file)),
                  function(i) {
                    return(
                      get_top_results_by_pvalue(
                        file.path(base_dir, key_file$Dir[i], 'BP.sig.tsv'),
                        num_results, key_file$Gene[i], key_file$Filtered[i])
                    )
                  }
         )
  )
# drop levels
subset_results <- droplevels(subset_results)
subset_results$log10pval <- -log10(subset_results$pval)

# get un-truncated terms using the GOTERMS from GO.db
full_terms <- sapply(as.character(subset_results[['GO.ID']]),
                      function(term){
                        desc <- as.character(subset_results$Term[ subset_results[['GO.ID']] == term ][1])
                        if (grepl("\\.\\.\\.$", desc)) {
                          # try and get untruncated version from GO.db
                          go_term <- GOTERM[[term]]
                          if (is.null(go_term)) {
                            if (debug) {
                              cat(term, "\n")
                            }
                            return(desc)
                          } else {
                            return(Term(go_term))
                          }
                        } else {
                          return(desc)
                        }
                      }
  )
subset_results[['Term']] <- full_terms
y_labels <- unique(subset_results$Term)
names(y_labels) <- unique(subset_results$GO.ID)

# make new factor of Gene and Filtered
subset_results$category <- interaction(subset_results$Filtered, subset_results$Gene)

# subset by Gene and plot by Filtered status
results_by_gene <- split(subset_results, subset_results$Gene)
results_by_gene <- lapply(results_by_gene, function(x){ droplevels(x) })
# sort by filtered status and then pvalue
results_by_gene <- lapply(results_by_gene,
                          function(results_df){
                            # make sure df is ordered by Filtered (filtered, then unfiltered)
                            # and then pvalue
                            results_df <- results_df[ order(as.character(results_df$Filtered), results_df$pval), ]
                            # set levels
                            # reverse levels so that filtered gets plotted at the top
                            results_df$GO.ID <- factor(results_df$GO.ID,
                                                       levels = rev(unique(results_df$GO.ID)))
                            return(results_df)
                          })

plot_list <- lapply(results_by_gene,
                    function(results_df){
                      labels <- unique(results_df$Term)
                      names(labels) <- unique(results_df$GO.ID)
                      bubble_plot(results_df, x = "Filtered", y = "GO.ID",
                                  size = "FoldEnrichment", fill = "log10pval",
                                  y_labels = labels) + labs(title = results_df$Gene[1])
                    }
                  )
pdf(file.path('plots', 'ko_response_filtering.pdf'),
    width = 12, height = 10)
invisible(lapply(plot_list, function(plot){ print(plot) }))
dev.off()

# redo for Fcho2 and Hira (not chosen yet)
# reopen results files and get numbers on where each of the top 20 terms is
# in the other list (if present in the other list at all)

# plot Fcho2 as an .eps file
genes <- c('Fcho2', 'Hira')
for (gene in genes) {
  dir_num <- key_file[ key_file$Gene == gene & key_file$Filtered == 'unfiltered', 'Dir']
  unfiltered_results_file <- file.path(base_dir, dir_num, 'BP.sig.tsv')
  unfiltered_results <- read.delim(unfiltered_results_file)
  
  dir_num <- key_file[ key_file$Gene == gene & key_file$Filtered == 'filtered', 'Dir']
  filtered_results_file <- file.path(base_dir, dir_num, 'BP.sig.tsv')
  filtered_results <- read.delim(filtered_results_file)
  
  results_df <- results_by_gene[[gene]]
  find_position_in_file <- function(term, results_to_search){
    posn <- which(as.character(results_to_search$GO.ID) == as.character(term))
    if (length(posn) == 0) {
      return(NA)
    } else if (length(posn) == 1) {
      return(posn[1])
    } else {
      # complain
      stop("Term appears more then once in results file\n",
           paste(term, unfiltered_results_file) )
    }
  }
  # go through term by term and get position in unfiltered set
  posns_in_unfilt <- sapply(results_df$GO.ID,
                            find_position_in_file, unfiltered_results)

  posns_in_filt <- sapply(results_df$GO.ID,
                          find_position_in_file, filtered_results)
  
  num_rows <- nrow(results_df)
  posn_df <- data.frame(
    GO.ID = factor(rep(results_df$GO.ID, 2),
                   levels = levels(results_df$GO.ID)),
    Filtered = factor(rep(c('unfiltered', 'filtered'), each = num_rows),
                       levels = c('unfiltered', 'filtered')),
    posns = c(posns_in_unfilt, posns_in_filt)
  )

  # output plot_data
  plot_data <- merge(results_df, posn_df, all.y = TRUE)
  plot_data <- unique(plot_data)
  plot_data$Gene <- rep(gene, nrow(plot_data))
  plot_data$Term <- sapply(as.character(plot_data$GO.ID),
         function(term){ return(plot_data$Term[ plot_data$GO.ID == term &
                                               !is.na(plot_data$Term) ][1] ) })
  write.table(plot_data,
    file = file.path('output', paste0('ko_response_filtering-', gene, '.tsv')),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

  labels <- unique(results_df$Term)
  names(labels) <- unique(results_df$GO.ID)
  plot <- bubble_plot(results_df, x = "Filtered", y = "GO.ID",
                      size = "FoldEnrichment", fill = "log10pval",
                      y_labels = labels, stroke = 0.2)
  plot <- plot +
    geom_text(data = posn_df, aes(x = Filtered, GO.ID, label = posns),
              nudge_x = -0.3)
              
  postscript(file = file.path('plots', paste0('ko_response_filtering-', gene, '.eps')),
             width = 7.5, height = 6.5, paper = 'special', horizontal = FALSE)
  print(plot + guides(size = 'none', fill = 'none'))
  dev.off()

  # plot a wider one so that the legend fits in
  postscript(file = file.path('plots', paste0('ko_response_filtering-', gene, '-legend.eps')),
             width = 11, height = 7, paper = 'special', horizontal = FALSE)
  print(plot)
  dev.off()
}


# collapse results based on GO graph
extracted_subtrees <- extract_subtrees(ontology_graph, as.character(subset_results$GO.ID))

# annotate each term with its parent
graph_levels <- topGO::buildLevels(ontology_graph)
root_node_lookup <- new.env(hash = T, parent = emptyenv())
for (subgraph in extracted_subtrees) {
  has_parents <- sapply(edges(subgraph),
                        function(edges){ length(edges) == 0 })
  root_node <- names(has_parents)[has_parents]
  
  for (node in nodes(subgraph)) {
    if (exists(node, envir = root_node_lookup)) {
      # check it's the same root node
      current_root_node <- root_node_lookup[[node]]
      if (length(current_root_node) == 1) {
        if (current_root_node != root_node) {
          if(graph_levels$nodes2level[[current_root_node]] ==
              graph_levels$nodes2level[[root_node]]) {
            # cat both terms together
            root_node_lookup[[node]] <- c(current_root_node, root_node)
          } else {
            # choose the term with the lowest level (highest number)
            # if new root node has higher number replace current root node
            # otherwise, leave as is
            if(graph_levels$nodes2level[[root_node]] >
                graph_levels$nodes2level[[current_root_node]]) {
              # change root node
              root_node_lookup[[node]] <- root_node
            }
          }
        }
      } else {
        # check whether root node is in the list of root nodes
        if ( !(root_node %in% current_root_node) ) {
          if (graph_levels$nodes2level[[current_root_node[1]]] ==
              graph_levels$nodes2level[[root_node]]) {
            # cat terms together
            root_node_lookup[[node]] <- c(current_root_node, root_node)
          } else {
            # choose the term with the lowest level (highest number)
            # if new root node has higher number replace current root node
            # otherwise, leave as is
            if(graph_levels$nodes2level[[root_node]] >
                graph_levels$nodes2level[[current_root_node[1]]]) {
              # change root node
              root_node_lookup[[node]] <- root_node
            }
          }
        }
      }
    } else {
      root_node_lookup[[node]] <- root_node
    }
  }
}

# add parent term to subset_results df
root_terms <- character(length = nrow(subset_results))
root_terms_description <- character(length = nrow(subset_results))
for (i in seq_len(nrow(subset_results))) {
  term_id <- as.character(subset_results$GO.ID)[i]
  root_term <- root_node_lookup[[term_id]]
  if (length(root_term) == 1) {
    root_terms[i] <- root_term
    root_terms_description[i] <- Term(GOTERM[[root_term]])
  } else {
    root_terms[i] <- paste0(root_term, collapse = '/')
    root_terms_description[i] <- paste0(sapply(root_term, function(x){ Term(GOTERM[[x]]) }), collapse = '/')
  }
}
subset_results$parent_term <- root_terms
subset_results$parent_term_description <- root_terms_description
# aggreagate terms to parent term
aggregated_results <-
  ddply(subset_results,
        .(Gene, Filtered, category, parent_term, parent_term_description),
        summarise,
        Count = length(GO.ID),
        Max_log10pval = max(log10pval),
        Max_FoldEnrichment = max(FoldEnrichment)
  )

aggregated_results$parent_term <- factor(aggregated_results$parent_term,
                                         levels = unique(aggregated_results$parent_term))
y_labels <- unique(aggregated_results$parent_term_description)
names(y_labels) <- unique(aggregated_results$parent_term)
pdf(file.path('plots', 'ko_response_filtering_aggregated.pdf'),
    width = 15, height = 50)
bubble_plot(aggregated_results, x = "category", y = "parent_term",
            size = "Count", fill = "Max_log10pval",
            y_labels = y_labels)
dev.off()

# terms don't aggregate particularly well
