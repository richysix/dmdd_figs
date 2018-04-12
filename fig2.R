#!/usr/bin/env Rscript

library('optparse')

option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option("--debug", action="store_true", default=FALSE,
              help="Working directory [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'fig2.R',
    usage = "Usage: %prog [options] input_file" ),
  positional_arguments = 6
)

#cmd_line_args <- list(
#  options = list(directory = 'cwd',
#                 debug = TRUE,
#                 verbose = FALSE ),
#  args = c('data/go_results.tsv',
#           'output/KOs_ordered_by_delay.txt',
#           'data/emap_results.all.tsv',
#           'output/duplicated_terms-edited.tsv',
#           'data/sig_gene_counts.tsv',
#           'output/mrna_abnormal-jaccard-all.rda'
# )
#)

debug <- cmd_line_args$options[['debug']]
verbose <- cmd_line_args$options[['verbose']]

if( verbose ){
  cat( "Working directory:",
      cmd_line_args$options[['directory']], "\n", sep=" " )
}

plots_dir <- file.path(getwd(), 'plots')

packages <- c('ggplot2', 'viridis', 'reshape2', 'ontologyIndex', 'ontologyPlot',
              'plyr', 'svglite', 'GO.db', 'devtools', 'cowplot', 'ggdendro',
              'ggrepel', 'grid')
for( package in packages ){
  library(package, character.only = TRUE)
}
# load biovisr and miscr. If not installed, install from github repo
packages <- c('biovisr', 'miscr')
for( package in packages ){
  library( package, character.only = TRUE )
}

# GO enrichment summary
go_results_file <- cmd_line_args$args[1]
go_results <- read.delim(go_results_file)

# read in KO ordering
ko_order_file <- cmd_line_args$args[2]
ko_order <- read.table(ko_order_file)
names(ko_order) <- c('Gene', 'Delay Category')

# set levels of genes
go_results$Gene <- factor(go_results$Gene,
                          levels = ko_order[[1]])

# get un-truncated terms using the GOTERMS from GO.db
full_terms <- sapply(as.character(go_results[['GO.ID']]),
                      function(term){
                        desc <- as.character(go_results$Term[ go_results[['GO.ID']] == term ][1])
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
go_results[['Term']] <- full_terms

go_results_by_domain <- split(go_results, go_results$Domain)

# FUNCTION
# top_shared_terms - get the specified number of terms that appear the most times
top_shared_terms <- function(go_results, num_terms) {
  # get counts of number of KOs a term appears in
  go_term_counts <- table(go_results$GO.ID)
  go_term_counts <- go_term_counts[ order(go_term_counts, decreasing = TRUE) ]
  # get top most shared terms
  top_terms <- names(go_term_counts[ go_term_counts >= go_term_counts[num_terms]])
  top_terms_by_shared <-
    do.call(rbind,
            lapply( top_terms,
                    function(term){
                      go_results[ go_results$GO.ID == term, ]
                    }
            )
          )
  
  # order levels of GO.ID by number mutants
  # use reverse order because first level is plotted at the bottom (i.e. y = 1)
  top_terms_by_shared$GO.ID <-
    factor(top_terms_by_shared$GO.ID, levels = rev(top_terms) )
  ## return plot and data table
  #return( list(data = top_shared_terms_table,
  #             plot = top_shared_go_results_plot) )
  return(top_terms_by_shared)
}

# calculate overlap by terms
terms_overlap <- function(results, mut1, mut2, term_col_name) {
  jaccard_index =
    length( intersect(results[ results$Gene == mut1, term_col_name ],
                            results[ results$Gene == mut2, term_col_name ]) ) /
    length( union(results[ results$Gene == mut2, term_col_name ],
                  results[ results$Gene == mut1, term_col_name ]) )
  return( jaccard_index )
}

# use hierarchical clustering to order by the overlap in terms
cluster_by_jaccard <- function(results, term_col_name){
  # cluster genes by jaccard distance
  genes <- unique(results$Gene)
  jaccard_index <- numeric(length = length(genes))
  i <- 1
  for ( mut1 in genes) {
    for ( mut2 in genes) {
      jaccard_index[i] <- terms_overlap(results, mut1, mut2, term_col_name)
      i <- i + 1
    }
  }
  
  jaccard_index <- matrix(jaccard_index, nrow = length(genes),
                              dimnames = list(genes, genes) )
  # cluster
  jaccard_dist <- as.dist(1 - jaccard_index)
  clust_by_overlap <- hclust(jaccard_dist, method = "ward.D2")
  
  # reorder matrix by clustering
  jaccard_index <- jaccard_index[ clust_by_overlap$order,
                                          clust_by_overlap$order ]
  
  # plot reordered matrix as heatmap
  jaccard_index_m <- melt(jaccard_index)
  # reverse levels of y axis variable so that it plots sensibly
  jaccard_index_m$Var2 <- factor( jaccard_index_m$Var2,
                              levels = rev(levels(jaccard_index_m$Var2))
    )
  jaccard_index_heatmap <- ggplot(data = jaccard_index_m) + 
    geom_raster( aes( x = Var1, y = Var2, fill = value ) ) +
    scale_x_discrete( position = 'top') +
    scale_fill_viridis(name = "Jaccard\nIndex") + 
    theme_minimal() +
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title = element_blank() )
  
  # return hclust object, jaccard matrix, heatmap object
  return(
    list(
      jaccard_matrix = jaccard_index,
      clust_obj = clust_by_overlap,
      heatmap = jaccard_index_heatmap
    )
  )
}

# input is ko_response only
results_set <- 'ko_response'
for ( domain in names(go_results_by_domain) ) {
  go_results_subset <- go_results_by_domain[[domain]]
  
  cluster_results <- cluster_by_jaccard(go_results_subset, 'GO.ID')
  
  # plot the tree as a dendrogram
  tree_file <-
    file.path(plots_dir,
              paste(results_set, domain, 'mutants-clust-by-go',
                    'tree', 'pdf', sep=".") )
  pdf(tree_file)
  plot(cluster_results$clust_obj)
  dev.off()
  
  # output reordered matrix to file
  matrix_file <-
    file.path('output',
              paste(results_set, domain,
                    'mutants-clust-by-go', 'tsv', sep=".") )
  write.table(cluster_results$jaccard_matrix, file=matrix_file, quote=FALSE,
              row.names=TRUE, col.names=NA, sep="\t")
  
  mat_plot_file <-
    file.path(plots_dir,
              paste(results_set, domain,
                    'mutants-clust-by-go', 'png', sep=".") )
  png(filename = mat_plot_file,
       width = 480, height = 480, units = "px", pointsize = 12,
       bg = "white")
  print(cluster_results$heatmap)
  dev.off()

  # reorder levels of Gene by tree and merge to root terms
  go_results_subset$Gene <-
    factor(go_results_subset$Gene,
           levels = row.names(cluster_results$jaccard_matrix) )
  go_results_by_domain[[domain]] <- go_results_subset
  
  top_terms_by_shared <- top_shared_terms(go_results_subset, num_terms = 50)
  go_term_counts <- table(top_terms_by_shared$GO.ID)
  
  # get descriptions in the correct order
  term_labels <-
    sapply(levels(top_terms_by_shared$GO.ID),
          function(id){
            top_terms_by_shared$Term[ top_terms_by_shared[['GO.ID']] == id ][1]
          })
  
  top_shared_go_results_plot <-
    bubble_plot(top_terms_by_shared,
                x = 'Gene', y = 'GO.ID', size = 'Fold.Enrichment',
                fill = 'X.log10.pvalue.', x_labels = NULL,
                y_labels = term_labels )
    
  postscript(file = file.path(plots_dir,
                    paste('top_shared_go_results', domain, 'eps', sep = ".")),
              width = 8, height = 10, paper = 'special')
  print(top_shared_go_results_plot)
  dev.off()
  
  postscript(file = file.path(plots_dir,
                    paste('top_shared_go_results-no-y-text', domain, 'eps', sep = ".")),
              width = 5.6, height = 10, paper = 'special')
  print(top_shared_go_results_plot +
          theme(axis.text.y = element_blank(),
                axis.text.x = element_text(colour = 'black'),
                panel.grid.major =
                  ggplot2::element_line(colour = 'grey80', linetype = 'dashed', size = 0.2))
  )
  dev.off()
  
  # make table to reflect the plot
  top_shared_terms_table <- dcast(top_terms_by_shared, GO.ID + Term ~ Gene,
                                  value.var = 'X.log10.pvalue.')
  top_shared_terms_table <- merge(data.frame( GO.ID = names(go_term_counts),
                                             Count = as.integer(go_term_counts) ),
                                  top_shared_terms_table)
  # sort by Occurence
  top_shared_terms_table <-
    top_shared_terms_table[ order(top_shared_terms_table$Count, decreasing = TRUE), ]
  
  # output plot data as table
  # write to file
  write.table(top_shared_terms_table,
              file = file.path(wd, 'output',
                               paste('top_50_go_p', domain, 'tsv', sep = ".")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

################################################################################
## get top 50 by pvalue
#top_hits_by_pvalue <- function(terms_df, num = 50,
#                               pvalue_col = 'pval',
#                               term_col = 'GO.ID') {
#  top_by_pvalue <- terms_df[ order(terms_df[[pvalue_col]]), ]
#  top_GO.ID_by_pvalue <- as.character(unique(top_by_pvalue[[term_col]])[seq_len(num)])
#  top_terms_by_pvalue <-
#    do.call(rbind,
#            lapply( top_GO.ID_by_pvalue,
#                    function(term, terms_df){
#                      terms_df[ terms_df[[term_col]] == term, ]
#                    },
#                    terms_df
#            )
#          )
#  return(top_terms_by_pvalue)
#}
#
#
#top_terms_by_pvalue <- top_hits_by_pvalue(go_results, num = 50)
#
#top_terms_by_pvalue$GO.ID <-
#  factor(top_terms_by_pvalue$GO.ID, levels = rev(top_GO.ID_by_pvalue) )
## get descriptions in the correct order
#term_labels <-
#  sapply(levels(top_terms_by_pvalue$GO.ID),
#         function(id){ term_description_hash[[id]] } )
#
#top_pvalue_go_results_plot <- GO_bubble_plot(top_terms_by_pvalue, term_labels)
#
#pdf(file = 'plots/top_pvalue_go_results.pdf', width = 10, height = 10)
#print(top_pvalue_go_results_plot)
#dev.off()
#
#################################################################################
## and top 50 by FE
top_hits_fe <- function(go_results, num_terms = 50){
  top_by_FE <-
    go_results[ order(go_results$Fold.Enrichment, decreasing = TRUE), ]
  top_GO.ID_by_FE <- as.character(unique(top_by_FE$GO.ID)[seq_len(num_terms)])
  top_terms_by_FE <-
    do.call(rbind,
            lapply( top_GO.ID_by_FE,
                    function(term){
                      go_results[ go_results$GO.ID == term, ]
                    }
            )
          )
  
  top_terms_by_FE$GO.ID <-
    factor(top_terms_by_FE$GO.ID, levels = rev(top_GO.ID_by_FE) )
  term_labels <-
    sapply(levels(top_terms_by_FE$GO.ID),
           function(id){ top_terms_by_FE$Term[ top_terms_by_FE[['GO.ID']] == id ][1] } )
  
  top_fe_go_results_plot <-
    bubble_plot(top_terms_by_FE,
                x = 'Gene', y = 'GO.ID', size = 'Fold.Enrichment',
                fill = 'X.log10.pvalue.', x_labels = NULL,
                y_labels = term_labels )
  
  # make table to reflect the plot
  top_fe_terms_table <- dcast(top_terms_by_FE, GO.ID + Term ~ Gene,
                                  value.var = 'X.log10.pvalue.')
  
  # return plot and data table
  return( list(data = top_fe_terms_table,
               plot = top_fe_go_results_plot) )
}

# get top terms by fold enrichment
# do all sets
for ( domain in names(go_results_by_domain) ) {
  top_results <- top_hits_fe(go_results_by_domain[[domain]],
                              num_terms = 50)
  
  #pdf(file = file.path(plots_dir,
  #                  paste('top_FE_go_results', domain, 'pdf', sep = ".")),
  #    width = 10, height = 10)
  png(file = file.path(plots_dir,
                    paste('top_FE_go_results', domain, 'png', sep = ".")),
      width = 800, height = 800)
  print(top_results[['plot']] +
        scale_fill_viridis(direction = -1,
                           guide = guide_colourbar(title = expression(paste('-', log[10], "[pvalue]", sep = '')))) +
        scale_size(name = 'Fold Enrichment')
        )
  dev.off()
  
  # output plot data as table
  # write to file
  write.table(top_results[['data']],
              file = file.path(wd, 'output',
                               paste('top_50_go_fe', domain, 'tsv', sep = ".")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

#################################################################################


filtering_data <- read.delim('data/go_results-filtering.tsv', header = TRUE)
# set levels of Filtered
filtering_data$Filtered <- factor(filtering_data$Filtered,
                                  levels = c('unfiltered', 'filtered'))

filtering_counts <-
    ddply(filtering_data,
          .(Gene, Comparison, Domain, Filtered), summarise,
          num_sig_terms = length(GO.ID)
         )

filtering_counts_by_gene_by_domain <-
  split(filtering_counts,
        f = list(filtering_counts$Gene, filtering_counts$Domain))

# If a mutant has no sig terms it won't appear in filtering_counts
# Add a row with num_sig_terms = 0 for each of these
filtering_counts_edited <-
  do.call(rbind,
    lapply(filtering_counts_by_gene_by_domain,
            function(subset_by_gene_by_domain){
              if (nrow(subset_by_gene_by_domain) < 2) {
                by_gene_by_domain_list <- vector('list', length = 2)
                i <- 1
                for (filter_class in c("unfiltered", "filtered")) {
                  if (any(subset_by_gene_by_domain$Filtered == filter_class)) {
                    by_gene_by_domain_list[[i]] <- subset_by_gene_by_domain[ subset_by_gene_by_domain$Filtered == filter_class, ]
                  } else {
                    by_gene_by_domain_list[[i]] <-
                      data.frame(
                        Gene = subset_by_gene_by_domain$Gene[1],
                        Comparison = subset_by_gene_by_domain$Comparison[1],
                        Domain = subset_by_gene_by_domain$Domain[1],
                        Filtered = filter_class,
                        num_sig_terms = 0
                      )
                  }
                  i <- i + 1
                }
                return(do.call(rbind, by_gene_by_domain_list))
              } else {
                return(subset_by_gene_by_domain)
              }
            }
          )
        )

filtering_counts_w <- dcast(filtering_counts_edited, Gene + Domain ~ Filtered,
                            value.var = 'num_sig_terms')

# Plot Difference in sig genes between Filtered and unfiltered
# read in sig genes file
sig_genes_file <- cmd_line_args$args[5]
sig_genes <- read.delim(file = sig_genes_file)

# get sig_genes data for delayed mutants
sig_genes_delayed <- 
  do.call(rbind,
          lapply(levels(filtering_counts_w$Gene),
                 function(gene_name){
                   return(sig_genes[ sig_genes$Gene == gene_name, ])
                 })
        )
sig_genes_delayed <- sig_genes_delayed[ sig_genes_delayed$Set == 'ko_response', ]

sig_genes_delayed_subset <-
  do.call(rbind,
    lapply(unique(as.character(sig_genes_delayed$Gene)),
              function(gene_name, sig_genes_delayed) {
                gene_subset <-
                  sig_genes_delayed[ sig_genes_delayed$Gene == gene_name, ]
                for (comp in c('hom_vs_het_wt', 'hom_vs_het',
                                'het_vs_wt')) {
                  sig_genes_delayed_count <-
                    gene_subset[ gene_subset$Comparison == comp, ]
                  if(nrow(sig_genes_delayed_count) == 2) {
                    return(sig_genes_delayed_count)
                  }
                }
              },
              sig_genes_delayed
    )
  )

sig_genes_delayed_subset$Type <- factor(sig_genes_delayed_subset$Type,
                                        levels = c('unfiltered', 'filtered'))
sig_genes_delayed_subset_w <- dcast(sig_genes_delayed_subset, Gene ~ Type,
                            value.var = 'Count')

filtering_genes_plot <- ggplot() +
  geom_point(data = sig_genes_delayed_subset,
             aes(x = Type, y = Count,
                 colour = Gene)) +
  geom_segment(data = sig_genes_delayed_subset_w,
               aes(y = unfiltered, yend = filtered, colour = Gene,
                   x = 'unfiltered', xend = 'filtered')) +
  geom_text_repel(data = sig_genes_delayed_subset[ sig_genes_delayed_subset$Type == 'unfiltered', ],
                  aes(x = Type, y = Count, label = Gene),
                  xlim = c(NA, 1) )

# print to pdf
pdf(file = file.path('plots', 'filtering.pdf'))
print(filtering_genes_plot)
print(filtering_genes_plot + scale_y_log10())

# subset each to BP
for (domain in levels(filtering_counts_w$Domain)) {
  filtering_counts_subset <- filtering_counts[ filtering_counts$Domain == domain, ]
  filtering_counts_w_subset <- filtering_counts_w[ filtering_counts_w$Domain == domain, ]
  
  filtering_plot <- ggplot() +
    geom_point(data = filtering_counts_subset,
               aes(x = Filtered, y = num_sig_terms,
                   colour = Gene)) +
    geom_segment(data = filtering_counts_w_subset,
                 aes(y = unfiltered, yend = filtered, colour = Gene,
                     x = 'unfiltered', xend = 'filtered')) +
    geom_text_repel(data = filtering_counts_subset[ filtering_counts_subset$Filtered == 'unfiltered', ],
                    aes(x = Filtered, y = num_sig_terms, label = Gene),
                    xlim = c(NA, 1) ) +
    labs(title = paste0('Effect of Filtering on the number of GO enrichments (domain = ', domain, ')')) +
    theme(axis.title.x = element_blank())
  
  print(filtering_plot)
  print(filtering_plot + scale_y_log10())

  # plot sig genes against sig terms
  merged_data <- merge(sig_genes_delayed_subset, filtering_counts_subset,
           by.x = c('Gene', 'Comparison', 'Type'),
           by.y = c('Gene', 'Comparison', 'Filtered'))
  scatter_plot <- ggplot(data = merged_data) +
    geom_point(aes(x = Count, y = num_sig_terms, colour = Gene, shape = Type)) +
    geom_smooth(aes(x = Count, y = num_sig_terms),
                formula = y ~ x, method = 'lm') + 
    scale_x_log10() + scale_y_log10()
  print(scatter_plot)
}
dev.off()


################################################################################
# EMAPA
emap_results_file <- cmd_line_args$args[3]
emap_results <- read.delim(file = emap_results_file)
# order genes by delay order
emap_results$Gene <- factor(emap_results$Gene,
                            levels = ko_order[[1]])
emap_results$Gene <- droplevels(emap_results$Gene)

# split by results set
emap_results_list <- split(emap_results, emap_results$Set)

# LOAD root_df
load(file.path('output', 'root_terms.rda'))

# read in edited file
edited_file <- cmd_line_args$args[4]
duplicate_terms <- read.delim(file = edited_file)

duplicates_to_remove <-
  duplicate_terms[ duplicate_terms$to_use == 0, c('Term.ID', 'parent_id') ]

to_remove <-
  Reduce('|', lapply(row.names(duplicates_to_remove),
       function(x, root_df){
        term1 <- as.character(duplicates_to_remove[x, 'Term.ID'])
        term2 <- as.character(duplicates_to_remove[x, 'parent_id'])
        to_remove <- root_df$Term.ID == term1 & root_df$parent_id == term2
        return(to_remove)
      },
      root_df)
  )
root_df <- root_df[ !to_remove, ]

# calculate overlap by terms
terms_overlap <- function(results, mut1, mut2){
  jaccard_index =
    length( intersect(results$Term.ID[ results$Gene == mut1 ],
                            results$Term.ID[ results$Gene == mut2 ]) ) /
    length( union(results$Term.ID[ results$Gene == mut2 ],
                  results$Term.ID[ results$Gene == mut1 ]) )
  return( jaccard_index )
}

# go through each results set and plot a bubble plot
plot_list <- vector('list', length = length(emap_results_list))
names(plot_list) <- names(emap_results_list)
clust_objs <- vector('list', length = length(emap_results_list))
names(clust_objs) <- names(emap_results_list)

# calculate maximum log10(pvalue)
max_log10_pvalue <- 0
max_pvalue_set <- ''
for (set in names(emap_results_list)) {
  current_max <- max( emap_results_list[[set]][['X.log10.pvalue.']] )
  if (current_max > max_log10_pvalue) {
    max_log10_pvalue <- current_max
    max_pvalue_set <- set
  }
}

for (results_set in names(emap_results_list)) {
  if (debug) {
    cat(results_set, "\n")
  }
  results <- emap_results_list[[results_set]]
  # reorder levels of Gene by tree and merge to root terms
  results$Gene <- factor(results$Gene,
                         levels = unique(results$Gene) )
  results_summarise <- merge(results, root_df, all.x = TRUE)
  #dim(unique(results_summarise[ is.na(results_summarise$name),
  #                             c('Term.ID', 'Description')]))
  # remove NAs
  results_summarise_filtered <-
      results_summarise[ !is.na(results_summarise$name), ]
  
  # aggregate terms that have the same parent id
  results_aggregated <-
    ddply(results_summarise_filtered,
          .(Gene, parent_id, parent_name), summarise,
          num_sig_terms = length(Term.ID),
          min_p_value = min(Adjusted.p.value),
          max_log10_p = max(X.log10.pvalue.)
         )
  # droplevels
  results_aggregated$Gene <- droplevels(results_aggregated$Gene)
  # cluster genes by jaccard distance
  genes <- unique(results_aggregated$Gene)
  jaccard_index <- numeric(length = length(genes))
  i <- 1
  for ( mut1 in genes) {
    for ( mut2 in genes) {
      jaccard_index[i] <- terms_overlap(results, mut1, mut2)
      i <- i + 1
    }
  }
  
  jaccard_index <- matrix(jaccard_index, nrow = length(genes),
                              dimnames = list(genes, genes) )
  # cluster
  jaccard_dist <- as.dist(1 - jaccard_index)
  clust_by_overlap <- hclust(jaccard_dist, method = "ward.D2")
  clust_objs[[results_set]] <- clust_by_overlap
  
  # reorder levels of Gene by clustering
  results_aggregated$Gene <- factor(results_aggregated$Gene,
                                    levels = levels(results_aggregated$Gene)[clust_by_overlap$order] )

  # plot a dendrogram
  tree_file <-
    file.path(plots_dir,
              paste(results_set, 'mutants-clust-tree', 'pdf', sep=".") )
  pdf(tree_file)
  plot(clust_by_overlap)
  dev.off()
  
  # reorder matrix by clustering
  jaccard_index <- jaccard_index[ clust_by_overlap$order,
                                          clust_by_overlap$order ]
  # output reordered matrix to file
  matrix_file <-
    file.path('output',
              paste(results_set, 'mutants-clust-mat', 'tsv', sep=".") )
  write.table(jaccard_index, file=matrix_file, quote=FALSE,
              row.names=TRUE, col.names=NA, sep="\t")
  
  # plot reordered matrix as heatmap
  jaccard_index_m <- melt(jaccard_index)
  # reverse levels of y axis variable so that it plots sensibly
  jaccard_index_m$Var2 <- factor( jaccard_index_m$Var2,
                              levels = rev(levels(jaccard_index_m$Var2))
    )
  jaccard_index_heatmap <- ggplot(data = jaccard_index_m) + 
    geom_raster( aes( x = Var1, y = Var2, fill = value ) ) +
    scale_x_discrete( position = 'top') +
    scale_fill_viridis(name = "Jaccard\nIndex") + 
    theme_minimal() +
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title = element_blank() )
  
  mat_plot_file <-
    file.path(plots_dir,
              paste(results_set, 'mutants-clust_by_emapa', 'png', sep=".") )
  png(filename = mat_plot_file,
       width = 480, height = 480, units = "px", pointsize = 12,
       bg = "white")
  print(jaccard_index_heatmap)
  dev.off()
  
  # reverse levels for plotting
  results_aggregated$parent_id <-
    factor(results_aggregated$parent_id,
           levels = rev(as.character(levels(results_aggregated$parent_id))))
  results_aggregated$parent_name <-
    factor(results_aggregated$parent_name,
           levels = rev(as.character(levels(results_aggregated$parent_name))))
  
  emapa_bubble_plot <-
    bubble_plot(results_aggregated, x = 'Gene',
      y = 'parent_id', size = 'num_sig_terms',
      fill = 'max_log10_p', y_labels = levels(results_aggregated$parent_name),
      stroke = 0.2) +
      scale_fill_viridis(limits = c(1,max_log10_pvalue),
                          direction = -1) +
      theme(legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.position = 'bottom', legend.direction = 'horizontal',
            legend.box = 'vertical' )
  
  plot_list[[results_set]] <- emapa_bubble_plot
  
  postscript(file = file.path(plots_dir,
                       paste0(results_set, '.emap_aggregated_summary.eps')),
      width = 7.5, height = 8)
  print(emapa_bubble_plot)
  dev.off()
}

# get legend of max p set
overall_legend <- get_gg_legend(plot_list[[max_pvalue_set]])

plot_list <- lapply(plot_list,
                    function(plot){
                      plot + 
                        labs(fill = expression(paste('Max(-', log[10], '[pvalue])', sep = '')),
                             size = 'Number of Significant Terms') +
                        theme_void() +
                        theme(
                          axis.text.x = element_text(size = 12, colour = 'black', angle = 90,
                                                              hjust = 0, debug = FALSE),
                          axis.text.y = element_text(size = 12, colour = 'black', angle = 0,
                                                              hjust = 1, debug = FALSE),
                          panel.grid.major = element_line(colour = 'grey80', linetype = 'dotted'),
                          legend.title = element_text(size = 14),
                          legend.text = element_text(size = 12),
                          legend.position = 'none')
                    })

# plot
save_plot(file.path(plots_dir, paste0("emap-bubble_plots-by_results_set.eps")),
          plot_grid(plotlist = plot_list,
          nrow = 2, ncol = 2, align = 'vh', axis = 'tblr'),
          nrow = 2, ncol = 2, device = 'eps',
          base_width = 7.5, base_height = 6)
postscript(file = file.path(plots_dir, paste0("emap-bubble_plots-by_results_set-legend.eps")),
           width = 7.5, height = 2, paper = 'special', horizontal = FALSE)
grid.draw(overall_legend)
dev.off()

# Gene list overlaps
# plot tree and overlap matrix for mrna_abnormal
results_set <- 'mrna_abnormal'
  
# load cluster and heatmap plot and output both together with annotations
load(cmd_line_args$args[6])
tree_plot_data <- dendro_data(jaccard_clust)
delay_plot_data <- label(tree_plot_data)
delay_plot_data$y <- rep(-0.1, nrow(delay_plot_data))

# add annotations for delay category
delay_cat <- ko_order$`Delay Category`
names(delay_cat) <- ko_order$Gene
delay_plot_data$Delay <- factor(delay_cat[as.character(delay_plot_data$label)],
                                levels = c('None', 'Slight', 'Moderate', 'Severe'))

colour_palette <- 
   c( 'blue' = rgb(0,0.45,0.7),
      'yellow' = rgb(0.95, 0.9, 0.25),
      'sky_blue' = rgb(0.35, 0.7, 0.9),
      'purple' = rgb(0.8, 0.6, 0.7)
)
names(colour_palette) <- levels(delay_plot_data$Delay)

# get sig genes data
sig_genes_subset <- sig_genes[ sig_genes$Set == results_set &
                                sig_genes$Type == 'filtered', ]
sig_genes_subset <-
  do.call(rbind,
    lapply(unique(as.character(sig_genes_subset$Gene)),
              function(gene_name, sig_genes_subset) {
                gene_subset <-
                  sig_genes_subset[ sig_genes_subset$Gene == gene_name, ]
                for (comp in c('hom_vs_het_wt', 'hom_vs_het',
                                'het_vs_wt')) {
                  sig_genes_count <-
                    gene_subset[ gene_subset$Comparison == comp, ]
                  if(nrow(sig_genes_count) == 1) {
                    return(sig_genes_count)
                  }
                }
              },
              sig_genes_subset
    )
  )
# check that gene names are unique
if (length(unique(sig_genes_subset$Gene)) != nrow(sig_genes_subset)) {
  stop(paste0('mrna_abnormal gene-list overlap tree and matrix plot: ',
              'Gene names in significant genes data set are not unique' ))
}
row.names(sig_genes_subset) <- sig_genes_subset$Gene
sig_genes_subset <- sig_genes_subset[as.character(delay_plot_data$label), ]
sig_genes_subset$x <- label(tree_plot_data)$x
sig_genes_subset$y <- rep(0, nrow(sig_genes_subset))

dendro_plot <- ggplot() +
  geom_point(data = delay_plot_data, aes(x = label, y = y, fill = Delay),
             size = 4, shape = 23) +
  geom_segment(data = segment(tree_plot_data), size = 0.3, lineend = 'square',
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = sig_genes_subset, aes(x = Gene, y = y, label = Count),
            size = 4.2, nudge_x = 0.4, nudge_y = 0.3,
            angle = 90, hjust = 1) +
  scale_fill_manual(values = colour_palette) +
  theme_void() + 
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'top')

# plot both tree and bubble plot together
save_plot(file.path(plots_dir, paste0(results_set, "-gene_list-overlaps.eps")),
          plot_grid(dendro_plot, clust_heatmap_plot,
          nrow = 2, ncol = 1, rel_heights = c(1,2), align = 'v', axis = 'tb'),
          ncol = 1, device = 'eps',
          base_height = 8, base_aspect_ratio = 0.95)

################################################################################
# save plot objects
save.image(file = file.path(wd, 'output', 'fig2.RData'))

################################################################################

