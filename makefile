outputs: setup.done \
plots/tissue_vs_WE_heatmap.eps output/highly_expressed-tissues.tsv \
output/tissue_vs_WE_Venn_numbers.tsv output/tissue_only-counts.tsv \
output/human-mim.tsv \
plots/emap-bubble_plots-by_results_set.eps \
output/duplicated_terms.tsv \
output/mrna_abnormal-jaccard-all.rda \
output/mean_by_mut_gt.counts.tsv output/samples_by_mut_by_gt.txt \
plots/mrna_abnormal.emap_aggregated_summary.eps \
plots/mrna_abnormal-gene_list-overlaps.eps \
plots/morc2a-repeats.eps \
plots/fig5c.eps \
plots/fig5b.eps \
plots/fig4-heatmap.eps

# Suppl Fig 3
plots/num_repeats_by_line.eps plots/num_repeats_by_type_combined.eps \
plots/repeats-volcano.eps plots/repeats-length.eps \
plots/repeats-counts_vs_fc.eps: SupplFig-repeats.R output/num_repeats_by_line.tsv \
output/repeats_by_line_by_type.tsv output/repeats-for_volcano_plot.tsv \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Morc2a/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Morc2a/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.tsv \
output/repeats-enriched_families.tsv
	/software/R-3.3.0/bin/Rscript SupplFig-repeats.R \
	output/num_repeats_by_line.tsv \
	output/repeats_by_line_by_type.tsv \
	output/repeats-for_volcano_plot.tsv \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Morc2a/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Morc2a/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.tsv \
	output/repeats-enriched_families.tsv
		   
# fig5d
plots/morc2a-repeats.eps: \
fig5d.R data/fig5d_repeats_de_all.tsv \
data/fig5d_repeats_de.tsv data/fig5d_repeats_location.tsv
	/software/R-3.3.0/bin/Rscript fig5d.R \
	data/fig5d_repeats_de_all.tsv \
	data/fig5d_repeats_de.tsv \
	data/fig5d_repeats_location.tsv

# fig5c
plots/fig5c.eps: \
fig5c.R output/fig5c-for-enrichment-test.tsv \
output/fig5c-dhx35-repeats-introns-genes-sig.tsv
	/software/R-3.3.0/bin/Rscript fig5c.R \
	output/fig5c-for-enrichment-test.tsv \
	output/fig5c-dhx35-repeats-introns-genes-sig.tsv

# fig5b
plots/fig5b.eps: \
fig5b.R data/repeats-location_vs_genes.tsv
	/software/R-3.3.0/bin/Rscript fig5b.R \
	data/repeats-location_vs_genes.tsv

# fig 4b
plots/fig4-heatmap.eps: fig4.R data/fig4b_log2fc.tsv
	/software/R-3.3.0/bin/Rscript fig4.R \
	data/fig4b_log2fc.tsv

# Mean counts by KO gene and genotype for PCA
output/mean_by_mut_gt.counts.tsv output/samples_by_mut_by_gt.txt: \
mean_expression_by_mut_gt.R output/all_samples_merged-counts.tsv \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/samples-minus-outliers.txt \
output/KOs_ordered_by_delay.txt
	/software/R-3.3.0/bin/Rscript mean_expression_by_mut_gt.R \
	output/all_samples_merged-counts.tsv \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/samples-minus-outliers.txt \
	output/KOs_ordered_by_delay.txt

# Overlaps by genes of mrna_abnormal lists
output/mrna_abnormal-jaccard-all.rda: cluster_by_overlap.R \
output/mrna_abnormal-hom_vs_het_wt-sig_genes.out
	/software/R-3.3.0/bin/Rscript cluster_by_overlap.R \
	--output_base plots/mrna_abnormal-PC3-jaccard-all \
	--cluster_methods ward.D2 \
	--output_data_file output/mrna_abnormal-jaccard-all.rda \
	output/mrna_abnormal-hom_vs_het_wt-sig_genes.out

# aggregate EMAPA terms to umbrella terms
output/duplicated_terms.tsv: data/emap_results.all.tsv \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/emap/EMAPA_Nov_17_2017.obo \
data/root_terms.txt
	/software/R-3.3.0/bin/Rscript process_emap.R data/emap_results.all.tsv \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/emap/EMAPA_Nov_17_2017.obo \
	data/root_terms.txt

# Fig3 - sample overview
plots/Sample_overview.eps: fig3.R \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/dmdd-genes.txt \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/samples-minus-outliers.txt \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/ko_expr/ko_expr.tsv \
data/Mm_GRCm38_e88_baseline.rda data/sig_gene_counts.tsv \
output/human-mim-edited.tsv
	export R_LIBS_USER=.R/lib:/software/team31/R-3.3.0/:/software/team31/R; \
	/software/R-3.3.0/bin/Rscript fig3.R \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/dmdd-genes.txt \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/samples-minus-outliers.txt \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/ko_expr/ko_expr.tsv \
	data/Mm_GRCm38_e88_baseline.rda \
	data/sig_gene_counts.tsv \
	output/human-mim-edited.tsv

# get data from OMIM
output/human-mim.tsv: get_mim_data.R \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/dmdd-genes.txt
	/software/R-3.3.0/bin/Rscript get_mim_data.R \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/dmdd-genes.txt

# Fig2
plots/emap-bubble_plots-by_results_set.eps: fig2.R \
data/go_results.tsv /lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/dmdd-genes.txt \
output/KOs_ordered_by_delay.txt data/sig_gene_counts.tsv \
data/emap_results.all.tsv output/duplicated_terms-edited.tsv \
output/mrna_abnormal-jaccard-all.rda
	export R_LIBS_USER=.R/lib:/software/team31/R:/software/team31/R-3.3.0; \
	/software/R-3.3.0/bin/Rscript fig2.R \
	data/go_results.tsv \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/dmdd-genes.txt \
	output/KOs_ordered_by_delay.txt \
	data/sig_gene_counts.tsv \
	data/emap_results.all.tsv \
	output/duplicated_terms-edited.tsv \
	output/mrna_abnormal-jaccard-all.rda

# suppl fig1 - tissues vs whole embryos
plots/tissue_vs_WE_heatmap.eps output/highly_expressed-tissues.tsv \
output/tissue_vs_WE_Venn_numbers.tsv output/tissue_only-counts.tsv: \
fig1.R data/Mm_GRCm38_e88_baseline.rda \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/PRJEB4513-E8.25/downsample/deseq2/counts.txt \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/PRJEB4513-E8.25/downsample/deseq2/samples.txt
	export R_LIBS_USER=.R/lib:/software/team31/R-3.3.0/:/software/team31/R; \
	/software/R-3.3.0/bin/Rscript fig1.R \
	data/PRJEB4513-E8.25/4567_somites-counts.tsv \
	data/PRJEB4513-E8.25/4567_somites-samples.tsv \
	data/PRJEB4513-E8.25/tissues-counts.tsv \
	data/PRJEB4513-E8.25/tissues-samples.tsv

setup.done: setup.R
	/software/R-3.3.0/bin/Rscript setup.R
