outputs: data/Mm_GRCm38_e88_baseline.rda \
plots/tissue_vs_WE_heatmap.eps output/highly_expressed-tissues.tsv \
output/human-mim.tsv \
plots/Figure2.eps \
output/duplicated_terms.tsv \
output/mrna_abnormal-jaccard-all.rda \
output/mean_by_mut_gt.counts.tsv output/samples_by_mut_by_gt.txt \
plots/mrna_abnormal.emap_aggregated_summary.eps \
plots/mrna_abnormal-gene_list-overlaps.eps \
plots/morc2a-repeats.eps

# fig5c
plots/morc2a-repeats.eps: \
data/fig5c_repeats_de_all.tsv \
data/fig5c_repeats_de.tsv \
data/fig5c_repeats_location.tsv
	/software/R-3.3.0/bin/Rscript fig5c.R \
	data/fig5c_repeats_de_all.tsv \
	data/fig5c_repeats_de.tsv \
	data/fig5c_repeats_location.tsv

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
plots/mrna_abnormal.emap_aggregated_summary.eps \
plots/mrna_abnormal-gene_list-overlaps.eps: fig2.R \
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
plots/tissue_vs_WE_heatmap.eps output/highly_expressed-tissues.tsv: \
fig1.R data/Mm_GRCm38_e88_baseline.rda \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/PRJEB4513-E8.25/downsample/deseq2/counts.txt \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/PRJEB4513-E8.25/downsample/deseq2/samples.txt
	export R_LIBS_USER=.R/lib:/software/team31/R-3.3.0/:/software/team31/R; \
	/software/R-3.3.0/bin/Rscript fig1.R data/Mm_GRCm38_e88_baseline.rda \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/PRJEB4513-E8.25/downsample/deseq2/counts.txt \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/PRJEB4513-E8.25/downsample/deseq2/samples.txt

# baseline data
data/Mm_GRCm38_e88_baseline.rda: baseline.R /lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/baseline-grandhet/all.tsv \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/baseline-grandhet/deseq2/samples.txt
	/software/R-3.3.0/bin/Rscript baseline.R

