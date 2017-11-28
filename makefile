output/fig2.RData plots/embryo_stage_colour.pdf \
plots/embryo_stage_colour_zero_white.pdf \
plots/embryo_stage_colour_zero_white.eps \
plots/embryo_stage_size_colour.pdf \
plots/embryo_stage_by_gene_by_gt.pdf \
plots/embryo_ko_expr_plot.pdf \
plots/mouse_baseline_heatmap.pdf \
plots/mouse_baseline_theiler_stage_heatmap.pdf \
plots/mouse_baseline_theiler_stage_log10_heatmap.pdf \
plots/sig_genes_heatmap.pdf: fig2.R data/Mm_GRCm38_e88_baseline.rda \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/samples-minus-outliers.txt \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/KO_expr.tsv
	/software/R-3.3.0/bin/Rscript fig2.R \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/samples-minus-outliers.txt \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/KO_expr.tsv \
	data/Mm_GRCm38_e88_baseline.rda \
	data/Dr_GRCz10_e90_baseline.rda \
	data/sig_gene_counts.tsv

data/Mm_GRCm38_e88_baseline.rda: /lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/baseline-grandhet/all.tsv \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/baseline-grandhet/deseq2/samples.txt
	/software/R-3.3.0/bin/Rscript baseline.R

