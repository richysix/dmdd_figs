outputs: data/Mm_GRCm38_e88_baseline.rda \
plots/tissue_vs_WE_heatmap.eps output/highly_expressed-tissues.tsv \
output/human-mim.tsv \
plots/Figure2.eps

# Fig2
plots/Figure2.eps: fig2.R data/Mm_GRCm38_e88_baseline.rda \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/samples-minus-outliers.txt \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/ko_expr/ko_expr.tsv 
	/software/R-3.3.0/bin/Rscript fig2.R \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/samples-minus-outliers.txt \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/ko_expr/ko_expr.tsv \
	data/Mm_GRCm38_e88_baseline.rda \
	data/sig_gene_counts.tsv \
	output/human-mim-edited.tsv \
	data/Dr_GRCz10_e90_baseline.rda

# get data from OMIM
output/human-mim.tsv: get_mim_data.R \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/dmdd-genes.txt
	/software/R-3.3.0/bin/Rscript get_mim_data.R \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/dmdd-genes.txt

# suppl fig1 - tissues vs whole embryos
plots/tissue_vs_WE_heatmap.eps output/highly_expressed-tissues.tsv: \
fig1.R data/Mm_GRCm38_e88_baseline.rda /lustre/scratch117/maz/team31/projects/mouse_DMDD/PRJEB4513-E8.25/all.tsv
	export R_LIBS_USER=.R/lib:/software/team31/R-3.3.0/:/software/team31/R; \
	/software/R-3.3.0/bin/Rscript fig1.R data/Mm_GRCm38_e88_baseline.rda \
	/lustre/scratch117/maz/team31/projects/mouse_DMDD/PRJEB4513-E8.25/all.tsv

# baseline data
data/Mm_GRCm38_e88_baseline.rda: baseline.R /lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/baseline-grandhet/all.tsv \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/baseline-grandhet/deseq2/samples.txt
	/software/R-3.3.0/bin/Rscript baseline.R

