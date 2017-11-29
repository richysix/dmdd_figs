plots/Figure2.eps: fig2.R data/Mm_GRCm38_e88_baseline.rda \
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

