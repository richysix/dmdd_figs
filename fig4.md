```
export ROOT=/lustre/scratch117/maz/team31/projects/mouse_DMDD
```

Fig. 5a
KO response plot
plot is produced by fig4.R

Fig. 5b - Heatmap
Table (fig5_table.txt) is copied from an Excel file
Convert table to long format for heatmap script

```
perl -F"\t" -lane 'if($. == 1){
print join("\t", qw{Gene KO log2fc});
for($i = 1; $i < scalar @F; $i++){
  $col_header[$i] = $F[$i];
} }
else{
  for($i = 1; $i < scalar @col_header; $i++){
    $value = $F[$i] eq "" ? "NA" : $F[$i];
	print join("\t", $F[0], $col_header[$i], $value,);
  }
}' data/fig5b_table.txt > output/fig5b_table.tsv
```

output as pdf, svg and eps
```
for suffix in pdf eps svg
do
export R_LIBS_USER='.R/lib'
/software/R-3.3.0/bin/Rscript \
fig5.R --output_file plots/fig5b-heatmap-truncated.$suffix \
output/fig5b_table.tsv
done

```

Fig. 5d - Zkscan17 heatmap

```
# genes to plot are in data/fig5d_data_go.tsv
mut=Zkscan17
countsFile=$ROOT/lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv
head -n1 $countsFile > output/fig5d.tsv
for geneId in $( cut -f1 data/fig5d_data_go.tsv | grep -v Gene )
do
grep $geneId $countsFile
done >> output/fig5d.tsv

# melt data to long format
# and centre and scale counts
echo 'library("reshape2")
Args            <- commandArgs()
data <- read.delim(Args[6])
count_col_names <- colnames(data)[ grepl("normalised\\.count$", colnames(data)) ]
# subset data to counts
count_data <- as.matrix(data[, count_col_names])
# centre and scale by row
count_data_scaled <- t(scale(t(count_data)))
# add back gene name and gene id columns
scaled_data <- cbind(data[, c("Gene.ID", "Name")], count_data_scaled)
colnames(scaled_data) <- gsub("\\.normalised\\.count$", "", colnames(scaled_data))
count_col_names <- gsub("\\.normalised\\.count$", "", count_col_names)
# melt
data_m <- melt(scaled_data, id.vars = c("Gene.ID", "Name"),
                measure.vars = count_col_names,
                variable.name = "Sample",
                value.name = "Normalised.Counts.Scaled")
write.table(data_m, file = Args[7], quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
' > reshape-norm_counts.R

/software/R-3.3.0/bin/Rscript reshape-norm_counts.R output/fig5d.tsv output/fig5d_long_scaled.tsv

# Remove Zkscan17_ from sample names
mv output/fig5d_long_scaled.tsv output/fig5d_long_scaled.tmp
sed -e 's|Zkscan17_||' output/fig5d_long_scaled.tmp > output/fig5d_long_scaled.tsv
rm output/fig5d_long_scaled.tmp

# plot heatmap
export R_LIBS_USER=.R/lib
/software/R-3.3.0/bin/Rscript \
~rw4/checkouts/team31/scripts/matrix_heatmap_plot.R \
--output_file plots/fig5d_heatmap.eps \
--x_column Sample --y_column Gene.ID --y_labels_column Name \
--data_column Normalised.Counts.Scaled --data_axis_label 'Expression Level' \
--colour_palette diverging --colour_legend_position bottom \
--width 7.5 --height 5.5 \
--reverse_y_axis --x_axis_position top --header \
--output_data_file output/fig5d-heatmap.rda \
output/fig5d_long_scaled.tsv

# combine heatmap with plot of categories
export R_LIBS_USER=.R/lib
/software/R-3.3.0/bin/Rscript fig5d.R \
data/fig5d_data_go.tsv output/fig5d-heatmap.rda

```

Fig 5f
Count Plots

```
# make count file for plots
# Alas2 = ENSMUSG00000025270
# Hbb-y = ENSMUSG00000052187
mut=Nadk2
countsFile=$ROOT/lane-process/$mut/deseq2-baseline-grandhet-blacklist-adj-gt-adj-sex-stage-nicole-definite-maybe-outliers/hom_vs_het_wt.tsv

head -n1 $countsFile > output/fig5f-counts.tsv
grep -E 'ENSMUSG000000(25270|52187)' $countsFile >> output/fig5f-counts.tsv

# samples file
# make new samples file with categories baseline, het_wt and hom to match plots in fig. 2
# one of the "wts" is actually a het, but this doesn't matter as we are combining hets and wts into het_wt.
perl -F"\t" -lane 'if($. == 1){print join("\t", @F, "category"); }
else{ $category = $F[1] eq "het" ? "het_wt" : $F[1] eq "wt" ? "het_wt" : $F[1];
print join("\t", @F, $category, ); }' \
$ROOT/lane-process/$mut/deseq2-baseline-grandhet-blacklist-adj-gt-adj-sex-stage-nicole-definite-maybe-outliers/samples.txt \
 > output/Nadk2-samples.txt

# rerun counts script
Rscript ~rw4/checkouts/bio-misc/graph_rnaseq_counts.R output/fig5f-counts.tsv \
output/Nadk2-samples.txt \
plots/fig5f.eps default category stage

```
