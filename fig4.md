# Fig.4 Overview of KO reponse

```
export ROOT=/lustre/scratch117/maz/team31/projects/mouse_DMDD
```

## Fig. 4a
KO response plot
plot is produced by fig2.R

Fig. 4b - Heatmap

Get genes that appear in four or more mutants
```
# create a file of which genes are DE in which lines
for mut in $( cut -f1 output/KOs_ordered_by_delay.txt )
do
for comparison in hom_vs_het_wt het_vs_wt hom_vs_het
do
file=data/mutant_response/$mut-deseq2-blacklist-adj-gt-adj-sex-outliers-$comparison-mutant_response.sig.tsv
if [[ -e $file ]]; then
  lines=$(wc -l $file | awk '{print $1}')
  if [[ $lines -eq 1 ]]; then
    echo "$file: NO SIG GENES" 1>&2
  else
    grep ENSMUSG $file | \
    perl -F"\t" -lane 'print join("\t", $F[0], "'$mut'", "1");'
  fi
  break
fi
done
done > output/mutant_response-sig_genes-upset.tmp

# subset to ciliopathy mutants
grep -E 'B9d2|Cc2d2a|Ift140|Kif3b|Kifap3|Nek9|Rpgripl1' output/mutant_response-sig_genes-upset.tmp > output/mutant_response-sig_genes-BCIKKNR-upset.tmp

# reshape from long format to wide
for set in BCIKKNR
do
Rscript reshape-long_to_wide.R output/mutant_response-sig_genes-$set-upset.tmp output/mutant_response-sig_genes-$set-upset.tsv
done

# get genes that are DE in at least 4 lines
Rscript get_intersection_genes.R \
output/mutant_response-sig_genes-$set-upset.tsv \
output/mutant_response-sig_genes-BCIKKNR-intersection-4.txt \
4

# get log2fc data
echo -e "Line\tGene ID\tGene Name\tlog2fc\tClass" > data/fig4b_log2fc.tsv
for mut in B9d2 Cc2d2a Rpgripl1 Kif3b Kifap3 Ift140 Nek9
do
for gene in $( cat output/mutant_response-sig_genes-BCIKKNR-intersection-4.txt )
do
line=$(grep $gene data/mutant_response/$mut*tsv)
num_lines=$(echo $line | grep -c ^ENS)
if [[ $num_lines -eq 0 ]];then
  echo -e "$mut\t$gene\tNA\tNA\tNA"
else
echo $line | perl -lane 'BEGIN{
%class_for = (
Cyp26c1 => "Downstream of Shh signalling",
Dmrt2 => "Shh signalling interactors",
Emx2 => "Downstream of Shh signalling",
Foxa1 => "Downstream of Shh signalling",
Foxa2 => "Downstream of Shh signalling",
Gli1 => "Downstream of Shh signalling",
Gm3764 => "Novel",
Gm38103 => "Novel",
"Mir9-3hg" => "Downstream of Shh signalling",
"Nkx2-1" => "Downstream of Shh signalling",
"Nkx2-2" => "Downstream of Shh signalling",
"Nkx2-9" => "Downstream of Shh signalling",
Pax6 => "Downstream of Shh signalling",
Phox2b => "Downstream of Shh signalling",
Pitx2 => "Downstream of Shh signalling",
Shh => "Downstream of Shh signalling",
Slit1 => "Downstream of Shh signalling",
Sox21 => "Novel",
Stmn2 => "Shh signalling interactors",
Tbx15 => "Downstream of Shh signalling",
Wnt8b => "Shh signalling interactors",
); }
{ print join("\t", "'$mut'", @F[0,9,3], $class_for{$F[9]}, ); }'
fi
done
done >> data/fig4b_log2fc.tsv
```

Produce heatmap
```
Rscript fig4.R \
data/fig4b_log2fc.tsv
```

Fig. 4d - Zkscan17 heatmap

```
# genes to plot are in data/fig4d_data_go.tsv
mut=Zkscan17
countsFile=data/mutant_response/Zkscan17-deseq2-blacklist-adj-gt-adj-sex-outliers-hom_vs_het_wt-mutant_response.sig.tsv
head -n1 $countsFile > output/fig4d.tsv
for geneId in $( cut -f1 data/fig4d_data_go.tsv | grep -v Gene )
do
grep $geneId $countsFile
done >> output/fig4d.tsv

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

Rscript reshape-norm_counts.R output/fig4d.tsv output/fig4d_long_scaled.tsv

# Remove Zkscan17_ from sample names
mv output/fig4d_long_scaled.tsv output/fig4d_long_scaled.tmp
sed -e 's|Zkscan17_||' output/fig4d_long_scaled.tmp > output/fig4d_long_scaled.tsv
rm output/fig4d_long_scaled.tmp

# download heatmap script
curl -LO https://github.com/richysix/bioinf-gen/raw/master/matrix_heatmap_plot.R

# plot heatmap
Rscript matrix_heatmap_plot.R \
--output_file plots/fig4d_heatmap.eps \
--x_column Sample --y_column Gene.ID --y_labels_column Name \
--data_column Normalised.Counts.Scaled --data_axis_label 'Expression Level' \
--colour_palette diverging --colour_legend_position bottom \
--width 7.5 --height 5.5 \
--reverse_y_axis --x_axis_position top --header \
--output_data_file output/fig4d-heatmap.rda \
output/fig4d_long_scaled.tsv

# combine heatmap with plot of categories
Rscript fig4d.R \
data/fig4d_data_go.tsv output/fig4d-heatmap.rda

```

Fig 6c
Count Plots
Make count file for plots
Merge baseline counts with Nadk2 sample counts

```
# Alas2 = ENSMUSG00000025270
# Hbb-y = ENSMUSG00000052187
echo "data/counts/Nadk2-deseq2-blacklist-adj-gt-adj-sex-outliers.tsv
output/Mm_GRCm38_e88_baseline.tsv" > output/Nadk2-counts-files.txt

# run script to make count files
perl merge_deseq_counts.pl \
output/Nadk2-counts-files.txt \
 > output/Nadk2-baseline-counts.tsv
```

Make a sample file for Nadk2 plus baseline samples
```
grep -E 'condition|Nadk2' data/counts/samples-gt-gender-stage-somites.txt | \
perl -F"\t" -lane 'if($. == 1){print join("\t", @F,); }
else{ $category = $F[1] eq "het" ? "het_wt" : $F[1] eq "wt" ? "het_wt" : $F[1];
print join("\t", $F[0], $category, @F[2..4], ); }' > output/Nadk2-samples.txt

grep -E 'condition|Nadk2' data/counts/samples-gt-gender-stage-somites.txt | \
cut -f4 | grep -v stage | sort -u | grep -f - output/samples-Mm_GRCm38_e88_baseline.txt | \
awk -F"\t" 'BEGIN{OFS = "\t"} { print $1, "baseline", $3, $4, $5 }' \
 >> output/Nadk2-samples.txt

# normalise RNA-seq
Rscript normalise_rnaseq.R \
output/Nadk2-baseline-counts.tsv output/Nadk2-samples.txt \
output/Nadk2-baseline-normalised-counts.tsv
```

Get counts for the genes in Fig. 6c
```
for gene in '^Gene' ENSMUSG00000025270 ENSMUSG00000052187
do
grep $gene output/Nadk2-baseline-normalised-counts.tsv
done > output/Nadk2-counts.tsv
```

Run count plot script
```
Rscript graph_rnaseq_counts.R \
output/Nadk2-counts.tsv output/Nadk2-samples.txt \
plots/Nadk2_count_plots.eps default condition stage
```
