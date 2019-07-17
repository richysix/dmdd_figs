# Fig.5 Overview of KO reponse

```
export ROOT=/lustre/scratch117/maz/team31/projects/mouse_DMDD
```

## Fig. 5a
KO response plot
plot is produced by fig2.R

Fig. 5b - Heatmap

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
echo -e "Line\tGene ID\tGene Name\tlog2fc\tClass" > data/fig5b_log2fc.tsv
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
done >> data/fig5b_log2fc.tsv
```

Produce heatmap
```
Rscript fig5b.R \
data/fig5b_log2fc.tsv
```

Fig. 5c - Zkscan17 heatmap

```
# genes to plot are in data/fig5c_data_go.tsv
mut=Zkscan17
countsFile=data/mutant_response/Zkscan17-deseq2-blacklist-adj-gt-adj-sex-outliers-hom_vs_het_wt-mutant_response.sig.tsv
head -n1 $countsFile > output/fig5c.tsv
for geneId in $( cut -f1 data/fig5c_data_go.tsv | grep -v Gene )
do
grep $geneId $countsFile
done >> output/fig5c.tsv

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

Rscript reshape-norm_counts.R output/fig5c.tsv output/fig5c_long_scaled.tsv

# Remove Zkscan17_ from sample names
mv output/fig5c_long_scaled.tsv output/fig5c_long_scaled.tmp
sed -e 's|Zkscan17_||' output/fig5c_long_scaled.tmp > output/fig5c_long_scaled.tsv
rm output/fig5c_long_scaled.tmp

# download heatmap script
curl -LO https://github.com/richysix/bioinf-gen/raw/master/matrix_heatmap_plot.R

# plot heatmap
Rscript matrix_heatmap_plot.R \
--output_file plots/fig5c_heatmap.eps \
--x_column Sample --y_column Gene.ID --y_labels_column Name \
--data_column Normalised.Counts.Scaled --data_axis_label 'Expression Level' \
--colour_palette diverging --colour_legend_position bottom \
--width 7.5 --height 5.5 \
--reverse_y_axis --x_axis_position top --header \
--output_data_file output/fig5c-heatmap.rda \
output/fig5c_long_scaled.tsv

# combine heatmap with plot of categories
Rscript fig5c.R \
data/fig5c_data_go.tsv output/fig5c-heatmap.rda

```
