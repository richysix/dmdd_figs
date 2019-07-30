# Fig.2

## Delay PCA

Make a merged count file for all genes

```
# first create the list of all files
find data/counts/ | grep tsv$ | sort > output/all_mutants-counts-files.txt

# download merge script
curl -LO https://github.com/richysix/bioinf-gen/raw/master/merge_deseq_counts.pl

# run script to make count files
perl merge_deseq_counts.pl \
output/all_mutants-counts-files.txt \
 > output/all_samples_merged-counts.tsv

# make samples file containing just somite number for PCA script
cut -f1,4 output/sample_info.txt | sed -e 's|stage|condition|' > output/samples-somites.txt
```

Download PCA script from another repository and run PCA
The script produces a series of plots of successive principal components plotted
against each other. In the PC3 vs PC2 plot the association with somite number can clearly be seen

````
curl -LO https://github.com/iansealy/bio-misc/raw/master/pca_rnaseq.R

gene=50000
transform=vst
Rscript pca_rnaseq.R \
output/all_samples_merged-counts.tsv \
output/samples-somites.txt \
plots/all_mutants-pca.pdf \
output/all_mutants.$transform.$genes $transform $genes 1 0
```

To reproduce the image in Fig. 2a.
Produce PCA plot of PC2 vs PC3 using the [PCA plot](https://richysix.shinyapps.io/pca_plot)
Shiny app with genotype as shape and somite_number as fill colour and save as rda file.
Then edit the plot in an interactive R session.

```
R
```

```R
library('ggplot2')
library('viridis')
load('output/pca-all_genes-top50000-somite_number-pc3_pc2.rda')
# plot
plot_data <- pca_plot$data
# remove samples with NA as somite number
plot_data <- plot_data[ !is.na(plot_data$somite_number), ]

# plot homs separately from hets & wts
# set limits to same for each plot
min_x <- floor(min(plot_data$PC3))
max_x <- ceiling(max(plot_data$PC3))
min_y <- floor(min(plot_data$PC2))
max_y <- ceiling(max(plot_data$PC2))

# make new factor for hom vs het-wt
plot_data$gt_group <- as.character(plot_data$condition)
plot_data$gt_group[ plot_data$gt_group == 'wt' ] <- 'wt-het'
plot_data$gt_group[ plot_data$gt_group == 'het' ] <- 'wt-het'
plot_data$gt_group <- factor(plot_data$gt_group, level = c('hom', 'wt-het'))

postscript(file = file.path('plots', paste0('pca-all_genes-top50000-somite_number-hom_vs_het_wt.eps') ),
            width = 15, height = 7, paper = 'special')
print(ggplot(data = plot_data) +
    geom_point(aes(x = PC3, y = PC2, fill = somite_number, shape = condition), size = 2, stroke = 0.2) +
    scale_fill_viridis(name = 'Somite number', guide = guide_colourbar(order = 1) ) +
    scale_shape_manual(values = c(21,22,23), name = 'Genotype', guide = guide_legend(order = 2)) +
    labs(x = 'PC3 (6.3% variance)', y = 'PC2 (7.7% variance)') +
    facet_wrap( ~ gt_group, nrow = 1) +
    theme_minimal() +
    theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14))
            )
dev.off()

# plot PC values against somite_number
pdf(file = file.path('plots', 'pca-all_genes-top50000-pcs_vs_somite_number.pdf'))
for (pc in paste0('PC', 1:8)) {
    cor_plot <- ggplot(data = plot_data,
                        aes_(x = as.name(pc), y = quote(somite_number))) +
        geom_point() +
        geom_smooth(method = "lm")
    print(cor_plot)
}
dev.off()

# test all components for relationship with somite number
sink('output/pca-all_genes-top50000-pcs_vs_somite_number.txt')
for (pc in paste0('PC', 1:8)) {
    pc_somite_lm <- lm(somite_number ~ get(pc), data = plot_data)
    cat(pc, ' vs somite number\n')
    print(summary(pc_somite_lm))
    cat('\n')
}
sink()

# calculate Pearson cor coefficient for PC3
cor(plot_data$somite_number, plot_data$PC3, method = 'spearman')
[1] -0.8648007

# load PC3 data
pc3_var <- read.delim(file.path('output', 'all_mutants.vst.50000-PC3.tsv'), check.names = FALSE)
pc3_var <- pc3_var[ , c('Gene ID', '% variance explained', 'Name', 'Description',
                        'Chr', 'Start', 'End', 'Strand', 'Biotype')]
pdf(file = file.path('plots', 'all_mutants.vst.50000-PC3.variance-explained.pdf'))
ggplot(data = pc3_var, aes(x = 1:nrow(pc3_var), y = cumsum(`% variance explained`), group = 1) ) +
 geom_line() + labs(x = 'Gene', y = 'Cumulative % Var explained') +
dev.off()
pc_threshold <- 50
cutoff <- which(cumsum(pc3_var[['% variance explained']]) > pc_threshold)[1]
top_genes <- pc3_var[ 1:cutoff, ]
write.table(top_genes, file = file.path('output', 'pc3-top_genes.tsv'),
quote = FALSE, row.names = FALSE, sep = "\t")
```

## Fig 2c

Merge baseline counts with Brd2 sample counts
```
echo "data/counts/Brd2-deseq2-blacklist-adj-gt-adj-sex-outliers.tsv
output/Mm_GRCm38_e88_baseline.tsv" > output/Brd2-counts-files.txt

# run script to make count files
perl merge_deseq_counts.pl \
output/Brd2-counts-files.txt \
 > output/Brd2-baseline-counts.tsv
```

Make a sample file for Brd2 plus stage-matched baseline samples
```
grep -E 'condition|Brd2' data/counts/samples-gt-gender-stage-somites.txt \
 > output/Brd2-samples.txt

grep -E 'condition|Brd2' data/counts/samples-gt-gender-stage-somites.txt | \
cut -f4 | grep -v stage | sort -u | grep -f - output/samples-Mm_GRCm38_e88_baseline.txt | \
awk -F"\t" 'BEGIN{OFS = "\t"} { print $1, "baseline", $3, $4, $5 }' \
 >> output/Brd2-samples.txt
```

Download normalisation script and normalise Brd2 and baseline counts together
```
curl -LO https://github.com/iansealy/bio-misc/raw/master/normalise_rnaseq.R
Rscript normalise_rnaseq.R \
output/Brd2-baseline-counts.tsv output/Brd2-samples.txt \
output/Brd2-baseline-normalised-counts.tsv
```

Get counts for the genes in Fig. 2c
```
head -n1 output/Brd2-baseline-normalised-counts.tsv > output/Brd2-counts.tsv
for gene in '^Gene' ENSMUSG00000032607 ENSMUSG00000022454 ENSMUSG00000035835 ENSMUSG00000021279 
do
grep $gene output/Brd2-baseline-normalised-counts.tsv
done >> output/Brd2-counts.tsv
```

Download count plot script from another repository and run
```
curl -LO https://github.com/iansealy/bio-misc/raw/master/graph_rnaseq_counts.R

Rscript graph_rnaseq_counts.R \
output/Brd2-counts.tsv output/Brd2-samples.txt \
plots/Brd2_count_plots.eps default condition stage
```

## Supplementary Fig. 2

Plot Nell2 (ENSMUSG00000022454) counts with all baseline samples for Supplementary
Figure 2.
Make a sample file for Brd2 plus all baseline samples
```
grep -E 'condition|Brd2' data/counts/samples-gt-gender-stage-somites.txt \
 > output/Brd2-all-samples.txt

grep -v condition output/samples-Mm_GRCm38_e88_baseline.txt | \
awk -F"\t" 'BEGIN{OFS = "\t"} { print $1, "baseline", $3, $4, $5 }' \
 >> output/Brd2-all-samples.txt
```

Normalise counts together and make count plot

```
Rscript normalise_rnaseq.R \
output/Brd2-baseline-counts.tsv output/Brd2-all-samples.txt \
output/Brd2-all-baseline-normalised-counts.tsv

# get Nell2 counts
grep -E '^Gene|ENSMUSG00000022454' output/Brd2-all-baseline-normalised-counts.tsv \
 > output/Brd2-Nell2-counts.tsv

Rscript graph_rnaseq_counts.R \
output/Brd2-Nell2-counts.tsv output/Brd2-all-samples.txt \
plots/Brd2_Nell2_count_plot.eps default condition stage
```

Download the sig gene lists
```
echo "deseq2-blacklist-adj-gt-adj-sex-outliers-mutant_response-sig.tgz https://ndownloader.figshare.com/files/12207071
deseq2-blacklist-adj-gt-adj-sex-outliers-delay-sig.tgz https://ndownloader.figshare.com/files/12207335
deseq2-blacklist-adj-gt-adj-sex-outliers-no_delay-sig.tgz https://ndownloader.figshare.com/files/12207476
deseq2-blacklist-adj-gt-adj-sex-outliers-discard-sig.tgz https://ndownloader.figshare.com/files/12207617" > output/sig-downloads.txt

ROOT=$( pwd )
for set in mutant_response delay no_delay discard
do
mkdir data/$set
file=$( grep "\-$set\-" output/sig-downloads.txt | awk '{print $1}' )
link=$( grep "\-$set\-" output/sig-downloads.txt | awk '{print $2}' )
curl -L --output data/$set/$file $link
cd data/$set
tar -xzvf $file
cd $ROOT
done
```

## EMAPA

Collect together all EMAPA results

```
echo "emapa-mutant_response-results.tgz https://ndownloader.figshare.com/files/12244730
emapa-delay-results.tgz https://ndownloader.figshare.com/files/12244829
emapa-no_delay-results.tgz https://ndownloader.figshare.com/files/12244892
emapa-discard-results.tgz https://ndownloader.figshare.com/files/12244976" > output/emapa-downloads.txt

# make directories
mkdir data/emap/
ROOT=$( pwd )
for set in mutant_response delay no_delay discard
do
mkdir data/emap/$set
file=$( grep "\-$set\-" output/emapa-downloads.txt | awk '{print $1}' )
link=$( grep "\-$set\-" output/emapa-downloads.txt | awk '{print $2}' )
curl -L --output data/emap/$set/$file $link
cd data/emap/$set
tar -xzvf $file
cd $ROOT
done

# collect data for all mutants for each set
# make config file and run with collate_emap_results
echo -e "Mutant\tComparison\tSet\tFile" > output/sets_to_use.tsv
for set in mutant_response delay no_delay discard
do
for file in $( find data/emap/$set/ -type f -name "*tsv" )
do
mut=$( echo $file | sed -e 's|^.*/||; s|\-.*$||' )
comparison=$( echo $file | sed -e 's|^.*/.*outliers\-||; s|\-.*$||' )
echo -e "$mut\t$comparison\t$set\t$file" >> output/sets_to_use.tsv
done
done

# run script to collect together results from appropriate comparisons
perl collate_emap_results.pl --header \
output/sets_to_use.tsv > output/emap_results.all.tsv

```

process EMAPA terms to collapse

```
curl -L --output data/emapa.obo http://purl.obolibrary.org/obo/emapa.obo

Rscript process_emap.R data/emap_results.all.tsv \
data/EMAPA_Nov_17_2017.obo data/root_terms.txt
```
The process_emap script outputs a file (output/duplicated_terms.tsv) of terms
that are children of two different parent terms, so add a column called "to_use"
of 1s and 0s indicating which parent to use.
Save it as output/duplicated_terms-edited.tsv


Get numbers of significant genes
```
data/sig_gene_counts.tsv
```

Calculate the pairwise gene overlap from the delay lists for all delayed mutants
```
mutants=( $(sort output/KOs_delayed.txt) )
for mut in ${mutants[@]}; do
file=data/delay/$mut-deseq2-blacklist-adj-gt-adj-sex-outliers-hom_vs_het_wt-delay.sig.tsv
if [[ -e $file ]]; then
  cut -f1 $file | grep -v '^Gene ID' > data/delay/$mut-delay.sig_genes.txt
fi
done

cat /dev/null > output/delay.sig_genes.err
for mut1 in $( seq 0 $(( ${#mutants[@]} - 1 )) )
do
for mut2 in $( seq 0 $(( ${#mutants[@]} - 1 )) )
do
if [ "$mut1" -lt "$mut2" ]; then
  file1=data/delay/${mutants[$mut1]}-delay.sig_genes.txt
  file2=data/delay/${mutants[$mut2]}-delay.sig_genes.txt
  if [[ -e $file1 && -e $file2 ]]; then
    file1_count=$( cat $file1 | wc -l )
    file2_count=$( cat $file2 | wc -l )
    if [[ $file1_count -gt 1 && $file2_count -gt 1 ]]; then
      perl ~/checkouts/team31/scripts/intersect_lists.pl --counts --jaccard \
        --header --output_no_header --no_venn \
        $file1 $file2
    else
      echo "ONE OF THE FILES HAS NO RESULTS: ${mutants[$mut1]}-hom_vs_het_wt ${mutants[$mut2]}-hom_vs_het_wt"  >> output/delay.sig_genes.err
    fi
  else
    echo "ONE OF THE FILES DOES NOT EXIST: ${mutants[$mut1]}-hom_vs_het_wt ${mutants[$mut2]}-hom_vs_het_wt" >> output/delay.sig_genes.err
  fi
fi
done
done | sed -e 's|data/delay/||g' | sed -e 's|-delay.sig_genes.txt||g' \
 > output/delay-hom_vs_het_wt-sig_genes.out

# cluster all based on gene list overlaps
Rscript cluster_by_overlap.R \
--output_base plots/delay-jaccard-all \
--cluster_methods ward.D2 \
--output_data_file output/delay-jaccard-all.rda \
output/delay-hom_vs_het_wt-sig_genes.out
```

Run fig2.R script

```
Rscript fig2.R \
data/go_results.tsv \
data/dmdd-genes.txt \
output/KOs_ordered_by_delay.txt \
data/sig_gene_counts.tsv \
data/emap_results.all.tsv \
output/duplicated_terms-edited.tsv \
output/delay-jaccard-all.rda
```

## Fig. S3

```
mut=Fcho2
echo "data/counts/$mut-deseq2-blacklist-adj-gt-adj-sex-outliers.tsv
output/Mm_GRCm38_e88_baseline.tsv" > output/$mut-counts-files.txt

# run script to make count files
perl merge_deseq_counts.pl \
output/$mut-counts-files.txt \
 > output/$mut-baseline-counts.tsv
```

Make a sample file for $mut plus stage-matched baseline samples
```
grep -E "condition|$mut" data/counts/samples-gt-gender-stage-somites.txt \
 > output/$mut-samples.txt

grep -E "condition|$mut" data/counts/samples-gt-gender-stage-somites.txt | \
cut -f4 | grep -v stage | sort -u | grep -f - output/samples-Mm_GRCm38_e88_baseline.txt | \
awk -F"\t" 'BEGIN{OFS = "\t"} { print $1, "baseline", $3, $4, $5 }' \
 >> output/$mut-samples.txt
```

Normalise Fcho2 and baseline counts together
```
Rscript normalise_rnaseq.R \
output/$mut-baseline-counts.tsv output/$mut-samples.txt \
output/$mut-baseline-normalised-counts.tsv
```

Get counts for the genes in Fig. S3a-b
```
head -n1 output/$mut-baseline-normalised-counts.tsv > output/$mut-counts.tsv
for gene in '^Gene' ENSMUSG00000034486 ENSMUSG00000020717
do
grep $gene output/$mut-baseline-normalised-counts.tsv
done >> output/$mut-counts.tsv
```

Run count plot script
```
Rscript graph_rnaseq_counts.R \
output/$mut-counts.tsv output/$mut-samples.txt \
plots/${mut}-count_plots.eps default condition stage
```
