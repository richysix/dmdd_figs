## Fig S4a

Count number of repeats in DE in each line

```
# Download data from figshare
mkdir data/repeats/repeatmasker
curl -L --output data/repeats/repeatmasker/deseq2-repeatmasker-all-adj-gt-adj-sex-outliers-repeats-sig.tgz https://ndownloader.figshare.com/files/16370093
cd data/repeats/repeatmasker
tar -xzvf deseq2-repeatmasker-all-adj-gt-adj-sex-outliers-repeats-sig.tgz
cd $ROOT

# count repeats in DE
grep -v Name data/repeats/repeatmasker/*-deseq2-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
sed 's|-deseq2-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv:.*$||g' | \
sed -e 's|data/repeats/repeatmasker/||' | \
sort | uniq -c | awk 'BEGIN{print "KO_line" "\t" "Num_repeats"} {print $2 "\t" $1}' \
 > output/num_repeats_by_line.tsv
```

## Fig S4b

Get number of repeats in DE in each line type and direction
```
echo -e "KO_Line\tType\tOrientation\tCount" > output/repeats_by_line_by_type.tsv
awk -F "\t" '($5 == $12) {print $8}' data/repeats/solely-introns-repeats-all.txt | \
grep -Ff - data/repeats/*-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
sed 's|^data/repeats/||; s|-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv:|\t|g' | \
cut -f 1,10 | awk -F "/" '{print $1}' | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2, $3, "sense", $1}' \
 >> output/repeats_by_line_by_type.tsv

awk -F "\t" '($5 != $12) {print $8}' data/repeats/solely-introns-repeats-all.txt | \
grep -Ff - data/repeats/*-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
sed 's|^data/repeats/||; s|-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv:|\t|g' | \
cut -f 1,10 | awk -F "/" '{print $1}' | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2, $3, "antisense", $1}' \
 >> output/repeats_by_line_by_type.tsv
```

## Fig S4d

Get genes that are enriched for repeats in de in introns and get log2fc and
padj for volcano plot

```
echo -e "Repeat_id\tType\tGene_id\tDhx35\tpval\tpadj\tlog2fc" > output/repeats-for_volcano_plot.tsv
sort -t$'\t' -k1,1 output/fig5c-plot_data.tsv | \
join -t$'\t' - <(sort -t$'\t' -k1,1 data/repeats/solely-introns-repeats-all.txt) | \
cut -f1,4,20 | sort -t$'\t' -k3,3 | join -t$'\t' -1 3 - \
<(sort -t$'\t' -k1,1 data/repeats/Dhx35-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv) | \
perl -F"\t" -lane 'next if $F[4] eq "NA"; if( $F[4] < 0.05 ){ print join("\t", @F[0,10,1..5]); }' | \
sort -uk1,1 >> output/repeats-for_volcano_plot.tsv
```

## Fig S4

Run SupplFig-repeats scripts

```
Rscript SupplFig-repeats.R \
output/num_repeats_by_line.tsv \
output/repeats_by_line_by_type.tsv \
output/repeats-for_volcano_plot.tsv \
data/repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv
data/repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.tsv
output/repeats-enriched_families.tsv
```
