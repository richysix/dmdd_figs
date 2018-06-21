Count number of repeats in DE in each line (Fig S4a)

```
grep : $ROOT/lane-process/*/deseq2-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv | \
cut -f 1,3,4,5,8,9 | sed 's/\/deseq2-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers\/hom_vs_het_wt.sig.tsv:/\t/g' | \
sed -e 's|/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/||' | \
cut -f 1 | sort | uniq -c | awk 'BEGIN{print "KO_line" "\t" "Num_repeats"} {print $2 "\t" $1}' \
 > output/num_repeats_by_line.tsv
```

get number of repeats in DE in each line type and direction (Fig S4b)
```
echo -e "KO_Line\tType\tOrientation\tCount" > output/repeats_by_line_by_type.tsv
awk -F "\t" '(($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/)) {print $8}' $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
grep -Ff - $ROOT/lane-process/*/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv | \
sed 's|^/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/||; s|/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers\/hom_vs_het_wt.sig.tsv:|\t|g' | \
cut -f 1,10 | awk -F "/" '{print $1}' | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2, $3, "sense", $1}' \
 >> output/repeats_by_line_by_type.tsv

awk -F "\t" '(($5 ~/+/ && $12 ~/-/) || ($5 ~/-/ && $12 ~/+/)) {print $8}' $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
grep -Ff - $ROOT/lane-process/*/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv | \
sed 's|^/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/||; s|/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers\/hom_vs_het_wt.sig.tsv:|\t|g' | \
cut -f 1,10 | awk -F "/" '{print $1}' | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2, $3, "antisense", $1}' \
 >> output/repeats_by_line_by_type.tsv
```


Fig S4d
```
# get genes that are enriched for repeats in de in introns and grep repeats
# grep repeats from sig file to get log2fc and padj
# plot volcano plot
echo -e "Repeat_id\tType\tGene_id\tDhx35\tpval\tpadj\tlog2fc" > output/repeats-for_volcano_plot.tsv
sort -t$'\t' -k1,1 output/fig5c-plot_data.tsv | \
join -t$'\t' - <(sort -t$'\t' -k1,1 $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt) | \
cut -f1,4,15 | sort -t$'\t' -k3,3 | join -t$'\t' -1 3 - \
<(sort -t$'\t' -k1,1 $ROOT/lane-process/Dhx35/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv) | \
perl -F"\t" -lane 'next if $F[4] eq "NA"; if( $F[4] < 0.05 ){ print join("\t", @F[0,10,1..5]); }' | \
sort -uk1,1 >> output/repeats-for_volcano_plot.tsv
```

