Get data for log2fc/padj for both repeats in DE in introns and the genes they are in
```
# Get gene ids for genes with repeats in introns that are in DE
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $8}' $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
grep -Ff - $ROOT/lane-process/Dhx35/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv | \
cut -f 1 | sort -u | grep -Ff - $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $1}' | sort -u > output/fig5c-genes-repeats-intron-de_dhx35.txt

# total number of repeats in introns
wc -l $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt
1639267 /lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/repenrich2/solely-introns-repeats-all.txt

# Number of repeats in gene
grep -Ff output/fig5c-genes-repeats-intron-de_dhx35.txt $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
cut -f1 | sort | uniq -c | awk 'BEGIN{OFS = "\t"} {print $2, $1}' > output/fig5c-num_repeats-introns.tsv

# number of repeats in introns in de
grep -Ff output/fig5c-genes-repeats-intron-de_dhx35.txt $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
cut -f 8 | grep -Ff - $ROOT/lane-process/Dhx35/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv | \
cut -f 1 | sort -u | grep -Ff - $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $1}' | sort | \
uniq -c | awk 'BEGIN{OFS = "\t"} {print $2, $1}' > output/fig5c-num_repeats_de-introns.tsv

# pvalue and log2fc of genes with repeats in introns in de
grep -Ff output/fig5c-genes-repeats-intron-de_dhx35.txt \
$ROOT/lane-process/Dhx35/deseq2-notranscriptome-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.tsv | \
cut -f 1,3,4 | sort -k1,1 > output/fig5c-genes_with_repeats_in_de-pval-log2fc.tsv

# two genes are missing because they are blacklisted
join -v1 output/fig5c-genes-repeats-intron-de_dhx35.txt output/fig5c-genes_with_repeats_in_de-pval-log2fc.tsv | \
awk '{print $1 "\tblacklist\tblacklist"}' >> output/fig5c-genes_with_repeats_in_de-pval-log2fc.tsv

# get gene ids for repeats in introns in de in all mutants tested
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $8}' $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
grep -Fhf - $ROOT/lane-process/*/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv | \
cut -f1 | sort -u | grep -Ff - $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $1}' | sort -u > output/fig5c-genes-repeats-intron-de_all.txt 

# get gene ids for repeats in introns in de in everything except Dhx35
for file in $(ls $ROOT/lane-process/*/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv | \
grep -v Dhx35)
do
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $8}' $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
grep -Fhf - $file | cut -f1 | sort -u | \
grep -Ff - $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt | \
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $1}'
done | sort -u > output/fig5c-genes-repeats-intron-de_all_but_dhx35.txt

# get genes with repeats in DE in Dhx35 only and ones in Dhx35 and others
comm -1 -3 output/fig5c-genes-repeats-intron-de_all_but_dhx35.txt output/fig5c-genes-repeats-intron-de_dhx35.txt | \
awk 'BEGIN{OFS = "\t"} {print $1, "Dhx35_only" }' > output/fig5c-Line-overlap.tsv
comm -1 -2 output/fig5c-genes-repeats-intron-de_all_but_dhx35.txt output/fig5c-genes-repeats-intron-de_dhx35.txt | \
awk 'BEGIN{OFS = "\t"} {print $1, "Dhx35_plus" }' >> output/fig5c-Line-overlap.tsv

# join files together
perl -le 'print join("\t", qw{gene_id repeats_intron repeats_intron_de source gene_padj gene_log2fc});' > output/fig5c-for-enrichment-test.tsv
join -t$'\t' output/fig5c-num_repeats-introns.tsv output/fig5c-num_repeats_de-introns.tsv | \
join -t$'\t' - <(sort -k1,1 output/fig5c-Line-overlap.tsv) | \
join -t$'\t' - <(sort -k1,1 output/fig5c-genes_with_repeats_in_de-pval-log2fc.tsv) >> output/fig5c-for-enrichment-test.tsv
```

Get repeat log2fc and gene_id info together
```
perl -le 'print join("\t", qw{repeat_id pval padj log2fc chr start end strand
name gene_id gene_chr gene_start gene_end gene_strand gene_biotype gene_name});' > output/fig5c-dhx35-repeats-introns-genes-sig.tsv
sort -t$'\t' -k1,1 $ROOT/lane-process/Dhx35/deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv | \
join -t$'\t' -2 8 - <(sort -t$'\t' -k8,8 $ROOT/lane-process/repenrich2/solely-introns-repeats-all.txt) | \
perl -F"\t" -lane 'if($F[55] eq $F[61]) {print join("\t", @F[0..8,51..57]);}' >> output/fig5c-dhx35-repeats-introns-genes-sig.tsv
```
