Fig. 7b
The data for the plot is in data/repeats-location_vs_genes.tsv

```
Rscript fig7b.R data/repeats-location_vs_genes.tsv
```

Fig. 7c
Download files for repeats analysis
```
mkdir data/repeats/
cd data/repeats/
curl -L --output deseq2-blacklist-adj-gt-adj-sex-outliers-notranscriptome-repeats-sig.tgz https://ndownloader.figshare.com/files/16370780
tar -xzvf deseq2-blacklist-adj-gt-adj-sex-outliers-notranscriptome-repeats-sig.tgz

curl -L --output deseq2-blacklist-adj-gt-adj-sex-outliers-notranscriptome-genes-sig.tgz https://ndownloader.figshare.com/files/11865350
tar -xzvf deseq2-blacklist-adj-gt-adj-sex-outliers-notranscriptome-genes-sig.tgz

cd $ROOT
```

The data on which repeats are in introns is in data/repeats/solely-introns-repeats-all.txt

Get gene ids for genes with repeats in introns that are in DE
```
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $8}' data/repeats/solely-introns-repeats-all.txt | \
grep -Ff - data/repeats/Dhx35-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
cut -f 1 | sort -u | grep -Ff - data/repeats/solely-introns-repeats-all.txt | \
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $1}' | sort -u > output/fig7c-genes-repeats-intron-de_dhx35.txt

# Count total number of repeats in introns
wc -l data/repeats/solely-introns-repeats-all.txt
1639267 data/repeats/solely-introns-repeats-all.txt
```

Count number of repeats in genes
```
grep -Ff output/fig7c-genes-repeats-intron-de_dhx35.txt data/repeats/solely-introns-repeats-all.txt | \
cut -f1 | sort | uniq -c | awk 'BEGIN{OFS = "\t"} {print $2, $1}' > output/fig7c-num_repeats-introns.tsv

# Count number of repeats in introns in de
grep -Ff output/fig7c-genes-repeats-intron-de_dhx35.txt data/repeats/solely-introns-repeats-all.txt | \
cut -f 8 | grep -Ff - data/repeats/Dhx35-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
cut -f 1 | sort -u | grep -Ff - data/repeats/solely-introns-repeats-all.txt | \
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $1}' | sort | \
uniq -c | awk 'BEGIN{OFS = "\t"} {print $2, $1}' > output/fig7c-num_repeats_de-introns.tsv
```

Get adjusted pvalues and log2fc of genes with repeats in introns in de
```
grep -Ff output/fig7c-genes-repeats-intron-de_dhx35.txt \
data/repeats/Dhx35-deseq2-notranscriptome-blacklist-adj-gt-adj-sex-outliers-hom_vs_het_wt.tsv | \
cut -f 1,3,4 | sort -k1,1 > output/fig7c-genes_with_repeats_in_de-pval-log2fc.tsv

# count genes
wc -l output/fig7c-genes_with_repeats_in_de-pval-log2fc.tsv
470 output/fig7c-genes_with_repeats_in_de-pval-log2fc.tsv

# two genes are missing because they are blacklisted
# join to gene list and note as blacklisted
join -v1 output/fig7c-genes-repeats-intron-de_dhx35.txt output/fig7c-genes_with_repeats_in_de-pval-log2fc.tsv | \
awk '{print $1 "\tblacklist\tblacklist"}' >> output/fig7c-genes_with_repeats_in_de-pval-log2fc.tsv

```

Get gene ids for repeats in introns in de in all mutants tested
```
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $8}' data/repeats/solely-introns-repeats-all.txt | \
grep -Fhf - data/repeats/*-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
cut -f1 | sort -u | grep -Ff - data/repeats/solely-introns-repeats-all.txt | \
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $1}' | sort -u > output/fig7c-genes-repeats-intron-de_all.txt 

# get gene ids for repeats in introns in de in everything except Dhx35
for file in $(ls data/repeats/*-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
grep -v Dhx35)
do
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $8}' data/repeats/solely-introns-repeats-all.txt | \
grep -Fhf - $file | cut -f1 | sort -u | \
grep -Ff - data/repeats/solely-introns-repeats-all.txt | \
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $1}'
done | sort -u > output/fig7c-genes-repeats-intron-de_all_but_dhx35.txt

# get genes with repeats in DE in Dhx35 only and ones in Dhx35 and others
comm -1 -3 output/fig7c-genes-repeats-intron-de_all_but_dhx35.txt output/fig7c-genes-repeats-intron-de_dhx35.txt | \
awk 'BEGIN{OFS = "\t"} {print $1, "Dhx35_only" }' > output/fig7c-Line-overlap.tsv
comm -1 -2 output/fig7c-genes-repeats-intron-de_all_but_dhx35.txt output/fig7c-genes-repeats-intron-de_dhx35.txt | \
awk 'BEGIN{OFS = "\t"} {print $1, "Dhx35_plus" }' >> output/fig7c-Line-overlap.tsv

# join files together
perl -le 'print join("\t", qw{gene_id repeats_intron repeats_intron_de source gene_padj gene_log2fc});' > output/fig7c-for-enrichment-test.tsv
join -t$'\t' output/fig7c-num_repeats-introns.tsv output/fig7c-num_repeats_de-introns.tsv | \
join -t$'\t' - <(sort -k1,1 output/fig7c-Line-overlap.tsv) | \
join -t$'\t' - <(sort -k1,1 output/fig7c-genes_with_repeats_in_de-pval-log2fc.tsv) >> output/fig7c-for-enrichment-test.tsv
```

Get repeat log2fc and gene_id info together
```
# get number of columns in sig file
file1_cols=$( head -n1 data/repeats/Dhx35-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
tr '\t' '\n' | wc -l )
# calculate col numbers for gene and repeat strand
gene_strand=$(( $file1_cols + 5 - 1 ))
repeat_strand=$(( $file1_cols + 11 - 1 ))

# join repeat and gene info together
perl -le 'print join("\t", qw{repeat_id pval padj log2fc chr start end strand
name gene_id gene_chr gene_start gene_end gene_strand gene_biotype gene_name});' > output/fig7c-dhx35-repeats-introns-genes-sig.tsv
sort -t$'\t' -k1,1 data/repeats/Dhx35-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
join -t$'\t' -2 8 - <(sort -t$'\t' -k8,8 data/repeats/solely-introns-repeats-all.txt) | \
perl -F"\t" -lane 'if($F['$gene_strand'] eq $F['$repeat_strand']) {
print join("\t", @F[0..8,'$file1_cols'..'$(( $file1_cols + 6 ))'] ); }' \
 >> output/fig7c-dhx35-repeats-introns-genes-sig.tsv
```

Run Fig. 7c script
```
Rscript fig7c.R output/fig7c-for-enrichment-test.tsv \
$( wc -l data/repeats/solely-introns-repeats-all.txt | awk '{print $1}' ) \
output/fig7c-dhx35-repeats-introns-genes-sig.tsv
```

Fig. 7d

All repeats counted by family
```
awk -F ":" '$1 !~/adjp/ { group = $2; sub(/\/.*$/, "", group); print $1 "\t" group "\tall"}' \
data/repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.tsv | \
sort | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' > output/all-repeats-count.tsv
```

Repeats in DE counted by family
```
awk -F ":" '$1 !~/adjp/ { group = $2; sub(/\/.*$/, "", group); print $1 "\t" group "\tde"}' \
data/repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
sort | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' > output/all-repeats-in-DE-count.tsv
```

All repeats longer than 4 kbp counted by family
```
awk '$3 !~/adjp/ {print $7-$6"\t"$0}' \
data/repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.tsv | \
awk -F "\t" '$1>=4000 {print $2}' | awk -F ":" '{ group = $2; sub(/\/.*$/, "", group); print $1 "\t" group "\tall"}' | \
sort | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' > output/all-repeats-count-gt=4000.tsv
```

Repeats in DE longer than 4 kbp counted by family

```
awk '$3 !~/adjp/ {print $7-$6"\t"$0}' \
data/repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
awk -F "\t" '$1>=4000 {print $2}' | awk -F ":" '{ group = $2; sub(/\/.*$/, "", group); print $1 "\t" group "\tde"}' | \
sort | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' > output/all-repeats-count-gt=4000-in-DE.tsv
```

Combine data
```
join -t$'\t' <( sort -t$'\t' -k1,1 output/all-repeats-count.tsv ) \
<( sort -t$'\t' -k1,1 output/all-repeats-in-DE-count.tsv ) | \
awk 'BEGIN{ OFS = "\t"; print "Family", "Group", "full_length", "repeats", "de", "not_de" }
{ print $1, $2, "all", $4, $7, $4 - $7 }' > output/fig7d_repeats_de.tsv

join -t$'\t' <( sort -t$'\t' -k1,1 output/all-repeats-count-gt=4000.tsv ) \
<( sort -t$'\t' -k1,1 output/all-repeats-count-gt=4000-in-DE.tsv ) | \
awk 'BEGIN{ OFS = "\t" }
{ print $1, $2, "full_length", $4, $7, $4 - $7 }' >> output/fig7d_repeats_de.tsv
```

Make sig files for heatmaps

```
# L1MdGf_I
grep -E 'L1MdGf_I:|adjp' data/notranscriptome-repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
awk '{print $1 "\t" $0}' | sed -e 's|^Name|Gene ID|' > output/L1MdGf_I_hom_vs_het_wt.sig.tsv

# MMERGLN-int
grep -E 'MMERGLN-int:|adjp' data/notranscriptome-repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
awk '{print $1 "\t" $0}' | sed -e 's|^Name|Gene ID|' > output/MMERGLN-int_hom_vs_het_wt.sig.tsv

# MMETn-int
grep -E 'MMETn-int:|adjp' data/notranscriptome-repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
awk '{print $1 "\t" $0}' | sed -e 's|^Name|Gene ID|' > output/MMETn-int_hom_vs_het_wt.sig.tsv

# samples file
grep -E 'condition|Morc2a' data/counts/samples-gt-gender-stage-somites.txt > output/Morc2a-samples.txt
```

The heatmaps were produced with the [geneExpr](https://richysix.shinyapps.io/geneexpr/) Shiny App
and the above count and sample files.

Locations of repeats
```
# total numbers of repeats are in data/repeat-location-all.tsv

# Exons only (one repeat can cover multiple exons, so uniq repeats)
cut -f1 data/notranscriptome-repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
grep -Ff - data/overlapping-exons-repeats-all.txt | \
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $8}' | \
sort -u | awk -F ":" '{print $1}' | sort | uniq -c | \
awk -F " " '{print "exon\t" $1 "\t" $2}' > output/repeat-location-exon.tsv

# Introns only
cut -f1 data/notranscriptome-repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.sig.tsv | \
grep -Ff - data/repeats/solely-introns-repeats-all.txt | \
awk -F "\t" '($5 ~/+/ && $12 ~/+/) || ($5 ~/-/ && $12 ~/-/) {print $8}' | \
sort -u | awk -F ":" '{print $1}' | sort | uniq -c | \
awk -F " " '{print "intron\t" $1 "\t" $2}' > output/repeat-location-intron.tsv

# Join files together and fill in missing values
join -t$'\t' -a1 -1 3 -2 3 <( sort -t$'\t' -k3,3 data/repeat-location-all.tsv ) \
<( sort -t$'\t' -k3,3 output/repeat-location-exon.tsv ) | \
awk 'BEGIN{OFS = "\t"} {if($4 == ""){ $4 = "NA"} if($5 == ""){ $5 = "NA" } print $0 }' | \
join -t$'\t' -a1 -2 3 - <( sort -t$'\t' -k3,3 output/repeat-location-intron.tsv ) | \
awk 'BEGIN{OFS = "\t"; print "Family", "exon", "intron", "intergenic" }
{if($6 == ""){ $6 = "NA"} if($7 == ""){ $7 = "NA" }
intergenic = $3; if( $5 != "NA" ){ intergenic -= $5 }
if( $7 != "NA" ){ intergenic -= $7 } print $1, $5, $7, intergenic }' \
 > output/fig7d_repeats_location.tsv
```

Run fig7d script

```
Rscript fig7d.R output/fig7d_repeats_de.tsv \
$( grep -v Name data/repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.tsv | wc -l ) \
$( awk 'BEGIN{ sum = 0; } { if( $3 !~/adjp/ && $7-$6 >= 4000 ){sum += 1} }
END{ print sum }' data/repeats/Morc2a-deseq2-notranscriptome-repeatmasker-all-adj-gt-adj-sex-outliers-hom_vs_het_wt.tsv ) \
output/fig7d_repeats_location.tsv
```
