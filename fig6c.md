# Fig 6c
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
