# Fig.4

```bash
# set up working directory
# change this if you are trying recreate the analysis
# everything else should then be relative to this directory
export ROOT=/lustre/scratch117/maz/team31/projects/mouse_DMDD
```

Run comparison of baseline

```
for mut in $( < ~/checkouts/mouse_dmdd/output/KOs_delayed.txt )
do
echo $mut
ROOT=/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/$mut
for comparison in $( find $ROOT/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers | \
grep done$ | sed -e 's|^.*/||; s|\.done$||' )
do
echo $comparison
rm -r $ROOT/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.baseline-comp/
perl ~/checkouts/team31/scripts/compare_sig_lists_baseline.pl \
--dir $ROOT/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.baseline-comp \
$ROOT/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.tsv \
$ROOT/deseq2-baseline-grandhet-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.tsv \
$ROOT/deseq2-baseline-grandhet-blacklist-adj-gt-adj-sex-stage-nicole-definite-maybe-outliers/$comparison.tsv
done
done
```

## KO lines summary

Get stage information for each Line.
data/dmdd-genes.txt has the mapping of directory names to gene names and Ensembl IDs

log2 Fold Change for each gene in homs (vs het_wt) and hets (vs wt) is in data/ko_expr.tsv

Numbers of significantly differentially expressed genes in each line is in data/sig_gene_counts.tsv

Get data from OMIM for mutant genes
```
Rscript get_mim_data.R data/dmdd-genes.txt
# produces output/human-mim.tsv
# Edit file to add an extra column called mim_short with the shortest abbreviation from description
```

Run R script to generate figure
```
Rscript fig4.R \
data/dmdd-genes.txt \
output/sample_info.txt \
output/KOs_ordered_by_delay.txt \
data/ko_expr.tsv \
data/Mm_GRCm38_e88_baseline.rda \
data/sig_gene_counts.tsv \
output/human-mim-edited.tsv

```
