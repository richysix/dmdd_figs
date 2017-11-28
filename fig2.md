# Fig.2

```bash
# set up working directory
# change this if you are trying recreate the analysis
# everything else should then be relative to this directory
export ROOT=/lustre/scratch117/maz/team31/projects/mouse_DMDD
```

## KO lines summary

Get stage information for each Line

```bash
for mut in $( cut -f2 $ROOT/lane-process/dmdd/deseq2/samples.txt  | grep _ | \
sed -E 's/_(wt|het|hom)//' | sort -u | grep -v Sh3pxd2a_i )
do
cat $ROOT/lane-process/$mut/deseq2-baseline-grandhet-blacklist-adj-gt-adj-sex-stage-nicole-definite-maybe-outliers/samples.txt
done | grep -vE 'condition|baseline' | \
perl -lane 'BEGIN{print join("\t", "", qw{condition group stage} ); }
{ print $_ }' > $ROOT/samples-minus-outliers.txt
```

Get log2 Fold Change for each gene in homs (vs het_wt) and hets (vs wt)

```bash
# lane-process/dmdd-genes.txt has the mapping of directory names to
# gene names and Ensembl IDs

# get log2fc for hom_vs_het_wt if it exists, then hom_vs_wt, then hom_vs_het
cat /dev/null > KO_expr.tsv
cat /dev/null > KO_expr.err
for mut in $( cut -f2 $ROOT/lane-process/dmdd/deseq2/samples.txt  | grep _ | \
sed -E 's/_(wt|het|hom)//' | sort -u | grep -v Sh3pxd2a_i )
do
geneId=$( grep -E "^$mut[[:space:]]" lane-process/dmdd-genes.txt | cut -f4 )
hom=0
for comparison in hom_vs_wt hom_vs_het
do
file="$ROOT/lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.tsv"
if [[ ! -e $file ]]
  then
    echo "$mut $file does not exist" >> KO_expr.err
  else
    grep -E "$geneId" $file | perl -F"\t" -lane '{print join("\t", $F[0], "'$mut'", "'$comparison'", $F[3], ); }'
    hom=1
    break
fi
done >> KO_expr.tsv
if [[ $hom -eq 0 ]]; then
    echo -e "$geneId\t$mut\thom\tNA" >> KO_expr.tsv
fi
for comparison in het_vs_wt
do
file="$ROOT/lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.tsv"
if [[ ! -e $file ]]
  then
    echo "$mut $file does not exist" >> KO_expr.err
    echo -e "$geneId\t$mut\thet\tNA" >> KO_expr.tsv
  else
    grep -E "$geneId" $file | grep ^ENS | perl -F"\t" -lane '{print join("\t", $F[0], "'$mut'", "'$comparison'", $F[3], ); }'
    break
fi
done
done >> KO_expr.tsv

# check numbers
awk '{print $3}' KO_expr.tsv | sort | uniq -c
# should be
# 72 het_vs_wt
#  2 hom_vs_het
# 51 hom_vs_wt
```

