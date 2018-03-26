# Fig.2

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

Get stage information for each Line

```bash
# lane-process/dmdd-genes.txt has the mapping of directory names to
# gene names and Ensembl IDs
$ROOT/lane-process/dmdd-genes.txt

# get individual sample info
for mut in $( cut -f2 $ROOT/lane-process/dmdd/deseq2/samples.txt  | grep _ | \
sed -E 's/_(wt|het|hom)//' | sort -u | grep -v Sh3pxd2a_i )
do
cat $ROOT/lane-process/$mut/deseq2-baseline-grandhet-blacklist-adj-gt-adj-sex-stage-nicole-definite-maybe-outliers/samples.txt
done | grep -vE 'condition|baseline' | \
perl -F"\t" -lane 'BEGIN{print join("\t", "", qw{condition group stage somite_number} ); }
{ $somite_num = $F[3]; $somite_num =~ s/somite[s]*//xms;
print join("\t", @F, $somite_num, ); }' > $ROOT/samples-minus-outliers.txt
```

Get log2 Fold Change for each gene in homs (vs het_wt) and hets (vs wt)

```bash
# get log2fc for hom_vs_het_wt if it exists, then hom_vs_wt, then hom_vs_het
mkdir ko_expr
cat /dev/null > ko_expr/ko_expr.tsv
cat /dev/null > ko_expr/ko_expr.err
ROOT=/lustre/scratch117/maz/team31/projects/mouse_DMDD
for mut in $( cut -f2 $ROOT/lane-process/dmdd/deseq2/samples.txt  | grep _ | \
sed -E 's/_(wt|het|hom)//' | sort -u | grep -vE 'Sh3pxd2a_i|Cenpl' )
do
geneId=$( grep -E "^$mut[[:space:]]" lane-process/dmdd-genes.txt | cut -f4 )
hom=0
counts=0
for comparison in hom_vs_wt hom_vs_het
do
# check for baseline comparison
file="$ROOT/lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.baseline-comp/ko_response.tsv"
if [[ -e $file ]]
  then
    grep -E "^Gene|^$geneId" $file > ko_expr/$geneId.expr.tsv
    counts=1
    grep -E "$geneId" $file | grep ^ENS | \
    perl -F"\t" -lane '{print join("\t", $F[0], "'$mut'", "'$comparison'", $F[2], ); }' \
        >> ko_expr/ko_expr.tsv
    hom=1
    break
  else
    echo "$mut $file does not exist" >> ko_expr/ko_expr.err
fi
# try without baseline
file="$ROOT/lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.tsv"
if [[ -e $file ]]
  then
    grep -E "^Gene|^$geneId" $file > ko_expr/$geneId.expr.tsv
    counts=1
    grep -E "$geneId" $file | grep ^ENS | \
    perl -F"\t" -lane '{print join("\t", $F[0], "'$mut'", "'$comparison'", $F[3], ); }' \
        >> ko_expr/ko_expr.tsv
    hom=1
    break
  else
    echo "$mut $file does not exist" >> ko_expr/ko_expr.err
fi
done
if [[ $hom -eq 0 ]]; then
    echo -e "$geneId\t$mut\thom\tNA" >> ko_expr/ko_expr.tsv
fi
for comparison in het_vs_wt
do
# baseline
file="$ROOT/lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.baseline-comp/ko_response.tsv"
if [[ -e $file ]]
  then
    if [[ $counts -eq 0 ]]; then
      grep -E "^Gene|^$geneId" $file > ko_expr/$geneId.expr.tsv
    fi
    grep -E "$geneId" $file | grep ^ENS | \
    perl -F"\t" -lane '{print join("\t", $F[0], "'$mut'", "'$comparison'", $F[2], ); }' \
        >> ko_expr/ko_expr.tsv
    break
  else
    echo "$mut $file does not exist" >> ko_expr/ko_expr.err
fi
file="$ROOT/lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.tsv"
if [[ -e $file ]]
  then
    if [[ $counts -eq 0 ]]; then
      grep -E "^Gene|^$geneId" $file > ko_expr/$geneId.expr.tsv
    fi
    grep -E "$geneId" $file | grep ^ENS | \
    perl -F"\t" -lane '{print join("\t", $F[0], "'$mut'", "'$comparison'", $F[3], ); }' \
        >> ko_expr/ko_expr.tsv
    break
  else
    echo "$mut $file does not exist" >> ko_expr/ko_expr.err
    echo -e "$geneId\t$mut\thet\tNA" >> ko_expr/ko_expr.tsv
fi
done
done

# check numbers
awk '{print $3}' ko_expr/ko_expr.tsv | sort | uniq -c
      2 het
     71 het_vs_wt
     20 hom
      2 hom_vs_het
     51 hom_vs_wt
```

Get numbers of significantly differentially expressed genes

```
perl -le 'print join("\t", qw{Gene Comparison Set Type Count});' \
 > data/sig_gene_counts.tsv
for comparison in hom_vs_het_wt het_vs_wt
do
cut -f2 $ROOT/lane-process/dmdd/deseq2/samples.txt  | grep _ | \
sed -E 's/_(wt|het|hom)//' | sort -u | grep -vE 'Sh3pxd2a_i|Cenpl' | \
perl ./get_number_sig_genes.pl \
--dir $ROOT/lane-process \
--comparison $comparison
done >> data/sig_gene_counts.tsv
echo 'Ift140
Oaz1' | \
perl ./get_number_sig_genes.pl \
--dir $ROOT/lane-process \
--comparison hom_vs_het >> data/sig_gene_counts.tsv

```
