Add repeat Groups to data file

```
mv data/fig5d_repeats_de.tsv data/fig5d_repeats_de.tmp
perl -F"\t" -lane 'if($. == 53){$group = "Unk"}
if($. < 53){$group = "SAT"}
if($. < 52){$group = "DNA"}
if($. < 50){$group = "LTR"}
if($. < 33){$group = "Line"}
if($. < 28){$group = "Unk"}
if($. < 27){$group = "SAT"}
if($. < 26){$group = "DNA"}
if($. < 24){$group = "LTR"}
if($. < 7 ){$group = "Line"}
if($. < 2){$group = "Group"}
print join("\t", $group, @F[0..3]);' data/fig5d_repeats_de.tmp > data/fig5d_repeats_de.tsv
rm data/fig5d_repeats_de.tmp
```

Download Dhx35/Morc2a data

```
curl -LO https://ndownloader.figshare.com/files/11865374
mkdir data/notranscriptome-repeats
mv 11865374 data/notranscriptome-repeats/notranscriptome-repeats_data.tgz
cd data/notranscriptome-repeats
tar -xzvf notranscriptome-repeats_data.tgz
cd ../../
```

Make sig files for heatmap

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