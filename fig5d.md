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

