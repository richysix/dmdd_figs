# Fig.4

## KO lines summary

data/dmdd-genes.txt has the mapping of directory names to gene names and Ensembl IDs

log2 Fold Change for each gene in homs (vs het_wt) and hets (vs wt) is in data/ko_expr.tsv

Numbers of significantly differentially expressed genes in each line is in data/sig_gene_counts.tsv

Get data from OMIM for mutant genes
```
Rscript get_mim_data.R data/dmdd-genes.txt
# produces output/human-mim.tsv
# Edit file to add an extra column called mim_short with the shortest abbreviation from description
awk -F"\t" '{if(NR == 1){ print $0 "\tmim_short" }
else{ if($9 == ""){ print $0 "\t" }
else{ array_length = split($9, info, ";"); minDescription = "";
for (i = 1; i <= array_length; i++){
description = info[i]; gsub(" ", "", description); sub("SYNDROME", "", description);
if(description == ""){ continue; }
if(minDescription == ""){ minDescription = description; }
else{ if(length(description) < length(minDescription)){ minDescription = description; } }
}
print $0 "\t" minDescription; } }
}' output/human-mim.tsv > output/human-mim-edited.tsv
```

Run R script to generate figure
```
Rscript fig4.R \
data/dmdd-genes.txt \
output/sample_info.txt \
output/KOs_ordered_by_delay.txt \
data/ko_expr.tsv \
data/Mm_GRCm38_e88_baseline.rds \
data/sig_gene_counts.tsv \
output/human-mim-edited.tsv
```
