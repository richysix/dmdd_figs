# Fig.4

```bash
# set up working directory
# change this if you are trying recreate the analysis
# everything else should then be relative to this directory
export ROOT=/lustre/scratch117/maz/team31/projects/mouse_DMDD
```

Get GO results

```
for domain in BP MF CC
do
cut -f2 $ROOT/lane-process/dmdd/deseq2/samples.txt  | grep _ | \
sed -E 's/_(wt|het|hom)//' | sort -u | grep -vE 'Sh3pxd2a_i|Cenpl' | \
perl get_go_results.pl \
--dir $ROOT/lane-process --go_domain $domain \
--baseline_file data/baseline_go_results.$domain.tsv \
--shared_file data/topgo-shared-terms.$domain.txt \
--observed_threshold 2 --fe_threshold 1
done > data/go_results.tmp

# remove header lines
head -n1 data/go_results.tmp > data/go_results.tsv
grep -v '^Gene' data/go_results.tmp >> data/go_results.tsv
rm data/go_results.tmp

# data/go_results.tsv is one of the input files for fig4.R
```



## EMAPA

Collect together all EMAPA results

```
# collect data for all mutants for each set
# make config file and run with collate_emap_results
cat /dev/null > $ROOT/mouse_dmdd_figs/data/sets_to_use.tsv
echo -e "Mutant\tComparison\tSet\tFile" >> $ROOT/mouse_dmdd_figs/data/sets_to_use.tsv
base='$ROOT/lane-process'
for set in ko_response ko_response mrna_as_wt not_used
do
for mut in $( < output/KOs_delayed.txt ) # delayed
do
# figure out comparison from samples file
num_wts=$( grep -v condition $base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/samples.txt | \
cut -f2 | grep wt | wc -l )
num_homs=$( grep -v condition $base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/samples.txt | \
cut -f2 | grep hom | wc -l )
if [[ $num_wts -eq 0 ]]; then
  comparison='hom_vs_het'
elif [[ $num_homs -eq 0 ]]; then
  comparison='het_vs_wt'
else
  comparison='hom_vs_het_wt'
fi
file=$base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.baseline-comp/$set.emap/results.tsv
echo -e "$mut\t$comparison\t$set\t$file" >> $ROOT/mouse_dmdd_figs/data/sets_to_use.tsv
done
done
# rest of mutants
for set in ko_response
do
for mut in $( cut -f2 $base/dmdd/deseq2/samples.txt  | grep _ | \
sed -E 's/_(wt|het|hom)//' | sort -u | grep -vE 'Sh3pxd2a_i|Cenpl' | \
grep -vf output/KOs_delayed.txt )
do
# figure out comparison from samples file
num_wts=$( grep -v condition $base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/samples.txt | \
cut -f2 | sort -u | grep wt | wc -l )
num_homs=$( grep -v condition $base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/samples.txt | \
cut -f2 | sort -u | grep hom | wc -l ) 
if [[ $num_wts -eq 0 ]]; then
  comparison='hom_vs_het'
elif [[ $num_homs -eq 0 ]]; then
  comparison='het_vs_wt'
else
  comparison='hom_vs_het_wt'
fi
file=$base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/$comparison.emap/results.tsv
echo -e "$mut\t$comparison\t$set\t$file" >> $ROOT/mouse_dmdd_figs/data/sets_to_use.tsv
done
done

# run script to collect together results from appropriate comparisons
echo -e "Gene\tComparison\tSet\tTerm ID\tDescription\tAnnotated\
\tExpected\tObserved\tFold Enrichment\tAdjusted p value\t-log10(pvalue)" \
 > $ROOT/mouse_dmdd_figs/data/emap_results.all.tsv

perl collate_emap_results.pl --header \
$ROOT/mouse_dmdd_figs/data/sets_to_use.tsv \
 >> $ROOT/mouse_dmdd_figs/data/emap_results.all.tsv

# get just mrna_abnormal results
head -n1 data/emap_results.all.tsv > output/emap_results.mrna_abnormal.tsv
grep mrna_abnormal data/emap_results.all.tsv | sort -t$'\t' -k5,5 >> output/emap_results.mrna_abnormal.tsv

 
export R_LIBS_USER='/lustre/scratch117/maz/team31/projects/mouse_DMDD/mouse_dmdd_figs/.R/lib:/software/team31/R-3.3.0'
/software/R-3.3.0/bin/Rscript ~rw4/checkouts/team31/scripts/bubble_plot.R \
--x_var='Gene' --y_var='Term ID' --y_labels 'Description' \
--x_levels=Edc4,Bap1,Ift140,Psph,Dhx35,Gpr107,Ehmt1,Zfyve20,Brd2,Nek9,Supt3,Col4a3bp,Hira,Fcho2,Fam160a,Kif3b,Pigf,Morc2a,Dennd4c,Fryl,Mad2l2,H13,Ccdc6,Nadk2 \
--size_var='Fold Enrichment' --fill_var='-log10(pvalue)' \
--reverse_y \
--width=12 --height=100 \
--output_file plots/emap_results.mrna_abnormal.pdf \
--output_data_file output/emap_results.mrna_abnormal.rda \
output/emap_results.mrna_abnormal.tsv

```




process EMAPA terms to collapse

```
/software/R-3.3.0/bin/Rscript process_emap.R data/emap_results.all.tsv \
/lustre/scratch117/maz/team31/projects/mouse_DMDD/emap/EMAPA_Nov_17_2017.obo \
data/root_terms.txt

tr '\r' '\n' < ~/sanger/ZMP/mouse\ DMDD/output/duplicated_terms.txt \
 > ~/sanger/ZMP/mouse\ DMDD/output/duplicated_terms-edited.tsv
rm ~/sanger/ZMP/mouse\ DMDD/output/duplicated_terms.txt
scp ~/sanger/ZMP/mouse\ DMDD/output/duplicated_terms-edited.tsv \
gen1:/lustre/scratch117/maz/team31/projects/mouse_DMDD/output/duplicated_terms-edited.tsv
```

Calculate the pairwise gene overlap from the ko_response lists for all delayed mutants
```
# create sig_genes files
base=/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process
cat /dev/null > output/sig_genes_files.txt
# delayed mutants
for mut in $(grep -v None output/KOs_ordered_by_delay.txt | cut -f1)
do
sig_file=$base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/ko_response.sig.tsv
if [[ -e $sig_file ]]
then
  sig_genes_file=$(echo $sig_file | sed -e 's|.tsv|_genes.txt|' )
  cut -f1 $sig_file > $sig_genes_file
else
  echo $sig_file does not exist
  sig_file=$base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/het_vs_wt.baseline-comp/ko_response.sig.tsv
  if [[ -e $sig_file ]]
  then
    sig_genes_file=$(echo $sig_file | sed -e 's|.tsv|_genes.txt|' )
    cut -f1 $sig_file > $sig_genes_file
  else
    echo $sig_file does not exist
    sig_file=$base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het.baseline-comp/ko_response.sig.tsv
    if [[ -e $sig_file ]]
    then
      sig_genes_file=$(echo $sig_file | sed -e 's|.tsv|_genes.txt|' )
      cut -f1 $sig_file > $sig_genes_file
    else
      echo $sig_file does not exist
    fi
  fi
fi
echo -e "$mut\t$sig_genes_file" >> output/sig_genes_files.txt
done
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Ift140/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/ko_response.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Ift140/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/het_vs_wt.baseline-comp/ko_response.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Coq4/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/ko_response.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Smg1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/ko_response.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Copg/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/ko_response.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Kif18b/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/ko_response.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Oaz1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/ko_response.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Oaz1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/het_vs_wt.baseline-comp/ko_response.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Timmdc1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/ko_response.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Mtor/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/ko_response.sig.tsv does not exist

# not delayed mutants
for mut in $(grep None output/KOs_ordered_by_delay.txt | cut -f1)
do
sig_file=$base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv
if [[ -e $sig_file ]]
then
  sig_genes_file=$(echo $sig_file | sed -e 's|.tsv|_genes.txt|')
  cut -f1 $sig_file > $sig_genes_file
else
  echo $sig_file does not exist
  sig_file=$base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/het_vs_wt.sig.tsv
  if [[ -e $sig_file ]]
  then
    sig_genes_file=$(echo $sig_file | sed -e 's|.tsv|_genes.txt|')
    cut -f1 $sig_file > $sig_genes_file
  else
    echo $sig_file does not exist
  fi
fi
echo -e "$mut\t$sig_genes_file" >> output/sig_genes_files.txt
done
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Gfm1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Dhodh/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Dpm1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Tmem30a/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Dctn4/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Pgap2/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Pigl/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Elac2/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Orc1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Smc3/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Dctn1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Atxn10/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Wrap53/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Crls/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.sig.tsv does not exist

# check grepping for mutant name only produces one file from output/sig_genes_files.txt
for mut in $(cut -f1 output/KOs_ordered_by_delay.txt)
do
echo -n "$mut "
grep -c $mut output/sig_genes_files.txt
done | awk '{if($2 != 1){print $0}}'

# overlap lists
mutants=( $(cut -f1 output/KOs_ordered_by_delay.txt) )
cat /dev/null > output/ko_response.sig_genes.err
base=/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process
for mut1 in $( seq 0 $(( ${#mutants[@]} - 1 )) )
do
for mut2 in $( seq 0 $(( ${#mutants[@]} - 1 )) )
do
if [ "$mut1" -lt "$mut2" ]; then
  file1=$( grep ${mutants[$mut1]} output/sig_genes_files.txt | cut -f2 )
  file2=$( grep ${mutants[$mut2]} output/sig_genes_files.txt | cut -f2 )
  if [[ -e $file1 && -e $file2 ]]; then
    perl ~/checkouts/team31/scripts/intersect_lists.pl --counts --jaccard \
    --header --output_no_header --no_venn \
    $file1 $file2
  else
    echo "ONE OF THE FILES DOES NOT EXIST: $file1 $file2" >> output/ko_response.sig_genes.err
  fi
fi
done
done > output/ko_response-hom_vs_het_wt-sig_genes.out.tmp

# edit output file
sed -e 's|/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/||g' \
output/ko_response-hom_vs_het_wt-sig_genes.out.tmp | \
sed -e 's|/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/[a-z._-/]*.txt||g' \
 > output/ko_response-hom_vs_het_wt-sig_genes.out

# plot heatmap
/software/R-3.3.0/bin/Rscript cluster_by_overlap.R \
--output_base plots/ko_response-gene_overlap \
--cluster_methods ward.D2 \
--output_data_file output/ko_response-gene_overlap.rda \
output/ko_response-hom_vs_het_wt-sig_genes.out

```

Calculate the pairwise gene overlap from the mrna_abnormal lists for all delayed mutants
```
mutants=( $(sort output/KOs_delayed.txt) )
cat /dev/null > output/mrna_abnormal.sig_genes.err
base=/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process
for mut1 in $( seq 0 $(( ${#mutants[@]} - 1 )) )
do
for mut2 in $( seq 0 $(( ${#mutants[@]} - 1 )) )
do
if [ "$mut1" -lt "$mut2" ]; then
  file1=$base/${mutants[$mut1]}/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt
  file2=$base/${mutants[$mut2]}/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt
  if [[ -e $file1 && -e $file2 ]]; then
    file1_count=$( cat $file1 | wc -l )
    file2_count=$( cat $file2 | wc -l )
    if [[ $file1_count -gt 1 && $file2_count -gt 1 ]]; then
      perl ~/checkouts/team31/scripts/intersect_lists.pl --counts --jaccard \
        --header --output_no_header --no_venn \
        $file1 $file2
    else
      echo "ONE OF THE FILES HAS NO RESULTS: ${mutants[$mut1]}-hom_vs_het_wt ${mutants[$mut2]}-hom_vs_het_wt"  >> output/mrna_abnormal.sig_genes.err
    fi
  else
    echo "ONE OF THE FILES DOES NOT EXIST: ${mutants[$mut1]}-hom_vs_het_wt ${mutants[$mut2]}-hom_vs_het_wt" >> output/mrna_abnormal.sig_genes.err
  fi
fi
done
done | sed -e 's|/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/||g' | \
sed -e 's|/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt||g' \
 > output/mrna_abnormal-hom_vs_het_wt-sig_genes.out

# cluster all based on gene list overlaps
/software/R-3.3.0/bin/Rscript cluster_by_overlap.R \
--output_base plots/mrna_abnormal-PC3-jaccard-all \
--cluster_methods ward.D2 \
--output_data_file output/mrna_abnormal-jaccard-all.rda \
output/mrna_abnormal-hom_vs_het_wt-sig_genes.out
```

get gene lists for specific overlaps
```
# get genes in each set
base=/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process
for mut in $( < output/KOs_delayed.txt )
do
file=$base/${mut}/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt
if [[ -e $file ]]; then
  lines=$(wc -l $file | awk '{print $1}')
  if [[ $lines -eq 1 ]]; then
    echo "$file: NO SIG GENES" 1>&2
  else
    grep ENSMUSG $file | \
    perl -F"\t" -lane 'print join("\t", $F[0], "'$mut'", "1");'
  fi
else
  echo "$file: DOES NOT EXIST" 1>&2
fi
done > output/mrna_abnormal-sig_genes-upset.tmp
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Ift140/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: DOES NOT EXIST
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Coq4/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: DOES NOT EXIST
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Pdzk1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: NO SIG GENES
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Smg1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: DOES NOT EXIST
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Zkscan14/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: NO SIG GENES
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Copg/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: DOES NOT EXIST
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Kif18b/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: DOES NOT EXIST
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Sh3pxd2a/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: NO SIG GENES
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Lztr/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: NO SIG GENES
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Oaz1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: DOES NOT EXIST
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Otud7b/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: NO SIG GENES
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Timmdc1/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: DOES NOT EXIST
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Mtor/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: DOES NOT EXIST
/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Smg9/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.sig_genes.txt: NO SIG GENES

# get gene lists from Ehmt1 Supt3 Zfyve20 Mad2l2 and Bap1 Hira intersection
grep -E 'Ehmt1|Supt3|Zfyve20|Mad2l2' output/mrna_abnormal-sig_genes-upset.tmp > output/mrna_abnormal-sig_genes-EMSZ-upset.tmp
grep -E 'Bap1|Hira' output/mrna_abnormal-sig_genes-upset.tmp > output/mrna_abnormal-sig_genes-BH-upset.tmp

# script to change the data from long to wide
echo 'library("reshape2")
Args <- commandArgs()
data <- read.delim(Args[6], header = FALSE)
names(data) <- c("gene", "mut", "present")
data_matrix <- dcast(data, gene ~ mut, value.var = "present", fill = 0)
write.table(data_matrix, file = Args[7], sep = "\t", quote = FALSE, row.names = FALSE)' > reshape-long_to_wide.R

for set in BH EMSZ
do
/software/R-3.3.0/bin/Rscript reshape-long_to_wide.R output/mrna_abnormal-sig_genes-$set-upset.tmp output/mrna_abnormal-sig_genes-$set-upset.tsv
done

# script to get list of genes in the intersection
echo 'Args <- commandArgs()
upset_wide <- read.delim(Args[6])
intersect_all <- upset_wide[ rowSums(upset_wide[,seq(2,length(upset_wide),1)]) == length(upset_wide) - 1, ]
write.table(as.data.frame(intersect_all$gene), file = Args[7], quote = FALSE,
row.names = FALSE, col.names = FALSE)' > get_intersection_genes.R

for set in BH EMSZ
do
/software/R-3.3.0/bin/Rscript get_intersection_genes.R \
output/mrna_abnormal-sig_genes-$set-upset.tsv \
output/mrna_abnormal-sig_genes-$set-intersection.txt
done
```

GO Enrichment

```
# make file of all genes that were tested in mrna_abnormal
for mut in $( < output/KOs_delayed.txt ) # delayed
do
file=lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.tsv
if [[ ! -e $file ]]; then
  file=lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/het_vs_wt.baseline-comp/mrna_abnormal.tsv
  if [[ ! -e $file ]]; then
	file=lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het.baseline-comp/mrna_abnormal.tsv
	if [[ ! -e $file ]]; then
      echo 1>&2 "$file DOES NOT EXIST"
	else
      cat $file
	fi
  else
    cat $file
  fi
else
  cat $file
fi
done | grep -v '^Gene' | perl -F"\t" -lane 'if($F[1] ne "NA"){ print $F[0] }' | \
sort -u > output/mrna_abnormal-all_genes.txt

# get genes from an all.tsv file for topgo
head -n1 output/all_samples_merged.counts.tsv > output/mrna_abnormal-geneU.tsv
sort -t$'\t' -k1,1 output/all_samples_merged.counts.tsv | \
join -t$'\t' - output/mrna_abnormal-all_genes.txt  >> output/mrna_abnormal-geneU.tsv

# run topgo
species=mus_musculus
ensembl=88
export R_LIBS_USER=/software/team31/R:.R/lib
for set in BH EMSZ
do
rm output/mrna_abnormal-$set-intersection.fisher/topgo.[oe]
bsub -o output/mrna_abnormal-$set-intersection.fisher/topgo.o \
-e output/mrna_abnormal-$set-intersection.fisher/topgo.e \
-R'select[mem>2000] rusage[mem=2000]' -M2000 \
perl -I/software/team31/packages/topgo-wrapper/lib \
/software/team31/packages/topgo-wrapper/script/run_topgo.pl \
--dir output/mrna_abnormal-$set-intersection.fisher \
--input_file output/mrna_abnormal-geneU.tsv \
--genes_of_interest_file output/mrna_abnormal-sig_genes-$set-intersection.txt \
--gene_field 1 \
--name_field 7 \
--description_field 8 \
--go_terms_file /software/team31/packages/topgo-wrapper/data/${species}_e${ensembl}_go.txt \
--header 1 \
--r_binary /software/R-3.1.2/bin/R
done

# aggregate terms using GO graph
domains=(BP MF CC)
ontologies=(biological_process molecular_function cellular_component)
for set in BH EMSZ
do
for idx in $(seq 0 2)
do
for level in 4
do
/software/R-3.1.2/bin/Rscript \
~rw4/checkouts/team31/scripts/topgo_getsubtrees.R \
--directory output --output_base mrna_abnormal-$set-intersection.fisher/${domains[${idx}]}-${level}-collapsed \
--genes_of_interest_file mrna_abnormal-sig_genes-BH-intersection.txt \
--ontology ${domains[${idx}]} --ontology_level $level --gene_id_column 'Gene.ID' \
--term_column GO.ID --term_name_column Term \
mrna_abnormal-$set-intersection.fisher/${ontologies[${idx}]}_all.gene2go.txt \
mrna_abnormal-geneU.tsv \
mrna_abnormal-$set-intersection.fisher/${domains[${idx}]}.sig.tsv
done
done
done

# some terms appear in more than one subtree
for set in BH EMSZ
do
echo -e -n "$set\t"
cut -f1 output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed.sig.tsv | sort | uniq -c | awk '{if($1 > 1){print $0}}' | wc -l
done
BH      193
EMSZ    129

# remove genes column and manually annotate which root term to keep
# get just those terms that appear more than once
for set in BH EMSZ
do
cut -f1 output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed.sig.tsv | \
sort | uniq -c | awk '{if($1 > 1){print $2}}' | xargs -I term \
grep term output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed.sig.tsv | \
cut -f1-10 > output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed.sig.tmp
done

# check genes have been filtered properly
for set in BH EMSZ
do
echo -e -n "$set\t"
grep -v remove output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed.sig.tmp | cut -f1 | sort -u | wc -l
done
BH      193
EMSZ    128

# get terms that appear only once and cat together with manually filtered terms
for set in BH EMSZ
do
cut -f1 output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed.sig.tsv | \
sort | uniq -c | awk '{if($1 > 1){print $2}}' | \
sort > output/mrna_abnormal-$set-intersection.fisher/BP-4-duplicate-terms.txt
# get non-duplicate terms
( sort -t$'\t' -k1,1 output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed.sig.tsv | \
join -t$'\t' -v 1 - output/mrna_abnormal-$set-intersection.fisher/BP-4-duplicate-terms.txt | cut -f1-10 ;
grep -v remove output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed.sig.tmp | cut -f1-10) | \
sort -t$'\t' -gk10,10 |  perl -F"\t" -lane 'if($. == 1){print join("\t", @F, qw{Fold.Enrichment log10.pvalue}); }
else{ print join("\t", @F, $F[7]/$F[8], -(log($F[9])/log(10)) ); }' > output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed-filtered.sig.tsv
done

# bubble plot of terms against umbrella terms
for set in BH EMSZ
do
/software/R-3.3.0/bin/Rscript bubble_plot.R --x_var='Root.Term' --y_var='GO.ID' \
--size_var='Fold.Enrichment' --fill_var='log10.pvalue' \
--y_labels=Term --x_labels=Root.Term.Name --reverse_y \
--output_file output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed-filtered-bubble.pdf \
--output_type=pdf --height 25 \
output/mrna_abnormal-$set-intersection.fisher/BP-4-collapsed-filtered.sig.tsv
done

# do Fcho2 as well
domains=(BP MF CC)
ontologies=(biological_process molecular_function cellular_component)
level=4
WD=/lustre/scratch117/maz/team31/projects/mouse_DMDD/lane-process/Fcho2/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp
for idx in $(seq 0 2)
do
/software/R-3.3.0/bin/Rscript \
~rw4/checkouts/team31/scripts/topgo_getsubtrees.R \
--directory $WD --output_base mrna_abnormal.fisher/${domains[${idx}]}-${level}-collapsed \
--ontology ${domains[${idx}]} --ontology_level $level \
--gene_id_column 'Gene.ID' --pvalue_column adjp \
--term_column GO.ID --term_name_column Term \
$WD/mrna_abnormal.fisher/${ontologies[${idx}]}_all.gene2go.txt \
$WD/mrna_abnormal.tsv \
$WD/mrna_abnormal.fisher/${domains[${idx}]}.sig.tsv
done

for domain in BP MF CC
do
echo -e -n "$domain\t"
cut -f1 $WD/mrna_abnormal.fisher/$domain-4-collapsed.sig.tsv | sort | uniq -c | awk '{if($1 > 1){print $0}}' | wc -l
done
BP      156
MF      3
CC      0

# remove genes column and manually annotate which root term to keep
# get just those terms that appear more than once
for domain in BP MF CC
do
cut -f1 $WD/mrna_abnormal.fisher/$domain-4-collapsed.sig.tsv | \
sort | uniq -c | awk '{if($1 > 1){print $2}}' | xargs -I term \
grep term $WD/mrna_abnormal.fisher/$domain-4-collapsed.sig.tsv | \
cut -f1-10 > $WD/mrna_abnormal.fisher/$domain-4-collapsed.sig.tmp
done

# check genes have been filtered properly
for domain in BP MF CC
do
echo -e -n "$domain\t"
grep -v remove $WD/mrna_abnormal.fisher/$domain-4-collapsed.sig.tmp | cut -f1 | sort -u | wc -l
done
BP      156
MF      3
CC      0

# get terms that appear only once and cat together with manually filtered terms
for domain in BP MF CC
do
cut -f1 $WD/mrna_abnormal.fisher/$domain-4-collapsed.sig.tsv | \
sort | uniq -c | awk '{if($1 > 1){print $2}}' | \
sort > $WD/mrna_abnormal.fisher/$domain-4-duplicate-terms.txt
# get non-duplicate terms
( sort -t$'\t' -k1,1 $WD/mrna_abnormal.fisher/$domain-4-collapsed.sig.tsv | \
join -t$'\t' -v 1 - $WD/mrna_abnormal.fisher/$domain-4-duplicate-terms.txt | cut -f1-10 ;
grep -v remove $WD/mrna_abnormal.fisher/$domain-4-collapsed.sig.tmp | cut -f1-10) | \
sort -t$'\t' -gk10,10 |  perl -F"\t" -lane 'if($. == 1){print join("\t", @F, qw{Fold.Enrichment log10.pvalue}); }
else{ print join("\t", @F, $F[7]/$F[8], -(log($F[9])/log(10)) ); }' > $WD/mrna_abnormal.fisher/$domain-4-collapsed-filtered.sig.tsv
done

# collect together results
level=4
for domain in BP
do
echo -e "Set\tTerm.ID\tDescription\tCount\tmaxlog10p" > output/mrna_abnormal-F_BH_EMSZ-fisher-level${level}-results.tsv
for set in BH EMSZ
do
grep -v '^GO.ID' output/mrna_abnormal-$set-intersection.fisher/$domain-$level-collapsed-filtered.sig.tsv | \
perl -F"\t" -lane 'BEGIN{ %terms = ();} {
if(!exists $terms{$F[1]}){ $terms{$F[1]} = {Description => $F[2], Count => 1, maxlog10p => $F[11]}; }
else{ $terms{$F[1]}->{'Count'}++;
$terms{$F[1]}->{'maxlog10p'} = $F[11] > $terms{$F[1]}->{'maxlog10p'} ? $F[11] : $terms{$F[1]}->{'maxlog10p'}; }
} END{ foreach $term (sort keys %terms ){
print join("\t", "'$set'", $term, $terms{$term}->{'Description'},
$terms{$term}->{'Count'}, $terms{$term}->{'maxlog10p'}, )} }' | sort -t$'\t' -grk4,4
done
set=F
grep -v '^GO.ID' $WD/mrna_abnormal.fisher/$domain-$level-collapsed-filtered.sig.tsv | \
perl -F"\t" -lane 'BEGIN{ %terms = ();} {
if(!exists $terms{$F[1]}){ $terms{$F[1]} = {Description => $F[2], Count => 1, maxlog10p => $F[11]}; }
else{ $terms{$F[1]}->{'Count'}++;
$terms{$F[1]}->{'maxlog10p'} = $F[11] > $terms{$F[1]}->{'maxlog10p'} ? $F[11] : $terms{$F[1]}->{'maxlog10p'}; }
} END{ foreach $term (sort keys %terms ){
print join("\t", "'$set'", $term, $terms{$term}->{'Description'},
$terms{$term}->{'Count'}, $terms{$term}->{'maxlog10p'}, )} }' | sort -t$'\t' -grk4,4
done >> output/mrna_abnormal-F_BH_EMSZ-fisher-level${level}-results.tsv

# filter less well represented terms
perl -F"\t" -lane 'if($F[3] > 3 || $. == 1){ print $_; }' \
output/mrna_abnormal-F_BH_EMSZ-fisher-level4-results.tsv \
 > output/mrna_abnormal-F_BH_EMSZ-fisher-level4-results-filtered.tsv

for domain in BP
do
for type in pdf svg
do
/software/R-3.3.0/bin/Rscript ~rw4/checkouts/team31/scripts/bubble_plot.R \
--x_var='Set' --y_var='Term.ID' \
--size_var='Count' --fill_var='maxlog10p' \
--y_labels=Description --reverse_y \
--output_file plots/GO-$domain-F_vs_BH_vs_EMSZ-level4-collapsed.$type \
--output_type=$type --width 8 \
output/mrna_abnormal-F_BH_EMSZ-fisher-level4-results-filtered.tsv
done
done

```

EMAPA

```
# Have a sig genes list
# Need an all list of just gene ids
cut -f1 output/mrna_abnormal-geneU.tsv > output/mrna_abnormal-geneU-gene_ids_only.txt

for set in BH EMSZ
do
mkdir output/mrna_abnormal-$set-intersection.emap
bsub -o output/mrna_abnormal-$set-intersection.emap/emap.o \
-e output/mrna_abnormal-$set-intersection.emap/emap.e \
-M1100 -R'select[mem>1100] rusage[mem=1100]' \
"perl ~/checkouts/team31/scripts/emapa_enrichment.pl \
--java /software/java/bin/java --Xms 1000m --Xmx 1000m \
--ontologizer ~rw4/bin/Ontologizer.jar \
--emapa_annotations /lustre/scratch117/maz/team31/projects/mouse_DMDD/emap/EMAPA-annotations.txt \
--obo_file /lustre/scratch117/maz/team31/projects/mouse_DMDD/emap/EMAPA_Nov_17_2017-edited.obo \
--input_file output/mrna_abnormal-geneU-gene_ids_only.txt \
--enriched_genes output/mrna_abnormal-sig_genes-$set-intersection.txt \
--dir output/mrna_abnormal-$set-intersection.emap \
--mtc Benjamini-Hochberg --calculation Parent-Child-Union \
--de_sig_level 0.05 --emap_sig_level 0.01 \
--gene_id_field 1 \
--comparison mrna_abnormal-$set-intersection"
done

# collapse terms
for set in BH EMSZ
do
for level in $(seq 3 4)
do
/software/R-3.3.0/bin/Rscript \
~/checkouts/team31/scripts/collapse_ontology_enrichment_results.R \
--output_directory output --output_base mrna_abnormal-$set-intersection.emap/subtrees-level${level} \
--term_column Term.ID --term_name_column Description --namespace anatomical_structure \
--ontology_level $level \
$ROOT/emap/EMAPA_Nov_17_2017.obo \
output/mrna_abnormal-$set-intersection.emap/results.tsv
done
done

# run on Fcho2 delay EMAPA enrichments as well
for level in $(seq 3 4)
do
/software/R-3.3.0/bin/Rscript \
~/checkouts/team31/scripts/collapse_ontology_enrichment_results.R \
--output_directory $ROOT/lane-process/Fcho2/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.emap/ \
--output_base results-collapsed-level${level} \
--term_column Term.ID --term_name_column Description --namespace anatomical_structure \
--ontology_level $level \
$ROOT/emap/EMAPA_Nov_17_2017.obo \
$ROOT/lane-process/Fcho2/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.emap/results.tsv
done
```

```
# plot terms vs collapsed term name as a bubble plot
set=BH
for level in $(seq 3 4)
do
/software/R-3.3.0/bin/Rscript bubble_plot.R --x_var='subtree_root_id' --y_var='Term.ID' \
--size_var='Fold.Enrichment' --fill_var='log10.pvalue' \
--y_labels=Description --x_labels=subtree_root_name --reverse_y \
--output_file output/mrna_abnormal-$set-intersection.emap/subtrees-level${level}-bubble.svg \
--output_type=svg --height 15 \
output/mrna_abnormal-$set-intersection.emap/subtrees-level${level}-results.tsv
done

set=EMSZ
for level in $(seq 3 4)
do
/software/R-3.3.0/bin/Rscript bubble_plot.R --x_var='subtree_root_id' --y_var='Term.ID' \
--size_var='Fold.Enrichment' --fill_var='log10.pvalue' \
--y_labels=Description --x_labels=subtree_root_name --reverse_y \
--output_file output/mrna_abnormal-$set-intersection.emap/subtrees-level${level}-bubble.svg \
--output_type=svg --height 25 \
output/mrna_abnormal-$set-intersection.emap/subtrees-level${level}-results.tsv
done

# Fcho2
EMAPA_DIR=$ROOT/lane-process/Fcho2/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.emap
for level in $(seq 3 4)
do
/software/R-3.3.0/bin/Rscript bubble_plot.R --x_var='subtree_root_id' --y_var='Term.ID' \
--size_var='Fold.Enrichment' --fill_var='log10.pvalue' \
--y_labels=Description --x_labels=subtree_root_name --reverse_y \
--output_file $EMAPA_DIR/results-collapsed-level${level}-bubble.svg \
--output_type=svg --height 25 \
$EMAPA_DIR/results-collapsed-level${level}-results.tsv
done
```

```
# cat together results and plot bubble plot
for level in $(seq 3 4)
do
echo -e "Set\tTerm.ID\tDescription\tCount\tmaxlog10p" > output/mrna_abnormal-BH_EMSZ-level${level}-results.tsv
for set in BH EMSZ
do
grep -v '^Term.ID' output/mrna_abnormal-$set-intersection.emap/subtrees-level${level}-results.tsv | \
perl -F"\t" -lane 'BEGIN{ %terms = ();} {
if(!exists $terms{$F[7]}){ $terms{$F[7]} = {Description => $F[8], Count => 1, maxlog10p => $F[10]}; }
else{ $terms{$F[7]}->{'Count'}++;
$terms{$F[7]}->{'maxlog10p'} = $F[10] > $terms{$F[7]}->{'maxlog10p'} ? $F[10] : $terms{$F[7]}->{'maxlog10p'}; }
} END{ foreach $term (sort keys %terms ){
print join("\t", "'$set'", $term, $terms{$term}->{'Description'},
$terms{$term}->{'Count'}, $terms{$term}->{'maxlog10p'}, )} }' | \
sort -t$'\t' -grk4,4 >> output/mrna_abnormal-BH_EMSZ-level${level}-results.tsv
done
done

export R_LIBS_USER='.R/lib'
for level in $(seq 3 4)
do
/software/R-3.3.0/bin/Rscript ~rw4/checkouts/team31/scripts/bubble_plot.R \
--x_var='Set' --y_var='Term.ID' \
--size_var='Count' --fill_var='maxlog10p' \
--y_labels=Description --reverse_y \
--output_file plots/EMAPA-BH_vs_EMSZ-level${level}-collapsed.svg \
--output_type=svg \
output/mrna_abnormal-BH_EMSZ-level${level}-results.tsv
done

# add in Fcho2
EMAPA_DIR=$ROOT/lane-process/Fcho2/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.baseline-comp/mrna_abnormal.emap
set=F
for level in $(seq 3 4)
do
grep -v '^Term.ID' $EMAPA_DIR/results-collapsed-level${level}-results.tsv | \
perl -F"\t" -lane 'BEGIN{ %terms = ();} {
if(!exists $terms{$F[7]}){ $terms{$F[7]} = {Description => $F[8], Count => 1, maxlog10p => $F[10]}; }
else{ $terms{$F[7]}->{'Count'}++;
$terms{$F[7]}->{'maxlog10p'} = $F[10] > $terms{$F[7]}->{'maxlog10p'} ? $F[10] : $terms{$F[7]}->{'maxlog10p'}; }
} END{ foreach $term (sort keys %terms ){
print join("\t", "'$set'", $term, $terms{$term}->{'Description'},
$terms{$term}->{'Count'}, $terms{$term}->{'maxlog10p'}, )} }' | \
sort -t$'\t' -grk4,4 | cat output/mrna_abnormal-BH_EMSZ-level${level}-results.tsv - \
 > output/mrna_abnormal-F_BH_EMSZ-level${level}-results.tsv
done

for level in $(seq 3 4)
do
/software/R-3.3.0/bin/Rscript ~rw4/checkouts/team31/scripts/bubble_plot.R \
--x_var='Set' --y_var='Term.ID' \
--size_var='Count' --fill_var='maxlog10p' \
--y_labels=Description --reverse_y \
--output_file plots/EMAPA-F_vs_BH_vs_EMSZ-level${level}-collapsed.svg \
--output_type=svg \
output/mrna_abnormal-F_BH_EMSZ-level${level}-results.tsv
done

# filter out terms that appear less than 3 times
perl -F"\t" -lane 'if($F[3] > 2 | $. == 1){ print $_; }' \
output/mrna_abnormal-F_BH_EMSZ-level4-results.tsv \
 > output/mrna_abnormal-F_BH_EMSZ-level4-results-filtered.tsv

for type in pdf svg eps
do
/software/R-3.3.0/bin/Rscript ~rw4/checkouts/team31/scripts/bubble_plot.R \
--x_var='Set' --y_var='Term.ID' \
--size_var='Count' --fill_var='maxlog10p' \
--y_labels=Description --reverse_y \
--output_file plots/EMAPA-F_vs_BH_vs_EMSZ-level4-collapsed.$type \
--output_type=$type --width 7 \
output/mrna_abnormal-F_BH_EMSZ-level4-results-filtered.tsv
done

# make smaller one with no legend
for type in eps
do
/software/R-3.3.0/bin/Rscript ~rw4/checkouts/team31/scripts/bubble_plot.R \
--x_var='Set' --y_var='Term.ID' \
--size_var='Count' --fill_var='maxlog10p' \
--y_labels=Description --reverse_y \
--output_file plots/EMAPA-F_vs_BH_vs_EMSZ-level4-collapsed-no_legend.$type \
--output_type=$type --width 5 \
output/mrna_abnormal-F_BH_EMSZ-level4-results-filtered.tsv
done

# make a version with sets on the y axis and terms on the x
for type in eps
do
/software/R-3.3.0/bin/Rscript ~rw4/checkouts/team31/scripts/bubble_plot.R \
--x_var='Term.ID' --y_var='Set' \
--size_var='Count' --fill_var='maxlog10p' \
--x_labels=Description --reverse_y \
--output_file plots/EMAPA-F_vs_BH_vs_EMSZ-level4-collapsed-rotated.$type \
--output_type=$type --height 5.5 --width 12 \
output/mrna_abnormal-F_BH_EMSZ-level4-results-filtered.tsv
done

```

PCA

```
# first create the list of all files
base=$ROOT/lane-process
for mut in $( grep -v Cenpl output/KOs_ordered_by_delay.txt )
do
file=$base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het_wt.tsv
if [[ ! -e $file ]]; then
  file=$base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/het_vs_wt.tsv
  if [[ ! -e $file ]]; then
	file=$base/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/hom_vs_het.tsv
	if [[ ! -e $file ]]; then
      echo 1>&2 "$file DOES NOT EXIST"
	else
      echo $file
	fi
  else
    echo $file
  fi
else
  echo $file
fi
done | grep -v '^Gene' | \
sort -u > output/all_mutants-counts-files.txt

# make count file
perl ~/checkouts/team31/scripts/merge_deseq_counts.pl \
output/all_mutants-counts-files.txt \
 > output/all_samples_merged.counts.tsv

# make merged samples file
echo -e "\tcondition\tmutant" > $ROOT/mouse_dmdd_figs/output/all_samples-no_Cenpl.tsv
for mut in $( grep -v Cenpl output/KOs_ordered_by_delay.txt | cut -f1 )
do
file=$ROOT/lane-process/$mut/deseq2-blacklist-adj-gt-adj-sex-nicole-definite-maybe-outliers/samples.txt
if [[ ! -e $file ]]; then
  echo 1>&2 "$file DOES NOT EXIST"
else
  cat $file
fi
done | grep -v condition | perl -F"\t" -lane '$mutant = $F[0]; $mutant =~ s/_.* \z//xms;
print join("\t", @F[0,1], $mutant)' \
 >> $ROOT/mouse_dmdd_figs/output/all_samples-no_Cenpl.tsv

# run PCA
for genes in 50000
do
for transform in vst
do
bsub -q basement -o output/pca.ko_response.$genes.$transform.o -e output/pca.ko_response.$genes.$transform.e \
-M20000 -R'select[mem>20000] rusage[mem=20000]' \
"export R_LIBS_USER=/software/team31/R-3.3.0
/software/R-3.3.0/bin/Rscript \
~rw4/checkouts/bio-misc/pca_rnaseq.R \
/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd/output/all_samples_merged.counts.tsv \
/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd/output/all_samples-no_Cenpl.tsv \
/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd/plots/all_mutants-pca.pdf \
/nfs/users/nfs_r/rw4/checkouts/mouse_dmdd/output/all_mutants.$transform.$genes $transform $genes"
done
done

# also make count and samples file with mean count for each mut gt combo
/software/R-3.3.0/bin/Rscript mean_expression_by_mut_gt.R \
output/all_samples_merged.counts.tsv \
$ROOT/samples-minus-outliers.txt \
output/KOs_ordered_by_delay.txt

# run pca script
for genes in 50000
do
for transform in vst
do
bsub -o output/pca.mean_by_mut_gt.$genes.$transform.o \
-e output/pca.mean_by_mut_gt.$genes.$transform.e \
-M5000 -R'select[mem>5000] rusage[mem=5000]' \
"export R_LIBS_USER=/software/team31/R-3.3.0
/software/R-3.3.0/bin/Rscript \
~rw4/checkouts/bio-misc/pca_rnaseq.R \
output/mean_by_mut_gt.counts.tsv \
output/samples_by_mut_by_gt.txt \
plots/mean_by_mut_gt-pca.pdf \
output/mean_by_mut_gt.$transform.$genes $transform $genes 1 0.001"
done
done
```

```
# Output plot from pca_plot Shiny app
# Save plot object as Rdata file
```

```R
# load plot object from rda file
load('output/pca_plot.2018-03-14.rda')
plot_data <- pca_plot$data
colour_blind_palette <- 
   c( 'blue' = rgb(0,0.45,0.7),
      'yellow' = rgb(0.95, 0.9, 0.25),
      'sky_blue' = rgb(0.35, 0.7, 0.9),
      'purple' = rgb(0.8, 0.6, 0.7)
)
names(colour_palette) <- levels(plot_data$`Delay Category`)
# print as eps file
postscript(file = 'plots/pca_plot.2018-02-23-pc3-pc4-no_names-edited.eps',
            width = 8, height = 7, horizontal = TRUE)
ggplot(data = plot_data) +
    geom_point(aes(x = PC3, y = PC4, fill = `Delay Category`, shape = Genotype), size = 4, stroke = 0.2) +
    scale_shape_manual(values = c(21,22,23)) +
    scale_fill_manual(values = colour_palette,
    guide = guide_legend(override.aes = list(shape = 21),order = 1) ) +
    theme_minimal() +
    theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14))
dev.off()

for (gt in c('wt', 'het', 'hom')) {
    postscript(file = file.path('plots', paste0('pca_plot-pc3_pc4-', gt, '-somite_number.eps') ),
                width = 8, height = 7, horizontal = TRUE)
    print(ggplot(data = plot_data[ plot_data$Genotype == gt, ]) +
        geom_point(aes(x = PC3, y = PC4, fill = Mean_somite_number, shape = `Delay Category`), size = 4, stroke = 0.2) +
        scale_shape_manual(values = c(21,22,23,24)) +
        scale_fill_viridis() +
        theme_minimal() +
        theme(axis.title = element_text(size = 14),
                axis.text = element_text(size = 12),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14))
                )
    dev.off()
}

load('output/pca-all_genes-top50000-somite_number-pc3_pc2.rda')
# plot
plot_data <- pca_plot$data
# plot homs separately from hets & wts
# set limits to same for each plot
min_x <- floor(min(plot_data$PC3))
max_x <- ceiling(max(plot_data$PC3))
min_y <- floor(min(plot_data$PC2))
max_y <- ceiling(max(plot_data$PC2))

# make new factor for hom vs het-wt
plot_data$gt_group <- as.character(plot_data$condition)
plot_data$gt_group[ plot_data$gt_group == 'wt' ] <- 'wt-het'
plot_data$gt_group[ plot_data$gt_group == 'het' ] <- 'wt-het'
plot_data$gt_group <- factor(plot_data$gt_group, level = c('hom', 'wt-het'))

postscript(file = file.path('plots', paste0('pca-all_genes-top50000-somite_number-hom_vs_het_wt.eps') ),
            width = 15, height = 7, paper = 'special')
print(ggplot(data = plot_data) +
    geom_point(aes(x = PC3, y = PC2, fill = somite_number, shape = condition), size = 4, stroke = 0.2) +
    scale_fill_viridis(name = 'Somite Number', guide = guide_colourbar(order = 1) ) +
    scale_shape_manual(values = c(21,22,23), name = 'Genotype', guide = guide_legend(order = 2)) +
    facet_wrap( ~ gt_group, nrow = 1) +
    theme_minimal() +
    theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14))
            )
dev.off()
```
