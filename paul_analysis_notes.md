# Steps used in analysis of indel variation in poplar branches
## General InDel Info
Steps:
1. Exploratory analysis

### Exploratory Analysis 
Goals: 
1. Extract info from .vcf files about the type (insertion or deletion) and \
size of the indels identified in each library.
2. Look at/for patterns of overall numbers and size distributions of indels \
in each of the libraries

#### R Script to generate vcftools commands get indel INFO
```
Rscript --vanilla /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/make_indel_info_commands.r` 
```

#### Get INFO about indels using vcftools
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info
chmod u+x get_indel_info.sh
./get_indel_info.sh
```

#### Make figures showing number and size distribution of Indels in each library
```
Rscript --vanilla /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/indel_gen_info_figs.r
```
Figures are found here:
`/home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info/figs/pop_branch_InDel_tally.png`
`/home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info/figs/pop_branch_size_hists_combo.pdf`

## Find Unique and Shared InDels
Steps:
1. Pairwise comparisons of libraries
2. Combine info for each library
3. Identify unique InDels

### Pairwise comparisons of libraries
Goals:
1. Find shared, overlapping, and non-shared indels in pairwise comparisons of \
libraries
2. Calculate closest indels to non-shared indels in each of the pairwise \
comparisons

#### Make VCFs that only contain chromosomes
Different  files would contain indels on different non-anchored contigs, 
so need to make vcf files that only contained positions anchored on chromosomes

Shell script to use vcftools to generate the chromosome-only .vcf files
```
cd /home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts
bash make_chrom_only_vcfs.sh
```


 
