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
Rscript --vanilla /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/make_indel_info_commands.r 
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

#### Get InDel info from chromosome-only VCFs
##### R Script to generate vcftools commands to get indel INFO
```
Rscript --vanilla /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/make_indel_info_chrom_only_commands.r
```
##### Get INFO about indels from the chromsome-only VCFs
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info
chmod u+x get_chrom_only_indel_info.sh
./get_chrom_only_indel_info.sh
```
#### Divide vcfs by Insertion and Deletion
Want to be able to analyze insertions and deletions separately, if need be

##### R script to generate files with positions of insertions and deletions
```
Rscript --vanilla /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/make_sep_INS_AND_DEL_commands.r
```
##### Transfer index files to directory with original vcfs
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info
tar -cvf branch_file_in_del_positions.tar *_positions.txt
cp branch_file_in_del_positions.tar ../struc_vcfs/
cd ../struc_vcfs
tar -xvf branch_file_in_del_positions.tar

``` 
##### Generate INS and DEL-specific .vcf files for each library
```
bash /home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/gen_INS_DEL_vcfs.sh
```
#### Pairwise comparisons of INS and DEL-specific vcf files
use vcftools to find shared, overlapping, and non-shared variation in pairwise \
comparisons

##### R script to generate vcf commands for pairwise comparisons
```
Rscript --vanilla /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/make_pw_comp_commands.r
```

##### Run pairwise comparisons using vcftools
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/shared_loci
bash share_diff_comp_commands.sh
bash share_diff_INS_comp_commands.sh

```

#### Calculate closest indel for non-shared indels from pairwise comparisons
##### R script with analysis
`/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/analyze_pairwise_diffs.r`

Note: this script is pretty slow and inefficient, so submit job to run it
##### shell script to execute R script
`/home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/make_diff_indel_files.sh`
##### run shell script to do analysis
```
cd /home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts
qsub -cwd -N pop_pair_diff -l h_vmem=4G -q all.q make_diff_indel_files.sh
```
### Combine Info from each library
Goals:
1. Make table for each sample that contains distance to next closest indel 
in each of the other samples for each of the indels found in the test sample

#### Combine the indel distances from each pairwise comparison file
For each sample, go through all it's pairwise `.diff_indels` files and pull 
out the indel distances when compared to each of the other samples.
Resulting tables have `.indel_dist_tot` suffix
#### R script for generating the tables
```
Rscript --vanilla /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/gen_samp_tot_dist_tabs.r
```
### Identify Unique Indels
Goals:
1. Calculate number of Indels that are not shared by any other samples and\
are a certain distance from any other indel in any sample
2. Look at size distribution of unique indels
#### Generate Figures Showing Number and Size Distributions of Unique Indels
```
Rscript --vanilla /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/unique_indel_figs.r
```
## Indel-based Distance Matrix
Steps:
1. Generate Distance Matrices Using different bp-distances
2. PCoA with distance matrices
### Generate Distance Matrix
#### Goals
* Calculate % shared indels for each library with the other libraries
* Calculate indel-based genetic distance using both directions of the pairwise
comparison. For example, samp1Vsamp2 and samp2Vsamp1
* PCoA of distance matrices
#### R Script
Generates distance matrices and produces PCoA figures
```
Rscript --vanilla /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/make_indel_dist_mat.r
```
### PCoA Figures
Distance cutoffs: exact same position; overlap; 100bp; 1kp; 5kb; 10kb
* All indels
  * https://github.com/grabowsp/poplar_branch_indels/blob/master/figs/indel_PCoA_combo.pdf
* Insertions-only
  * https://github.com/grabowsp/poplar_branch_indels/blob/master/figs/insertions_PCoA_combo.pdf
* Deletions-only
  * https://github.com/grabowsp/poplar_branch_indels/blob/master/figs/deletions_PCoA_combo.pdf
### Interpretation
* The libraries do NOT form clone-specific or lineage-specific clusters for
any of the comparisons (different distances; insertions, deletions, or both)
* It's possible that PCoA with so few samples isn't really the best way to go
about this
## Re-call Indels with Libraries Together
### Steps
* Sort .bam files with SAM Tools
* make mpileup file with SAM Tools
* Call Indels using mpile2indel in VarScan
### SAM Tools Steps
#### R Script to Generate Commands for Sorting .bams and making mpileup
```
Rscript --vanilla /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/make_samtools_pileup_command.r
```
#### Sort .bam files
##### Shell Script
```
/home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/sort_bams.sh
```
##### Run Script
```
qsub -N sort_pop_bams -l h_vmem=8G -q all.q /home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/sort_bams.sh
```
#### Generate mpileup
##### Shell Script
```
/home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/pop_samtools_mpileup.sh
```
##### Run Script
```
qsub -N gen_pop_mpile -l h_vmem=8G -q all.q /home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/pop_samtools_mpileup.sh
```

