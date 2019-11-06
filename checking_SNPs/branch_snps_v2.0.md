# Analysis of the second set of SNPs for the branches using new criteria

## Files from Sujan
### Locations in his directory
```
/home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/treewise_SNP_3steps/probability/tree04_pvalues_patterns.txt.bz2
/home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/treewise_SNP_3steps/probability/tree09_pvalues_patterns.txt.bz2
/home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/treewise_SNP_3steps/probability/tree13_pvalues_patterns.txt.bz2
/home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/treewise_SNP_3steps/probability/tree14_pvalues_patterns.txt.bz2
/home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/treewise_SNP_3steps/probability/tree15_pvalues_patterns.txt.bz2
```
### Copy to my directory
* /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519
```
cp /home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/treewise_SNP_3steps/probability/*.bz2 /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/
```
### Locations in my directory
```
/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree04_pvalues_patterns.txt.bz2
/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree09_pvalues_patterns.txt.bz2
/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree13_pvalues_patterns.txt.bz2
/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_pvalues_patterns.txt.bz2
/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree15_pvalues_patterns.txt.bz2
```
## Illumina-based VCF from Sujan
### Tree14
* Sujan's location
  * `/home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/treewise_SNP_3steps/probability/tree14.good_positions.v1.vcf.gz`
* My location
```
cp /home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/treewise_SNP_3steps/probability/tree14.good_positions.v1.vcf.gz /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/
```
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14.good_positions.v1.vcf`
### Tree13
* Sujan's location
  * `/home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/treewise_SNP_3steps/probability/tree13.good_positions.v1.vcf.gz`
* My location
```
cp /home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/treewise_SNP_3steps/probability/tree13.good_positions.v1.vcf.gz /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/
```
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree13

.good_positions.v1.vcf`

## Explore Tree14 data
### Make smaller file for Tree14


```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519
bzcat tree14_pvalues_patterns.txt.bz2 | head -10000 > tree14_10k.txt
```
### Look at data in R
* `/home/grabowsky/tools/workflows/poplar_branch_indels/checking_SNPs/tree14_snps_v2.0_analysis.r`

### Positions that pass filtering
* using 1-(1/10,000) as the cutoff for a good genotype
#### Positions for samtools mpileup
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_good_positions_v1.txt`
* Chr01 positions for testing
  * `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_good_pos_Chr01_v1.txt`
#### Make pseudo-vcf for Jerry's approach
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_good_pos_v1.vcf`
* Chr01 vcf for testing
  * `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_good_pos_Chr01_v1.vcf`
##### Steps for making VCF
* "Meat" of VCF
  * `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_good_pos_vcf_meat.txt`
* Header made with vim
  * `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_good_pos_vcf_head.txt`
* Make full file
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519
cat tree14_good_pos_vcf_head.txt tree14_good_pos_vcf_meat.txt > \
tree14_good_pos_v1.vcf
```
#### Extract SNP names and alleles from p-value file
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519
cut -f 1-4,9 tree14_pvalues_patterns.txt > tree14_pval_file_SNPinfo.txt

```
* File with position and allele info for Tree14 SNPs from Sujan
  * `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_pval_file_SNPinfo.txt`

### Make Pileup/VCF using positions
#### Test with Chr01 positions
* Copy command
```
cp /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth42_step2_mpileup_v2.sh /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_v2SNPs_mpileup_Chr01.sh
```
 * adjusted with vim
* Submit command
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519
qsub tree14_v2SNPs_mpileup_Chr01.sh
```
##### Extract genotypes
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519

cut -f 10-17 tree14_v2SNPs_PB.vcf | awk 'BEGIN{FS=":|\t"; OFS = "\t"} \
NR<=48 {next} \
{print $2,$4,$6,$8,$10,$12,$14,$16}' > tree14_v2SNPs_PB_Chr01_allele_counts.txt
```
* extract sample order
```
head -48 tree14_v2SNPs_PB.vcf | tail -n 1 | \
awk '{print $10,$11,$12,$13,$14,$15,$16,$17}' > tree14_PB_sample_order.txt
```
#### Run with full data set
* Copy command
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519

cp tree14_v2SNPs_mpileup_Chr01.sh tree14_v2SNPs_mpileup_full.sh
```
* adjusted with vim
* Submit command
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519
qsub tree14_v2SNPs_mpileup_full.sh
```
* Output file:
  * `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_v2SNPs_PB_full.vcf`
##### Extract genotypes
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519

cut -f 10-17 tree14_v2SNPs_PB_full.vcf | awk 'BEGIN{FS=":|\t"; OFS = "\t"} \
NR<=48 {next} \
{print $2,$4,$6,$8,$10,$12,$14,$16}' > tree14_v2SNPs_PB_full_allele_counts.txt
```

### Use Jerry's approach
* Directory with files
  `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method`
#### Test with Chr01 positions
##### Read-to-variant
* Copy command
```
cp /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/PAXL_t14_s2_Chr01_read_to_variant.sh /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method/PAXL_t14_v2_Chr01_read_to_variant.sh
```
* Make FOFN
```
echo /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_good_pos_Chr01_v1.vcf > /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method/tree14_test.fofn
```
* Submit job
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method
qsub PAXL_t14_v2_Chr01_read_to_variant.sh
```
* I had to adjust the command to remove using the temporary folder - there \
was an issue with saving and/or transfering the output from the temporary \
folder, so I'm running from the main directory
* Note: the output is saved to the directory that contains the BAM files
 
##### Call Ref and Alt alleles for each position 
* Copy command
```
cp /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/PAXL_t14_s2_Chr01_callRefAlt.sh /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method/PAXL_t14_v2_Chr01_callRefAlt.sh
```
* adjust with vim
* Submit job
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method
qsub PAXL_t14_v2_Chr01_callRefAlt.sh
```
##### Copy output to correct directory
```
cd /home/grabowsky_scratch/poplar_branch_files

mv PAXL.14.5v1.0Ref.ngmlr.sorted.withRG.outlier.dat \
PAXL_Tree14_Chr01.outlier.dat
mv PAXL.14.5v1.0Ref.ngmlr.sorted.withRG.phase_info.dat \
PAXL_Tree14_Chr01.phase_info.dat
mv PAXL.14.5v1.0Ref.ngmlr.sorted.withRG.read_to_variant.dat \
PAXL_Tree14_Chr01.read_to_variant.dat

mv /home/grabowsky_scratch/poplar_branch_files/*.dat \
/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method/
```

#### Full Tree14 SNP set for PAXL
##### Make FOFN
```
echo /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_good_pos_v1.vcf > /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method/tree14_full.fofn
```
##### Read-to-variant
* Copy and adjust command
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method
cp  PAXL_t14_v2_Chr01_read_to_variant.sh PAXL_t14_v2_full_read_to_variant.sh 
```
* Submit job
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method
qsub PAXL_t14_v2_full_read_to_variant.sh
```
##### Call Ref/Alt Alleles
* Make script and adjust with vim
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method
cp PAXL_t14_v2_Chr01_callRefAlt.sh PAXL_t14_v2_full_callRefAlt.sh
```
* submit job
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method
qsub PAXL_t14_v2_full_callRefAlt.sh
```
##### Copy output to correct directory
```
cd /home/grabowsky_scratch/poplar_branch_files

mv PAXL.14.5v1.0Ref.ngmlr.sorted.withRG.outlier.dat \
PAXL_Tree14_full.outlier.dat
mv PAXL.14.5v1.0Ref.ngmlr.sorted.withRG.phase_info.dat \
PAXL_Tree14_full.phase_info.dat
mv PAXL.14.5v1.0Ref.ngmlr.sorted.withRG.read_to_variant.dat \
PAXL_Tree14_full.read_to_variant.dat

mv /home/grabowsky_scratch/poplar_branch_files/*.dat \
/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method/
```

### Jerry's approach for remaining libraries
#### Generate Read-to-Variant scripts
```
bash
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method

for LIB in PBAU PBAW PBAT PAZF PAXN PAYK PAZH;
do 
sed 's/PAXL/'"$LIB"'/g' PAXL_t14_v2_full_read_to_variant.sh > \
$LIB'_t14_v2_full_read_to_variant.sh';
done
``` 
#### Submit Read-to-variant scripts
```
bash
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_jwj_method

for LIB in PBAU PBAW PBAT PAZF PAXN PAYK PAZH;
do qsub $LIB'_t14_v2_full_read_to_variant.sh';
done
```

## Generate SNP lists for Tree13 and Tree14
### Overview
* Use Illumina data to find SNPs with 1 high-confidence HET and 1 \
high-confidence HOM genotype in the branches
* Compare Illumina SNPs to genotypes from PacBio data
* Select SNPs that are 3-HOM:1-HET and have the same genotype calls in both \
datasets
### Tree14
#### P-value file based on Illumina read counts
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_pvalues_patterns.txt`
#### R script for selecting Illumina-PASS SNPs
* `/home/grabowsky/tools/workflows/poplar_branch_indels/checking_SNPs/tree14_snps_v2.0_analysis.r`
* Outputs from this script are `tree14_good_positions_v1.txt` and \
`tree14_good_positions_v1.dosages` described below
  * also generates a pseudo-VCF file to be used with Jerry's PacBio method
#### List of Illumina-PASS SNPS
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_good_positions_v1.txt`
* These SNPs have:
  * No missing data for the branches
  * Min 21 seq depth in all branches
  * 1+ branch with high-confidence HOM genotype call
  * 1+ branch with high-confidence HET genotype call
  * High-confidence = Conditional probability of genotype is 10,000X times \
higher than the other 2 genotype probabilities
#### Dosage genotypes for Illumina-PASS SNPs
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_good_positions_v1.dosages`
* These postions capture the true genotype a bit better than the VarScan \
method and match up much better with the PacBio genotypes
  * Maybe because VarScan considers MAF, which shouldn't matter in this case? \
I'm not quite sure what the reason is
#### Illumina VCF for Illumina-PASS SNPs
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14.good_positions.v1.vcf`
#### Generate pileup file for Illuina-PASS SNP positions from PacBio Data
PUT INFO HERE
#### R script for selecting "Good" Somatic SNPs
* `/home/grabowsky/tools/workflows/poplar_branch_indels/checking_SNPs/tree14_PB_pileup_analysis.r`
#### Positions of "Good" Somatic SNPS
*  `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_filtered_somatic_SNPs.positions`
#### Generate Illumina VCF for "Good" Somatic SNPs
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/
vcftools --vcf tree14.good_positions.v1.vcf -c \
--positions tree14_filtered_somatic_SNPs.positions --recode --recode-INFO-all \
> tree14.filtered_somatic_SNPs.v1.vcf
```
#### Generate Illumina per-sample and per-site depth file for Illumin-PASS SNPs
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519
vcftools --vcf tree14.good_positions.v1.vcf -c --geno-depth > \
tree14.good_positions.gdepth
```
#### Positions of "Expanded-Good" Somatic SNPs
* this includes SNPs that fit the expected branch order
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_filtered_somatic_SNPs_expanded.positions`
#### Generate Illumina VCF for "Expanded-Good" Somatic SNPs
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/
vcftools --vcf tree14.good_positions.v1.vcf -c \
--positions tree14_filtered_somatic_SNPs_expanded.positions --recode \--recode-INFO-all \
> tree14.filtered_somatic_SNPs.v2.vcf
```



### Tree13
#### P-value file based on Illumina read counts
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree13_pvalues_patterns.txt`
#### Extract SNP names and alleles from p-value file
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519
cut -f 1-4,9 tree13_pvalues_patterns.txt > tree13_pval_file_SNPinfo.txt

```
* File with position and allele info for Tree13 SNPs from Sujan
  * `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree13_pval_file_SNPinfo.txt`

#### R script for selecting Illumina-PASS SNPs
* `/home/grabowsky/tools/workflows/poplar_branch_indels/checking_SNPs/tree13_snps_v2.0_analysis.r`
* Outputs from this script are `tree13_good_positions_v1.txt` and \
`tree13_good_positions_v1.dosages` described below
  * also generates a pseudo-VCF file to be used with Jerry's PacBio method
#### List of Illumina-PASS SNPS
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree13_good_positions_v1.txt`
#### Dosage genotypes for Illumina-PASS SNPs
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree13_good_positions_v1.dosages`
#### Generate pileup file for Illumina-PASS positions in PacBio data
##### Generate Command
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519

cp tree14_v2SNPs_mpileup_full.sh tree13_v2SNPs_mpileup_full.sh
```
* adjusted with vim
* Submit command
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519
qsub tree13_v2SNPs_mpileup_full.sh
```
* Output file:
  * `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree13_v2SNPs_PB_full.vcf`
##### Extract genotypes
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519

cut -f 10-17 tree13_v2SNPs_PB_full.vcf | awk 'BEGIN{FS=":|\t"; OFS = "\t"} \
NR<=48 {next} \
{print $2,$4,$6,$8,$10,$12,$14,$16}' > tree13_v2SNPs_PB_full_allele_counts.txt
```
#### R script for selecting "Good" Somatic SNPs
* `/home/grabowsky/tools/workflows/poplar_branch_indels/checking_SNPs/tree13_PB_pileup_analysis.r`

#### Positions of "Good" Somatic SNPS
*  `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree13_filtered_somatic_SNPs.positions`
#### Generate Illumina VCF for "Good" Somatic SNPs
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/
vcftools --vcf tree13.good_positions.v1.vcf -c \
--positions tree13_filtered_somatic_SNPs.positions --recode --recode-INFO-all \
> tree13.filtered_somatic_SNPs.v1.vcf
```
#### Generate Illumina per-sample and per-site depth file for Illumin-PASS SNPs
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519
vcftools --vcf tree13.good_positions.v1.vcf -c --geno-depth > \
tree13.good_positions.gdepth
```
#### Positions of "Expanded-Good" Somatic SNPs
* this includes SNPs that fit the expected branch order
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree13_filtered_somatic_SNPs_expanded.positions`
#### Generate Illumina VCF for "Expanded-Good" Somatic SNPs
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/
vcftools --vcf tree13.good_positions.v1.vcf -c \
--positions tree13_filtered_somatic_SNPs_expanded.positions --recode \--recode-INFO-all \
> tree13.filtered_somatic_SNPs.v2.vcf
```

## Generating Down-sampled VCFs
* downsampling the Illumina Depth so all libraries have same depth
### Directory for files
* `/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/downsample`
  * SNP position files generated in the subsampling analysis R scripts
### Generate VCFs
#### Tree13
```
bash
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/downsample
for T13POS in tree13.somatic*;
do vcftools --vcf \
/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree13.good_positions.v1.vcf \
-c --positions $T13POS --recode --recode-INFO-all \
> $T13POS.vcf;
done
```
#### Tree14
```
bash
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/downsample
for T14POS in tree14.somatic*;
do vcftools --vcf \
/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14.good_positions.v1.vcf \
-c --positions $T14POS --recode --recode-INFO-all \
> $T14POS.vcf;
done
```
### zip files together
```
cd /home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/downsample
tar cfzv tree13.downsampled.vcfs.tar.gz tree13.somatic*vcf
tar cfzv tree14.downsampled.vcfs.tar.gz tree14.somatic*vcf


```
