# Use PacBio data to check Sujan's high-quality SNPs

## Overview
* There's a question about the SNPs
* Goal is to use PacBio data to try to verify the high-quality SNPs
### Plan
* Sujan generated a list of SNPs for me
  * I still need him to explain them to me
* Try using the samtools mpileup approach

```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb
samtools mpileup -l test_snp_list.txt /home/grabowsky_scratch/poplar_branch_files/PAXL.14.5v1.0Ref.ngmlr.sorted.withRG.bam

samtools mpileup -l test_snp_list.txt /home/grabowsky_scratch/poplar_branch_files/PAXL.14.5v1.0Ref.ngmlr.sorted.withRG.bam /home/grabowsky_scratch/poplar_branch_files/PAXN.14.5v1.0Ref.ngmlr.sorted.withRG.bam > test.mpileup


```

## Plan
* Make a pileup up using all 8 libraries for highest depth (42) SNPs in both /
trees

## SNPs from Sujan
### Location of SNPs
* /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/Paul_tree13_14.positions


## Pileup results
* Directory
  * `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results`
### DONT USE: Test with tree14 depth42 SNPs and samtools
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth42_step2_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree14_depth42_step2_mpileup.sh
```
* Pileup file is essentially impossible to readthrough by eye
#### USE THIS: Try making vcf with deprecated options
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth42_step2_mpileup_v2.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree14_depth42_step2_mpileup_v2.sh
```
* This one ends up outputting genotypes like the bcf, BUT it includes allele \
read counts, which I can use

### DONT USE: test with tree14 depth42 SNPs and bcftools
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth42_step2_bcftools_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree14_depth42_step2_bcftools_mpileup.sh
```
* I couldn't get the option to include allele reads to work, so I ended up \
needing to use the deprecated samtools vcf options

## Extract genotypes
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results

cut -f 10-17 tree14_depth42_v2.vcf | awk 'BEGIN{FS=":|\t"; OFS = "\t"} \
NR<=48 {next} \
{print $2,$4,$6,$8,$10,$12,$14,$16}' > tree14_depth42_allele_counts.txt
```
* extract sample order
```
head -48 tree14_depth42_v2.vcf | tail -n 1 | \
awk '{print $10,$11,$12,$13,$14,$15,$16,$17}' > sample_order.txt
```

## Make more pileup files
### Tree14, 42 Depth, Step 3
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth42_step3_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree14_depth42_step3_mpileup.sh
```
### Tree 13, 42 Depth, Step 2
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree13_depth42_step2_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree13_depth42_step2_mpileup.sh
```
### Tree 13, 42 Depth Step 3
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree13_depth42_step3_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree13_depth42_step3_mpileup.sh
```
### Tree 14, 24 Depth, Step 2
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth24_step2_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree14_depth24_step2_mpileup.sh
```
### Tree 14, 24 Depth, Step 3
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth24_step3_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree14_depth24_step3_mpileup.sh
```
### Tree 13, 24 Depth, Step 2
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree13_depth24_step2_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree13_depth24_step2_mpileup.sh
```
### Tree 13, 24 Depth, Step 3 
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree13_depth24_step3_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree13_depth24_step3_mpileup.sh
```

## Extract Genotype Info from pileups
### Tree 14
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results

cut -f 10-17 tree14_depth42_step3.vcf | awk 'BEGIN{FS=":|\t"; OFS = "\t"} \
NR<=48 {next} \
{print $2,$4,$6,$8,$10,$12,$14,$16}' > tree14_depth42_step3_allele_counts.txt

cut -f 10-17 tree14_depth24_step2.vcf | awk 'BEGIN{FS=":|\t"; OFS = "\t"} \
NR<=48 {next} \
{print $2,$4,$6,$8,$10,$12,$14,$16}' > tree14_depth24_step2_allele_counts.txt

cut -f 10-17 tree14_depth24_step3.vcf | awk 'BEGIN{FS=":|\t"; OFS = "\t"} \
NR<=48 {next} \
{print $2,$4,$6,$8,$10,$12,$14,$16}' > tree14_depth24_step3_allele_counts.txt
```

### Tree 13 - NEED TO DO

