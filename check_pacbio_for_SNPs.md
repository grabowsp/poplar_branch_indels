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
### Test with tree14 depth42 SNPs and samtools
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth42_step2_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree14_depth42_step2_mpileup.sh
```
* Pileup file is essentially impossible to readthrough by eye
#### Try making vcf with deprecated options
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth42_step2_mpileup_v2.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree14_depth42_step2_mpileup_v2.sh
```
### test with tree14 depth42 SNPs and bcftools
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth42_step2_bcftools_mpileup.sh`
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results
qsub tree14_depth42_step2_bcftools_mpileup.sh
```


