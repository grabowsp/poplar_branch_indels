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


## Check SNPs with Jerry's PB version of mpileup
### Make softlinks to bams
```
bash
cd /home/grabowsky_scratch/poplar_branch_files/
for i in *RG.bam;
do ln -s /home/grabowsky_scratch/poplar_branch_files/$i /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/$i;
done
```
### Link to reference
```
ln -s /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/Populus_trichocarpa_var_14.5.mainGenome.fasta
```

### Jerry's Test run
#### Copy Sujan's VCF
```
cp /home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/vcfs_to_send/treewise_snp/tree13_step2.vcf.gz /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/tree13_step2.vcf.gz
```
### Run with old Tree14 SNPs as a test
#### Copy Sujan's VCF
```
cp /home/smamidi_scratch/Ptricocarpa_V14.5_somatic_mutations/SNP_treewise/vcfs_to_send/treewise_snp/tree14_step2.vcf.gz /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/tree14_step2.vcf.gz
```
#### Make vcf.fofn
* 'file of file names'
  * file with the full path of the vcf used for this
```
echo \
/home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/tree14_step2.vcf \
> \
/home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/tree14_step2_vcf.fofn
```
#### Make submit script
```
cp /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree13_depth24_step3_mpileup.sh /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/PAXL_tree14_step2_jerryPB_look.sh
```
* adjust shell script to have commands that Jerry sent me
* NOTE: NEED TO ADD PATH to .py FILES
#### Submit job
```
cd /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps
qsub PAXL_tree14_step2_jerryPB_look.sh
```
* waiting for this to finish running

## Run on subsets of VCF
### Overview
* It's hard to run and process the full results, particularly for large sets \
of SNPs
* Should try breaking up the VCF into sub-files and running the mpileup \
and/or Jerry's approach separately on each of the sub-files
### Steps
* Find out best way to divide VCF
 * Will use VCF tools because retains Chromosome names in the process
 * I tried plink and it's fast but it renames the chromosomes in the \
vcf to numerical values rather than the original chromosome positions that \
were used for mapping
* Run mpileup and Jerry's approach on on sub-VCF from Sujan's old results
* Figure out how to process the outputs
### Divide VCF
#### Plink (old)
```
cd /home/grabowsky_scratch/poplar_branch_files/sub_vcfs
plink --vcf \
/home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/tree14_step2.vcf \
--keep-allele-order --make-bed --allow-extra-chr --out tree14_step2

plink --bfile tree14_step2 --allow-extra-chr --chr Chr01 --recode vcf --out \
tree14_step2_Chr01 

```
#### VCFtools (USE)
```
cd /home/grabowsky_scratch/poplar_branch_files/sub_vcfs

vcftools --vcf \
/home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/tree14_step2.vcf \
--chr Chr01 --out tree14_step2_try2 --recode --recode-INFO-all

vcftools --vcf \
/home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/tree14_step2.vcf \
--chr Chr01 --out tree14_step2_Chr01_tmp --kept-sites

sed '1d' tree14_step2_Chr01_tmp.kept.sites > tree14_step2_Chr01.positions
```

### mpileup for Chr01 SNPs
#### make submit script
```
cd /home/grabowsky_scratch/poplar_branch_files/sub_vcfs
cp /home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth24_step3_mpileup.sh tree14_step2_Chr01_mpileup.sh
```
* adjust with vim
#### submit
```
cd /home/grabowsky_scratch/poplar_branch_files/sub_vcfs
qsub tree14_step2_Chr01_mpileup.sh
```

### Run Jerry's program
#### Copy sub-VCF to directory with necessary files for testing
```
cp /home/grabowsky_scratch/poplar_branch_files/sub_vcfs/tree14_step2_Chr01.vcf \
/home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/
```
#### Make submission script for read-to-variant python command
* I had problems running the two commands together, so I'll try running \
them separately
```
cd /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/
cp PAXL_tree14_step2_jerryPB_look.sh PAXL_t14_s2_Chr01_read_to_variant.sh
```
* adjust with vim

#### make FOFN
```
echo /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/tree14_step2_Chr01.vcf > /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/tree14_step2_Chr01_vcf.fofn
```
#### Run read-to-variant
```
cd /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/
qsub PAXL_t14_s2_Chr01_read_to_variant.sh
```
#### make submission script for call-ref-alt command
```
cd /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/
cp PAXL_tree14_step2_jerryPB_look.sh PAXL_t14_s2_Chr01_callRefAlt.sh
```
* adjust with vim
#### Run call-ref-alt (NEED TO DO)
```
cd /home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/
qsub PAXL_t14_s2_Chr01_callRefAlt.sh
```


## Try to parse Jerry's output in R
```
test_file <- '/home/grabowsky_scratch/poplar_branch_files/pb_jerry_snps/PAXL.14.5v1.0Ref.ngmlr.sorted.withRG.phase_info.dat'

test <- scan(test_file, what = 'character', sep = '$')

dash_inds <- grep('---------', test)
read_name_inds <- dash_inds + 1
read_start_inds <- read_name_inds + 1
read_end_inds <- c(dash_inds[c(2:length(dash_inds))]-1, length(test))

info_list <- list()
#for(i in seq(length(read_name_inds))){
for(i in seq(100)){
  info_list[[i]] <- list()
  info_list[[i]][[1]] <-test[read_name_inds[i]] 
  info_list[[i]][[2]]<- test[read_start_inds[i]:read_end_inds[i]]
}

# need to adjust the indexing below...

read_sum_list <- lapply(info_list, function(x) data.frame(matrix(unlist(
  strsplit(info_list[[1]], split = '\t|;')), byrow = T, ncol = 5), 
  stringsAsFactors = F))

read_sum_list <- lapply(read_sum_list, function(x) x)

```
