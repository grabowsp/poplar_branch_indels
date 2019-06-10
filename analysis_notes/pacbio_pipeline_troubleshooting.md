# Code for troubleshooting the sample order-dependent likey artifact in \
the PacBio SV calling pipeline in my poplar data

## Find differnces between the single-sample VCF and the combined VCF
### Location of Files
#### Location of Original Files
* Combined VCF
  * `/home/f1p1/tmp/PBSV/Poplar14.5/LIBS`
* Individual VCFs
  * ex: PAXL
    * `/home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL`
#### Location of files in Work Directory
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller`
#### Location of chromosome-only files
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only`

### Look at differences between PBAW (13.2) and combo VCF
* PBAW (13.2) was last library in the command to make the combined VCF so \
will compare that to the combined VCF
#### VCFtools to find non-overlapping sites
```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller

vcftools --vcf /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.PBAW.vcf_chromOnly.vcf.recode.vcf --diff /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.ALLData.vcf_chromOnly.vcf.recode.vcf --diff-site --out /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PBAW_13_2_v_combo
```
#### Info from VCFtools
```
Found 32921 sites common to both files.
Found 9526 sites only in main file.
Found 16849 sites only in second file.
Found 9150 non-matching overlapping sites.

```
#### Look at output file in R
```
diff_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PBAW_13_2_v_combo.diff.sites_in_files'

diff_info <- read.table(diff_file, header = T, stringsAsFactors = F, sep = '\t')

diff_info$sv_len_1 <- nchar(diff_info$REF1)
diff_info$sv_len_2 <- nchar(diff_info$REF2)
diff_info$sv_len_3 <- nchar(diff_info$ALT1)
diff_info$sv_len_4 <- nchar(diff_info$ALT2)
diff_info$sv_len_max <- apply(
  diff_info[ , c('sv_len_1', 'sv_len_2', 'sv_len_3', 'sv_len_4')], 1, max)

insamp1_inds <- which(diff_info$IN_FILE == '1')

insamp1_100bp_inds <- intersect(insamp1_inds, 
  which(diff_info$sv_len_max >= 100))

# insapm1_100bp_inds[1] - same SV but off by 10bp - should check to make sure it is genotypes in PBAW in the combo file because is should be
## yes, it's genotyped

# insamp1_10obp_inds[2] - is a big INDEL, but it's right next to what seems
#  to be a translocation, at least it is in the PBAW file

# insamp1_100bp_inds[3] - similar to [1] same SV just off by a few bp - make
#   sure it's genotyped in the combo file
## yes, its genotyped

# others that are like [1] and [3]:
## 4, 5, 6, 7, 9, 10, 11, 13,14,15,16,17, 20, 21, 22, 23, 24,25,26,27,28,29,
### 30, 31,32,33,34,36, 39, 40,41,42,43,45,46,47, 48,49, 50

# indices with sequences that are hard to trust: 38

# insamp1_100bp_inds[8] : a 4490bp SV at position 1064395 is called as the same
##  SV as a 3984 SV at position 1064502 in the combo VCF. HOWEVER, PBAW ALSO
## has a 3982bp SV at position 1064506 which is the one that should be linked
## with the combo.vcf SV. So in the end, the 4490bp SV in PBAW is missing
## from the combo.vcf

# insamp1_100bp_inds[12]: a 247 SV at position 2177756 is just straight up 
##  missing from the combo.vcf; closest non-overlap SV positon in combo.vcf
## is >3kb away, and the sequence isn't overly repetitive

# insamp1_100bp_inds[18]: very similar pattern to [8]: a 978bp SV at position
## 3351906 in PBAW is called as overlap with a 102bp SV at position 3352587 in
## the combo.vcf. However, PBAW has a 102bp SV at position 3352587, so this
## SV is missing

## note: these mis-called overlaps are an issue with VCFtools and not the 
### PacBio pipeline, but they DO show that close-by SVs seem to be missed
### in PBAW

# insamp1_100bp_inds[19]: a 255 bp at position 5285244 in PBAW is just
## straight up missing in the combo.vcf. Closest non-overlap SV position in
## combo.vcf is 192 bp away, which is shorter than the length of the SV

# insamp1_100bp_inds[35]: a 107bp SV is just straight up missing in the 
#  combo.vcf. Closest non-overlapping position in combo.vcf is 300bp away

# insamp1_100bp_inds[37]: a 149bp SV at position 7377273 in PBAW is just
## missing from the combo.vcf, the sequence looks fine and there are plenty
## of reads in the PBAW solo vcf to support the SV; closest non-overlapping
## SV in combo.vcf is 7kb away

# insamp1_100bp_inds[44]: a 196bp SV at position 8651167 missing from combo.vcf
## closest non-overlapping positon in combo.vcf is 1.6kb away

```
### Look at rest of libraries using VCFtools

```
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller

vcftools --vcf /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.PAXL.vcf_chromOnly.vcf.recode.vcf --diff /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.ALLData.vcf_chromOnly.vcf.recode.vcf --diff-site --out /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PAXL_14_3_v_combo

vcftools --vcf /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.PAXN.vcf_chromOnly.vcf.recode.vcf --diff /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.ALLData.vcf_chromOnly.vcf.recode.vcf --diff-site --out /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PAXN_14_2_v_combo

vcftools --vcf /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.PAYZ.vcf_chromOnly.vcf.recode.vcf --diff /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.ALLData.vcf_chromOnly.vcf.recode.vcf --diff-site --out /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PAYZ_14_1_v_combo

vcftools --vcf /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.PAYK.vcf_chromOnly.vcf.recode.vcf --diff /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.ALLData.vcf_chromOnly.vcf.recode.vcf --diff-site --out /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PAYK_14_4_v_combo

vcftools --vcf /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.PAZH.vcf_chromOnly.vcf.recode.vcf --diff /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.ALLData.vcf_chromOnly.vcf.recode.vcf --diff-site --out /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PAZH_14_5_v_combo

#

vcftools --vcf /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.PBAU.vcf_chromOnly.vcf.recode.vcf --diff /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.ALLData.vcf_chromOnly.vcf.recode.vcf --diff-site --out /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PBAU_13_1_v_combo

vcftools --vcf /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.PBAT.vcf_chromOnly.vcf.recode.vcf --diff /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.ALLData.vcf_chromOnly.vcf.recode.vcf --diff-site --out /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PBAT_13_3_v_combo

vcftools --vcf /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.PAZG.vcf_chromOnly.vcf.recode.vcf --diff /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.ALLData.vcf_chromOnly.vcf.recode.vcf --diff-site --out /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PAZG_13_4_v_combo

vcftools --vcf /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.PAZF.vcf_chromOnly.vcf.recode.vcf --diff /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only/ref.ALLData.vcf_chromOnly.vcf.recode.vcf --diff-site --out /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/pb_pipe_troubleshoot/PAZF_13_5_v_combo

```

## Rerun the pacbio command with different sample orders
### Overview
* run the pacbio `call` program 3 more times but with different sample orders
* Want to try to re-create the bug using the different sample orders
* Analyses:
  * Filter SVs and call genotypes, then tally the different genotype classes
  * Look at overlap between the different VCFs if still see the bug
### shell scripts
* in `/home/f1p1/tmp/PBSV/Poplar14.5/LIBS`
  * `full.call.try2.sh`
  * `full.call.try3.sh`
  * `full.call.try4.sh`
### Analysis
#### Move files to working directory
* `/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller`
  * `ref.ALLData_try2.vcf`
  * `ref.ALLData.try3.vcf`
  * `ref.ALLData.try4.vcf`


