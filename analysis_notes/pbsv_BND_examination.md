# Examine why pbsv v2.2.0 generates so many BNDs

## Overview
### Approach
* Start looking at output for single library
  * PAXL
* Run `pbsv discover` without using the `--tandem-repeats` flag
* Run `pbsv call` on new svsig.gz file
* Compare number of SVs using ngmlr (or pbmm2?) aligner for pbsv v2.0, v2.1, \
and v2.2
### Notes
* v2.1 and origninal v2.2 results all used the tandem repeats flag

## Generate .svsig.gz file without using tandem repeats
### Script
* using ngmlr alignment
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PAXL.ngmlr.v2.2.disc.noTR.r01.sh`
### Submit script
```
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL
qsub PAXL.ngmlr.v2.2.disc.noTR.r01.sh
```

## Call SV's using .svsig.gz file made without tandem repeats flag
### Script
* using ngmlr alignment
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PAXL.ngmlr.v2.2.call.noTR.r01.sh`
### Submit script
```
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL
qsub PAXL.ngmlr.v2.2.call.noTR.r01.sh
```

## Compare v2.0, v2.1 and v2.2 outputs
### Location of Files for PAXL
* v2.0.1 output
  * `/home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL/ref.PAXL.vcf`
* v2.1.1
  * `/home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL/Ptr145v1.PAXL_v2_r1.vcf`
* v2.2.2 generated from .svsig.gz made using --tandem-repeats flag
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PtStettler14.ngmlr.pbsv_v2.2_PAXL_r01.vcf`
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PtStettler14.pbmm2.pbsv_v2.2_PAXL_r01.vcf`
* v2.2.2 generated without specifying tandem repeats
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PtStettler14.ngmlr.pbsv_v2.2_PAXL_r01.noTR.vcf`
### R script
```
v2_0_file <- '/home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL/ref.PAXL.vcf'
v2_0_vcf <- read.table(v2_0_file, header = F, stringsAsFactors = F, 
  sep = '\t')

v2_1_file <- '/home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL/Ptr145v1.PAXL_v2_r1.vcf'
v2_1_vcf <- read.table(v2_1_file, header = F, stringsAsFactors = F, 
  sep = '\t')

v2_2_n_file <- '/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PtStettler14.ngmlr.pbsv_v2.2_PAXL_r01.vcf'
v2_2_n_vcf <- read.table(v2_2_n_file, header = F, stringsAsFactors = F, 
  sep = '\t')

v2_2_p_file <- '/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PtStettler14.pbmm2.pbsv_v2.2_PAXL_r01.vcf'
v2_2_p_vcf <- read.table(v2_2_p_file, header = F, stringsAsFactors = F,
  sep = '\t')

v2_2_noTR_file <- '/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PtStettler14.ngmlr.pbsv_v2.2_PAXL_r01.noTR.vcf'
v2_2_noTR_vcf <- read.table(v2_2_noTR_file, header = F, stringsAsFactors = F,
  sep = '\t')

v2_0_raw_info <- strsplit(gsub('IMPRECISE;', '', v2_0_vcf[,8]), split = ';')
v2_0_type <- unlist(lapply(v2_0_raw_info,
  function(x) unlist(strsplit(x[[1]], split = '='))[2]))

v2_1_raw_info <- strsplit(gsub('IMPRECISE;', '', v2_1_vcf[,8]), split = ';')
v2_1_type <- unlist(lapply(v2_1_raw_info,
  function(x) unlist(strsplit(x[[1]], split = '='))[2]))

v2_2_n_raw_info <- strsplit(gsub('IMPRECISE;', '', v2_2_n_vcf[,8]), split = ';')
v2_2_n_type <- unlist(lapply(v2_2_n_raw_info,
  function(x) unlist(strsplit(x[[1]], split = '='))[2]))

v2_2_p_raw_info <- strsplit(gsub('IMPRECISE;', '', v2_2_p_vcf[,8]), split = ';')
v2_2_p_type <- unlist(lapply(v2_2_p_raw_info,
  function(x) unlist(strsplit(x[[1]], split = '='))[2]))

v2_2_noTR_raw_info <- strsplit(gsub('IMPRECISE;', '', v2_2_noTR_vcf[,8]), 
  split = ';')
v2_2_noTR_type <- unlist(lapply(v2_2_noTR_raw_info,
  function(x) unlist(strsplit(x[[1]], split = '='))[2]))

table(v2_0_type)
#   BND   DEL   INS   INV 
#  1066 26439 24022     7

table(v2_1_type)
#   BND   DEL   INS   INV 
# 10364 25273 25358    21

table(v2_2_n_type)
#   BND   cnv   DEL   DUP   INS   INV 
# 11418   623 26177  1686 27251    30

table(v2_2_p_type)
#   BND   cnv   DEL   DUP   INS   INV 
# 10598   422 25821  1425 30320    74

table(v2_2_noTR_type)
#   BND   cnv   DEL   DUP   INS   INV 
# 11418   623 26743  2044 23631    30
```

