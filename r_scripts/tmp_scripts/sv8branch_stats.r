# Consolidate 8-branch results for manuscript

ngmlr_indel_list <- list()

for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/', 
    'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, 
    '.8branch.vcf_SVstats.txt', sep = '')
  ngmlr_indel_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
}

ngmlr_indel_df <- data.frame(matrix(
  data = unlist(lapply(ngmlr_indel_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(ngmlr_indel_df) <- ngmlr_indel_list[[1]][,1]

apply(ngmlr_indel_df[, c(1:6)], 2, mean)
#  raw BND  raw DEL  raw DUP  raw INS  raw INV  raw cnv 
# 31187.00 29273.75  2888.00 25154.50    63.00  4816.25

apply(ngmlr_indel_df[, c(1:6)], 2, sd)
#   raw BND   raw DEL   raw DUP   raw INS   raw INV   raw cnv 
# 70.814311 28.616720  5.887841 73.776690  0.000000 11.412712

mean(ngmlr_indel_df[['Raw DEL > 100bp']])
# [1] 11527.25
sd(ngmlr_indel_df[['Raw DEL > 100bp']])
# [1] 22.18671

mean(ngmlr_indel_df[['Raw INS > 100bp']])
# [1] 12674
sd(ngmlr_indel_df[['Raw INS > 100bp']])
# [1] 9.092121


### Filtered counts
# DEL
mean(ngmlr_indel_df[['Number filtered Deletions']])
# [1] 10466.25
sd(ngmlr_indel_df[['Number filtered Deletions']])
# [1] 26.42442

mean(ngmlr_indel_df[['Number DELs > 100bp']])
# [1] 4539.75
sd(ngmlr_indel_df[['Number DELs > 100bp']])
# [1] 19.4315

# INS
mean(ngmlr_indel_df[['Number filtered Insertions']])
# [1] 6702.25
sd(ngmlr_indel_df[['Number filtered Insertions']])
# [1] 39.57588

mean(ngmlr_indel_df[['Number INSs > 100bp']])
# [1] 3920.75
sd(ngmlr_indel_df[['Number INSs > 100bp']])
# [1] 8.770215

### Indels shared by all 8 branches
1 - (mean(ngmlr_indel_df[['Number variable INDELs']]) / mean(ngmlr_indel_df[['Total number of filtered INDELs']]))
# [1] 0.9899106

### Variable InDels
# DEL
mean(ngmlr_indel_df[['N gain-of-het DELs var in 1 tree']])
# [1] 4.25
sd(ngmlr_indel_df[['N gain-of-het DELs var in 1 tree']])
# [1] 0.5

mean(ngmlr_indel_df[['N gain-of-het DELs > 100bp var in 1 tree']])
# [1] 4
sd(ngmlr_indel_df[['N gain-of-het DELs > 100bp var in 1 tree']])
# [1] 0

# INS
mean(ngmlr_indel_df[['N gain-of-het INSs var in 1 tree']])
# [1] 14
sd(ngmlr_indel_df[['N gain-of-het INSs var in 1 tree']])
# [1] 3.162278

mean(ngmlr_indel_df[['N gain-of-het INSs > 100bp var in 1 tree']])
# [1] 2
sd(ngmlr_indel_df[['N gain-of-het INSs > 100bp var in 1 tree']])
# [1] 1.154701

###############################

ngmlr_dup_list <- list()

for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/', 
    'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, 
    '.8branch.vcf_DUPstats.txt', sep = '')
  ngmlr_dup_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
}

ngmlr_dup_df <- data.frame(matrix(
  data = unlist(lapply(ngmlr_dup_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(ngmlr_dup_df) <- ngmlr_dup_list[[1]][,1]

mean(ngmlr_dup_df[['Raw DUP > 100bp']])
# 1894
sd(ngmlr_dup_df[['Raw DUP > 100bp']])
# 2.708013

### Filtered DUPs
mean(ngmlr_dup_df[['Number Filtered DUP']])
# [1] 645
sd(ngmlr_dup_df[['Number Filtered DUP']])
# [1] 6.582806

mean(ngmlr_dup_df[['Number Filtered DUP > 100bp']])
# [1] 297
sd(ngmlr_dup_df[['Number Filtered DUP > 100bp']])
# [1] 1.154701

# Duplications shared by all 8 branches
1 - (mean(ngmlr_dup_df[['Number variable Filtered DUP']])/mean(
  ngmlr_dup_df[['Number Filtered DUP']]))
# [1] 0.9906977

### Variable Duplications
mean(ngmlr_dup_df[['N decent variable DUPs']])
# [1] 3
sd(ngmlr_dup_df[['N decent variable DUPs']])
# [1] 1.154701

mean(ngmlr_dup_df[['N decent DUPs > 100bp']])
# [1] 2
sd(ngmlr_dup_df[['N decent DUPs > 100bp']])
# [1] 0

#################

ngmlr_inv_list <- list()

for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/',
    'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i,
    '.8branch.vcf_INVstats.txt', sep = '')
  ngmlr_inv_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
}

ngmlr_inv_df <- data.frame(matrix(
  data = unlist(lapply(ngmlr_inv_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(ngmlr_inv_df) <- ngmlr_inv_list[[1]][,1]

### Filtered INV
mean(ngmlr_inv_df[['Number Filtered INV']])
# [1] 3
sd(ngmlr_inv_df[['Number Filtered INV']])
# [1] 0

### Variable INV
mean(ngmlr_inv_df[['Number variable Filtered INV']])
# [1] 0
sd(ngmlr_inv_df[['Number variable Filtered INV']])
# [1] 0

########

ngmlr_bnd_list <- list()

for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/',
    'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i,
    '.8branch.vcf_BNDstats.txt', sep = '')
  ngmlr_bnd_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
}

ngmlr_bnd_df <- data.frame(matrix(
  data = unlist(lapply(ngmlr_bnd_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(ngmlr_bnd_df) <- ngmlr_bnd_list[[1]][,1]

#### Filtered BND
mean(ngmlr_bnd_df[['N BND-INS with matching 5 and 3prime genotypes']])
# [1] 46.25
sd(ngmlr_bnd_df[['N BND-INS with matching 5 and 3prime genotypes']])
# [1] 0.5

### Variable BND
mean(ngmlr_bnd_df[['N Variable Good BND-INS']])
# [1] 0.5
sd(ngmlr_bnd_df[['N Variable Good BND-INS']])
# [1] 0

### Filtered cnvs
mean(ngmlr_bnd_df[['N cnv POS with adequate coverage']])
# [1] 39
sd(ngmlr_bnd_df[['N cnv POS with adequate coverage']])
# [1] 0

############

ngmlr_indel_geno_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/', 
    'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, 
    '.8branch.vcf_goodINDEL_genos.rds', sep = '')
  ngmlr_indel_geno_list[[i]] <- readRDS(tmp_name)
}

indel_name_tab <- table(unlist(
  lapply(ngmlr_indel_geno_list, function(x) rownames(x))))

which(indel_name_tab == 4)
#   Chr01_38823408_INS_67  Chr01_5943571_DEL_5282 Chr05_17223380_DEL_9729 
#                       4                       8                      22 
#  Chr09_9139125_DEL_7970 
#                      38 

# Chr01_38823408_INS_67 is variable in 13.5 (PAZF)
### in individual vcf, is position 38823410 and length of 52
# Chr01_5943571_DEL_5282 does not actually pass the criteria 
# Chr05_17223380_DEL_9729 is variable in 13.5 (PAZF)
### in individual vcf, is same position and length of 9727
# Chr09_9139125_DEL_7970 is variable in 13.2 (PBAW)
### Need to check the individual vcf - for some reason, this one was never
###   generated

#### I want to check these actual genotypes to make sure they are ok since 
## the two deletions are 2 to 1 genotype transtitions

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
meta <- read.table(meta_in, stringsAsFactors = F, sep = '\t', header = T)

geno_1_labs <- c()
for(cln in colnames(ngmlr_indel_geno_list[[1]])){
  tmp_ind <- which(meta$lib_name == cln)
  geno_1_labs <- c(geno_1_labs, meta$branch_name[tmp_ind])
}

ngmlr_dup_geno_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/',
    'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i,
    '.8branch.vcf_goodDUP_genos.rds', sep = '')
  ngmlr_dup_geno_list[[i]] <- readRDS(tmp_name)
}

dup_name_tab <- table(unlist(
  lapply(ngmlr_dup_geno_list, function(x) rownames(x))))

which(dup_name_tab == 4)
# Chr06_25962893_DUP_28528   Chr14_8492556_DUP_7858 
#                        3                        5 

# Chr06_25962893_DUP_28528 is variable in 13.2 (PBAW)
# Chr14_8492556_DUP_7858 is variable in 14.5 (PAZH)

res1_vcf <- '/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r01.8branch.vcf'
raw_vcf <- read.table(res1_vcf, header = F, stringsAsFactors = F,
  sep = '\t')

# Next: want to see where these SV's map to

#############


