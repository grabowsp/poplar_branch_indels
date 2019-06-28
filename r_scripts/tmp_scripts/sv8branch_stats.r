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

mean(ngmlr_indel_df[['Number DELs > 1kbp']])
# [1] 841.25
sd(ngmlr_indel_df[['Number DELs > 1kbp']])
# [1] 8.732125

mean(ngmlr_indel_df[['Number DELs > 5kbp']])
# [1] 208
sd(ngmlr_indel_df[['Number DELs > 5kbp']])
# [1] 5.09902

mean(ngmlr_indel_df[['Number DELs > 10kbp']])
# [1] 83.25
sd(ngmlr_indel_df[['Number DELs > 10kbp']])
# [1] 1.707825

mean(ngmlr_indel_df[['Number DELs > 25kbp']])
# [1] 25.75
sd(ngmlr_indel_df[['Number DELs > 25kbp']])
# [1] 1.258306

mean(ngmlr_indel_df[['Number DELs > 50kbp']])
# [1] 11
sd(ngmlr_indel_df[['Number DELs > 50kbp']])
# [1] 0

# INS
mean(ngmlr_indel_df[['Number filtered Insertions']])
# [1] 6702.25
sd(ngmlr_indel_df[['Number filtered Insertions']])
# [1] 39.57588

mean(ngmlr_indel_df[['Number INSs > 100bp']])
# [1] 3920.75
sd(ngmlr_indel_df[['Number INSs > 100bp']])
# [1] 8.770215

mean(ngmlr_indel_df[['Number INSs > 1kbp']])
# 617.75
sd(ngmlr_indel_df[['Number INSs > 1kbp']])
# 3.86221

mean(ngmlr_indel_df[['Number INSs > 5kbp']])
# 28.25
sd(ngmlr_indel_df[['Number INSs > 5kbp']])
# 0.5

mean(ngmlr_indel_df[['Number INSs > 10kbp']])
# 0
sd(ngmlr_indel_df[['Number INSs > 10kbp']])
# 0

mean(ngmlr_indel_df[['Number INSs > 25kbp']])
# 0
sd(ngmlr_indel_df[['Number INSs > 25kbp']])
# 0

mean(ngmlr_indel_df[['Number INSs > 50kbp']])
# 0
sd(ngmlr_indel_df[['Number INSs > 50kbp']])
# 0

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

mean(ngmlr_dup_df[['Number Filtered DUP > 1kbp']])
# [1] 88.5
sd(ngmlr_dup_df[['Number Filtered DUP > 1kbp']])
# [1] 0.5773503

mean(ngmlr_dup_df[['Number Filtered DUP > 5kbp']])
# [1] 59
sd(ngmlr_dup_df[['Number Filtered DUP > 5kbp']])
# [1] 0.8164966

mean(ngmlr_dup_df[['Number Filtered DUP > 10kbp']])
# [1] 39
sd(ngmlr_dup_df[['Number Filtered DUP > 10kbp']])
# [1] 0

mean(ngmlr_dup_df[['Number Filtered DUP > 25kbp']])
# [1] 14
sd(ngmlr_dup_df[['Number Filtered DUP > 25kbp']])
# 0

mean(ngmlr_dup_df[['Number Filtered DUP > 50kbp']])
# [1] 5
sd(ngmlr_dup_df[['Number Filtered DUP > 50kbp']])
# 0

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

###########3
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

mean(ngmlr_inv_df[['Number Filtered INV']])
# 3
sd(ngmlr_inv_df[['Number Filtered INV']])
# 0

mean(ngmlr_inv_df[['N Filtered INV > 100bp']])
# 3
sd(ngmlr_inv_df[['N Filtered INV > 100bp']])
# 0

mean(ngmlr_inv_df[['N Filtered INV > 1kbp']])
# 1
sd(ngmlr_inv_df[['N Filtered INV > 1kbp']])
# 0

mean(ngmlr_inv_df[['N Filtered INV > 5kbp']])
# 0
sd(ngmlr_inv_df[['N Filtered INV > 5kbp']])
#0

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
# Look for loss-of-heterozygosity
all_indel_geno_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/',
    'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i,
    '.8branch.vcf_allVarINDEL_genos.rds', sep = '')
  all_indel_geno_list[[i]] <- readRDS(tmp_name)
}

all_indel_name_tab <- table(unlist(
  lapply(all_indel_geno_list, function(x) rownames(x))))

indel_overlap <- names(all_indel_name_tab)[which(all_indel_name_tab == 4)]

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
meta <- read.table(meta_in, stringsAsFactors = F, sep = '\t', header = T)

test_genos <- all_indel_geno_list[[1]][indel_overlap, ]
for(i in seq(ncol(test_genos))){
  test_lib <- colnames(test_genos)[i]
  meta_ind <- which(meta$lib_name == test_lib)
  colnames(test_genos)[i] <- paste('b_',meta$branch_name[meta_ind], sep = '')
}

branch_num_ord <- c(13.5,13.3,13.2,13.1,14.5,14.4,14.3,14.2)
branch_name_ord <- paste('b_', branch_num_ord, sep = '')

test_genos_2 <- test_genos[, branch_name_ord]
# nothing in tree 13...
test_genos_2[, c(5:8,1:4)]
# nothing in tree 14
## No loss-of-het shared by 2+ branches in a way that we would expect

all_dup_geno_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/',
    'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i,
    '.8branch.vcf_allVarDUP_genos.rds', sep = '')
  all_dup_geno_list[[i]] <- readRDS(tmp_name)
}

all_dup_name_tab <- table(unlist(
  lapply(all_dup_geno_list, function(x) rownames(x))))

dup_overlap <- names(all_dup_name_tab)[which(all_dup_name_tab == 4)]

dup_genos <- all_dup_geno_list[[1]][dup_overlap, ]
# only 2, and none show a characteristic loss-of-het pattern


