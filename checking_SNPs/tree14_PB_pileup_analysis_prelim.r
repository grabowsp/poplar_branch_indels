# Script for analyzing the mpileup results for Tree14

## LOAD FUNCTIONS AND PACKAGES
func_file <- paste('/home/grabowsky/tools/workflows/poplar_branch_indels/', 
  'checking_SNPs/', 'pb_SNPanalysis_funcs.r', sep = '')
source(func_file)

###############

## LOAD DATA

# load PacBio allele counts from mpileup
#t14_test_counts_in <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
#  'snps_v2/sujan_092519/', 'tree14_v2SNPs_PB_Chr01_allele_counts.txt', 
#  sep = '')
t14_test_counts_in <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree14_v2SNPs_PB_full_allele_counts.txt', 
  sep = '')
t14_test_counts_0 <- read.table(t14_test_counts_in, header = F, sep = '\t',
  stringsAsFactors = F)

# load overall PacBio-based VCF generated my mpileup
#t14_test_vcf_in <- paste('/home/grabowsky_scratch/poplar_branch_files/',
#  'snps_v2/sujan_092519/', 'tree14_v2SNPs_PB.vcf',
#  sep = '')
t14_test_vcf_in <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree14_v2SNPs_PB_full.vcf',
  sep = '')
t14_test_vcf_0 <- read.table(t14_test_vcf_in, stringsAsFactors = F)
# Need to remove INDELs from the VCF - do that below

# load Illunina SNP info based on file from Sujan
t14_Illum_snp_info_file <- paste('/home/grabowsky_scratch/', 
  'poplar_branch_files/snps_v2/sujan_092519/', 'tree14_pval_file_SNPinfo.txt',
  sep = '')
t14_Illum_snp_info <- read.table(t14_Illum_snp_info_file, header = F,
  stringsAsFactors = F)

# load Illumina-based VCF
t14_Illum_vcf_file <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree14.good_positions.v1.vcf', sep = '')
t14_Illum_vcf <- read.table(t14_Illum_vcf_file, stringsAsFactors = F)

# Load the Illumina-based genotype dosages
t14_dosage_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree14_good_positions_v1.dosages', sep = '')
t14_dosages <- read.table(t14_dosage_file, header = T, sep = '\t', 
  stringsAsFactors = F)

# load the sample order
samp_ord_in <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree14_PB_sample_order.txt',
  sep = '')
samp_ord <- unlist(
  read.table(samp_ord_in, header = F, stringsAsFactors = F))

# load info about mapping library to sample names
lib_map_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pop_branch_pb_lib_names.txt'
lib_map <- read.table(lib_map_in, header = T, sep = '\t', stringsAsFactors = F)

#################################
## ANALYSIS SCRIPT

# Remove INDELS and Scaffold SNPs from mpileup VCF and add snp_name
indel_inds <- grep('INDEL', t14_test_vcf_0[,8])
scaff_inds <- grep('scaff', t14_test_vcf_0[,1])
remove_inds <- union(indel_inds, scaff_inds)
t14_test_vcf <- t14_test_vcf_0[-remove_inds,]
t14_test_vcf$snp_name <- paste(t14_test_vcf[ ,1], t14_test_vcf[ ,2], sep = '_')

# Remove INDELS and Scaffold SNPs from the count matrix
t14_test_counts <- t14_test_counts_0[-remove_inds, ]

# Add the Illumina Ref and Alt alleles to the PB mpileup VCF
t14_Illum_snp_info$snp_name <- paste(t14_Illum_snp_info[,1], 
  t14_Illum_snp_info[,2], sep = '_')

illum_inds <- which(t14_Illum_snp_info$snp_name %in% t14_test_vcf$snp_name)

t14_Illum_sub_1 <- t14_Illum_snp_info[illum_inds, ]
illum_ord_1 <- order(t14_Illum_sub_1[,2])
illum_ord_2 <- order(t14_Illum_sub_1[illum_ord_1,1])

t14_Illum_sub_2 <- t14_Illum_sub_1[illum_ord_1[illum_ord_2], ]

# Check that the names match up, and if so, add Illumina Ref/Alt data to VCF
if( sum(t14_Illum_sub_2$snp_name == t14_test_vcf$snp_name) == 
 nrow(t14_test_vcf) ){
  t14_test_vcf$Illum_Ref <- t14_Illum_sub_2[,3]
  t14_test_vcf$Illum_Alt <- t14_Illum_sub_2[,4]
  } else {print('Illumina and PacBio snp_names do not match')}
#########

# Generate and all allele counts to the VCF

count_split_list <- apply(t14_test_counts, 1, function(x) 
  strsplit(x, split = ','))
count_mat_list <- lapply(count_split_list, function(x) 
  matrix(data = as.numeric(unlist(x)), nrow = 8, byrow = T))
count_sum_list <- lapply(count_mat_list, function(x) apply(x, 2, sum))
count_sum_string <- unlist(lapply(count_sum_list, function(x) 
  paste(x, sep = ':', collapse = ':')))

t14_test_vcf$PB_tot_al_counts <- count_sum_string

#############

# Link library to sample names so can add branch names to data matrices
lib_map$b_name <- gsub('.', '_', as.character(lib_map$branch_name), fixed = T)

l_ord <- c()
for(i in seq(length(samp_ord))){
  tmp_ind <- which(lib_map$lib_name == samp_ord[i])
  l_ord <- c(l_ord, tmp_ind)
}

b_ord <- lib_map$b_name[l_ord]

##########

# Find the number of alleles that pass the % threshold for each SNP 
t14_test_mat1 <- make_per_mat(counts_df = t14_test_counts, per_cut = 0.1)
colnames(t14_test_mat1) <- b_ord
rownames(t14_test_mat1) <- t14_test_vcf$snp_name

# Calculate the seq depth for each sample and SNP and 
##  find those that don't pass
t14_test_depth <- make_tot_depth_mat(counts_df = t14_test_counts)

min_depth <- 20

low_d_inds <- which(t14_test_depth < min_depth)

t14_test_mat2 <- t14_test_mat1
t14_test_mat2[low_d_inds] <- NA

na_counts <- apply(t14_test_mat2, 1, function(x) sum(is.na(x)))
t14_cols <- c(5:8)
t13_cols <- c(1:4)
t14_na_counts <- apply(t14_test_mat2[,t14_cols], 1, function(x) sum(is.na(x)))
t13_na_counts <- apply(t14_test_mat2[,t13_cols], 1, function(x) sum(is.na(x)))

t14_test_vcf$t14_PB_NAs <- t14_na_counts
t14_test_vcf$t13_PB_NAs <- t13_na_counts

############

# Make strings of genotypes for each tree

t14_geno_mat <- matrix(NA, nrow = nrow(t14_test_mat2), 
  ncol = ncol(t14_test_mat2))
t14_geno_mat[which(is.na(t14_test_mat2))] <- 'N'
for(i in seq(4)){
  tmp_inds <- which(t14_test_mat2 == i)
  t14_geno_mat[tmp_inds] <- as.character(i)
}

t14_geno_string <- apply(t14_geno_mat[ , c(8,7,6,5)], 1, function(x) 
  paste(x, collapse = ':'))
t13_geno_string <- apply(t14_geno_mat[ , c(4,3,2,1)], 1, function(x) 
  paste(x, collapse = ':'))

t14_test_vcf$tree14_PB_genos <- t14_geno_string
t14_test_vcf$tree13_PB_genos <- t13_geno_string

#############

# Add Illumina-based genotypes

t14_test_vcf$tree14_Illum_genos_1 <- t14_Illum_sub_2[,5]
tmp_genos_2 <- unlist(lapply(strsplit(t14_Illum_sub_2[,5], split = ''), 
  function(x) paste(x, collapse = ':')))
tmp_genos_2 <- gsub('R|V', '1', tmp_genos_2)
tmp_genos_2 <- gsub('H', '2', tmp_genos_2)

t14_test_vcf$tree14_Illum_genos_2 <- tmp_genos_2

###############

# Add Illumina-based Dosage genotypes
## need to order the SNPs so they sync up with the PB-based files
t14_dosages$snp_name <- paste(t14_dosages[,1], t14_dosages[,2], sep = '_') 

dose_ord_1 <- order(t14_dosages$pos)
dose_ord_2 <- order(t14_dosages$chr[dose_ord_1])

t14_dosages_1 <- t14_dosages[dose_ord_1[dose_ord_2], ]

## assign dosage-based discrete genotypes
t14_dosage_mat <- matrix(unlist(t14_dosages_1[ , c(3:6)]), ncol = 4, 
  nrow = nrow(t14_dosages_1), byrow = F)

t14_dosage_genos <- assign_geno_from_dosage(t14_dosage_mat, dose_dist_cut = 0.1)

# branch order needs to be reverse for the vcf
t14_dosages_1$discrete_dose_pattern_1 <- apply(t14_dosage_genos, 1, 
  function(x) paste(rev(x), collapse = ''))

t14_dosages_1$discrete_dose_pattern_2 <- gsub('R|V', '1', 
  t14_dosages_1$discrete_dose_pattern_1)
t14_dosages_1$discrete_dose_pattern_2 <- gsub('H', '2', t14_dosages_1$discrete_dose_pattern_2)

t14_dosages_1$discrete_dose_pattern_2 <- unlist(lapply(
  strsplit(t14_dosages_1$discrete_dose_pattern_2, split = ''), function(x) 
  paste(x, collapse = ':')
  ))

dose_inds <- which(t14_dosages_1$snp_name %in% t14_test_vcf$snp_name)

t14_test_vcf$t14_dose_pattern_1 <- (
  t14_dosages_1$discrete_dose_pattern_1[dose_inds])
t14_test_vcf$t14_dose_pattern_2 <- (   
  t14_dosages_1$discrete_dose_pattern_2[dose_inds])

########################

# Test the aggreement between PacBio genotypes and Illumina genotypes
#   both VarScan and discrete,dosage-based genotypes from Illumina data

sum(t14_test_vcf$tree14_PB_genos == t14_test_vcf$tree14_Illum_genos_2)
# 265

sum(t14_test_vcf$tree14_PB_genos == t14_test_vcf$t14_dose_pattern_2)
# 238

n_perms <- 10000

illum_geno_perm_vec <- rep(NA, times = n_perms)
for(i in seq(n_perms)){
  tmp_match <- sum(t14_test_vcf$tree14_PB_genos == 
    sample(t14_test_vcf$tree14_Illum_genos_2, size = nrow(t14_test_vcf)))
  illum_geno_perm_vec[i] <- tmp_match
}

summary(illum_geno_perm_vec)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  182.0   221.0   230.0   230.4   240.0   285.0

sum(illum_geno_perm_vec >= 265)
# 86

# pval for getting 265 = 86/10000 = 0.0086
# ie: expect to get 265+ matches from random = 0.86% of time
# if use mean value, then would say that ~35 matches happen more than expected


dose_geno_perm_vec <- rep(NA, times = n_perms)
for(j in seq(n_perms)){
  tmp_match <- sum(t14_test_vcf$tree14_PB_genos == 
    sample(t14_test_vcf$t14_dose_pattern_2, size = nrow(t14_test_vcf)))
  dose_geno_perm_vec[j] <- tmp_match
}

summary(dose_geno_perm_vec)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  128.0   173.0   182.0   181.6   190.0   245.0

sum(dose_geno_perm_vec >= 238)
# 1

# p-val = 1/10000 = 1e-4

# Compare: A) Tree14_PB_genos with no missing data vs B) Illumina/Dosage
#  genotypes with no missiong data

pb_full_inds <- which(t14_test_vcf$t14_PB_NAs == 0)
# 16,982

sum(t14_test_vcf$tree14_PB_genos[pb_full_inds] == 
  t14_test_vcf$tree14_Illum_genos_2[pb_full_inds])
# 265

n_perms <- 10000

sub_illum_perm_vec <- rep(NA, times = n_perms)
for(i in seq(n_perms)){
  tmp_match <- sum(t14_test_vcf$tree14_PB_genos[pb_full_inds] ==
    sample(t14_test_vcf$tree14_Illum_genos_2[pb_full_inds], 
    size = length(pb_full_inds)))
  sub_illum_perm_vec[i] <- tmp_match
}

summary(sub_illum_perm_vec)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  178.0   222.0   232.0   232.1   242.0   283.0

sum(sub_illum_perm_vec > 265)
# 99

dose_na_inds <- grep('N', t14_test_vcf$t14_dose_pattern_2)
# 5113
dose_full_inds <- setdiff(seq(nrow(t14_test_vcf)), dose_na_inds)
# 19,665

pb_dose_full_inds <- intersect(pb_full_inds, dose_full_inds)
# 13,684

sum(t14_test_vcf$tree14_PB_genos[pb_dose_full_inds] == 
  t14_test_vcf$t14_dose_pattern_2[pb_dose_full_inds])
# 236

sub_dose_perm_vec <- rep(NA, times = n_perms)
for(j in seq(n_perms)){
  tmp_match <- sum(t14_test_vcf$tree14_PB_genos[pb_dose_full_inds] ==
    sample(t14_test_vcf$t14_dose_pattern_2[pb_dose_full_inds], 
    size = length(pb_dose_full_inds)))
  sub_dose_perm_vec[j] <- tmp_match
}

summary(sub_dose_perm_vec)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  131.0   174.0   183.0   182.9   191.0   233.0

sum(sub_dose_perm_vec >= 236)
# 0

########################

# Look at SNPs where PB and Dosage genotypes match

pb_full_inds <- which(t14_test_vcf$t14_PB_NAs == 0)

dose_na_inds <- grep('N', t14_test_vcf$t14_dose_pattern_2)
dose_full_inds <- setdiff(seq(nrow(t14_test_vcf)), dose_na_inds)

pb_dose_full_inds <- intersect(pb_full_inds, dose_full_inds) 

match_full_inds <- intersect(pb_dose_full_inds, 
  which(t14_test_vcf$tree14_PB_genos == t14_test_vcf$t14_dose_pattern_2))

# Remove SNPs that are HET in Tree13
t13_2_inds <- grep('2', t14_test_vcf$tree13_PB_genos)
# 5999 with 0.2
# 10912 with 0.1

pb_dose_no13Het_inds <- setdiff(pb_dose_full_inds, t13_2_inds)
# 9327 with 0.2
# 6604 with 0.1

#t13_N_inds <- grep('N', t14_test_vcf$tree13_PB_genos)
# 10803

# length(setdiff(t13_N_inds, t13_2_inds))
# 9378

match_full_inds_2 <- setdiff(match_full_inds, t13_2_inds)
#83 with 0.2
# 124 with 0.1

table(t14_test_vcf$t14_dose_pattern_1[match_full_inds_2])
# 0.2: HHRH HHRR HHVH HRRH HRRR HVVV RHHH RHRH RHRR RRHH RRHR RRRH VVHV 
# 0.2:    3    1    1    1   11    4    2    2   22    3   18   14    1 

# 0.1: HHHR HHRH HHRR HRHH HRHR HRRH HRRR RHHH RHHR RHRH RHRR RRHH RRHR RRRH 
# 0.1:    1    6    1    3    2    2   13    5    1    2   37    6   22   23

sum(t14_test_vcf$tree14_PB_genos[pb_dose_no13Het_inds] == 
  t14_test_vcf$t14_dose_pattern_2[pb_dose_no13Het_inds])
# 83

n_perms <- 10000
filt_match_perm_vec <- rep(NA, times = n_perms)
for(j in seq(n_perms)){
  tmp_match <- sum(t14_test_vcf$tree14_PB_genos[pb_dose_no13Het_inds] ==
    sample(t14_test_vcf$t14_dose_pattern_2[pb_dose_no13Het_inds],
    size = length(pb_dose_no13Het_inds)))
  filt_match_perm_vec[j] <- tmp_match
}

summary(filt_match_perm_vec)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  13.00   26.00   29.00   29.28   33.00   53.00
# none are close to 83


################

# Look at Branch patterns with different levels of filtering

sing_het_string <- '2:1:1:1|1:2:1:1|1:1:2:1|1:1:1:2'
sing_ref_string <- '2:2:2:1|2:2:1:2|2:1:2:2|1:2:2:2'
no_het_string <- '1:1:1:1'
no_ref_string <- '2:2:2:2'

## Full Illumina Dosage-based calls
dose_full_sing_het <- grep(sing_het_string, t14_test_vcf$t14_dose_pattern_2)
length(dose_full_sing_het)
# 3010 single het call
length(dose_full_sing_het)/length(dose_full_inds)
# [1] 0.1530638 - 15.3%

dose_full_sing_ref <- grep(sing_ref_string, t14_test_vcf$t14_dose_pattern_2)
length(dose_full_sing_ref)
# 9996 single ref calls
length(dose_full_sing_ref)/length(dose_full_inds)
# [1] 0.5083143 - 50.8% 

dose_full_no_het <- grep(no_het_string, t14_test_vcf$t14_dose_pattern_2)
length(dose_full_no_het)
# 0

dose_full_no_ref <- grep(no_ref_string, t14_test_vcf$t14_dose_pattern_2)
length(dose_full_no_ref)
# 0

## Full PacBio genotypes

pb_full_no_het <- grep(no_het_string, t14_test_vcf$tree14_PB_genos)
length(pb_full_no_het)
# 11398
length(pb_full_no_het)/length(pb_full_inds)
# [1] 0.6711813

pb_full_no_ref <- grep(no_ref_string, t14_test_vcf$tree14_PB_genos)
length(pb_full_no_ref)
# 2504
length(pb_full_no_ref)/length(pb_full_inds)
# 14.7%

pb_full_sing_het <- grep(sing_het_string, t14_test_vcf$tree14_PB_genos)
length(pb_full_sing_het)
# 1175
length(pb_full_sing_het)/length(pb_full_inds)
# 6.9%

pb_full_sing_ref <- grep(sing_ref_string, t14_test_vcf$tree14_PB_genos)
length(pb_full_sing_ref)
# 1018
length(pb_full_sing_ref)/length(pb_full_inds)
# 6.0%

## Illumina Dosage without Tree13 PB-based HETs
dose_no13Het_inds <- setdiff(dose_full_inds, t13_2_inds)
length(dose_no13Het_inds)
# 14802
length(dose_no13Het_inds)/length(dose_full_inds)
# 75.3%

length(intersect(dose_full_sing_het, dose_no13Het_inds))
# 2515
length(intersect(dose_full_sing_het, 
  dose_no13Het_inds))/length(dose_no13Het_inds)
# 17.0%

length(intersect(dose_full_sing_ref, dose_no13Het_inds))
# 7158
length(intersect(dose_full_sing_ref, 
  dose_no13Het_inds))/length(dose_no13Het_inds)
# 48.4%

## PacBio Genotypes without Tree13 PB-based HETs
pb_no13Het_inds <- setdiff(pb_full_inds, t13_2_inds)
length(pb_no13Het_inds)
# 11,653
length(pb_no13Het_inds)/length(pb_full_inds)
# 68.6%

length(intersect(pb_full_no_het, pb_no13Het_inds))
# 10909 
length(intersect(pb_full_no_het, pb_no13Het_inds))/length(pb_no13Het_inds)
# 93.6%

length(intersect(pb_full_no_ref, pb_no13Het_inds))
# 5
length(intersect(pb_full_no_ref, pb_no13Het_inds))/length(pb_no13Het_inds)
# 0.04%

length(intersect(pb_full_sing_het, pb_no13Het_inds))
# 558
length(intersect(pb_full_sing_het, pb_no13Het_inds))/length(pb_no13Het_inds)
# 4.8%

length(intersect(pb_full_sing_ref, pb_no13Het_inds))
# 31
length(intersect(pb_full_sing_ref, pb_no13Het_inds))/length(pb_no13Het_inds)
# 0.26%

## SNPs with A) Dosages B) PB genotypes C) no Tree13 HETs
pb_dose_n13Het_inds <- intersect(pb_no13Het_inds, dose_no13Het_inds)
length(pb_dose_n13Het_inds)
# [1] 9327

### Illumina Dosage
length(intersect(dose_full_sing_het, pb_dose_n13Het_inds))
# 1662
length(intersect(dose_full_sing_het, 
  pb_dose_n13Het_inds))/length(pb_dose_n13Het_inds)
# 17.8%

length(intersect(dose_full_sing_ref, pb_dose_n13Het_inds))
# 4449
length(intersect(dose_full_sing_ref, 
  pb_dose_n13Het_inds))/length(pb_dose_n13Het_inds)
# 47.7%

### PB genotypes
length(intersect(pb_full_no_het, pb_dose_n13Het_inds))
# 8719
length(intersect(pb_full_no_het, 
  pb_dose_n13Het_inds))/length(pb_dose_n13Het_inds)
# 93.4%

length(intersect(pb_full_no_ref, pb_dose_n13Het_inds))
# 5
length(intersect(pb_full_no_ref, 
  pb_dose_n13Het_inds))/length(pb_dose_n13Het_inds)
# 0.05%

length(intersect(pb_full_sing_het, pb_dose_n13Het_inds))
# 461
length(intersect(pb_full_sing_het, 
  pb_dose_n13Het_inds))/length(pb_dose_n13Het_inds)
# 4.9%

length(intersect(pb_full_sing_ref, pb_dose_n13Het_inds))
# 21
length(intersect(pb_full_sing_ref, 
  pb_dose_n13Het_inds))/length(pb_dose_n13Het_inds)
# 0.23%

## Overlap between Dosage and PacBio
match_inds <- intersect(pb_dose_n13Het_inds, 
  which(t14_test_vcf$t14_dose_pattern_2 == t14_test_vcf$tree14_PB_genos))

length(intersect(dose_full_sing_het, match_inds))
# 70
length(intersect(dose_full_sing_het, match_inds))/length(match_inds)
# 84.3%

length(intersect(dose_full_sing_ref, match_inds))
# 6
length(intersect(dose_full_sing_ref, match_inds))/length(match_inds)
# 7.2%

table(t14_test_vcf$t14_dose_pattern_1[match_inds])
# HHRH HHRR HHVH HRRH HRRR HVVV RHHH RHRH RHRR RRHH RRHR RRRH VVHV 
#    3    1    1    1   11    4    2    2   22    3   18   14    1 

sing_het_full <- grep('2:1:1:1|1:2:1:1|1:1:2:1|1:1:1:2', 
  t14_test_vcf$tree14_PB_genos)

length(sing_het_full)/length(pb_full_inds)
# 6.9% show 1-H:3-R

sing_ref_full <- grep('2:2:2:1|2:2:1:2|2:1:2:2|1:2:2:2', 
  t14_test_vcf$tree14_PB_genos)

length(sing_ref_full)/length(pb_full_inds)
# 5.9% show 3-H:1-R

pb_no_het_inds <- grep('1:1:1:1', t14_test_vcf$tree14_PB_genos)

length(pb_no_het_inds)/length(pb_full_inds)
# 67.1% show 4-R

pb_all_het_inds <- grep('2:2:2:2', t14_test_vcf$tree14_PB_genos)

length(pb_all_het_inds)/length(pb_full_inds)
# 14.7%

length(c(pb_no_het_inds, pb_all_het_inds))/length(pb_full_inds)
# 81.8% of PacBio genotypes show no variation (though may be missing some 
#  singleton-hets in these data

###########################

# Selecting SNPs to visually inspect with IGV

# 3H-1R in Illumina, 4R with PacBio
tmp_inds <- sample(x = intersect(
  which(t14_test_vcf$t14_dose_pattern_1 == 'RHHH'),
  which(t14_test_vcf$tree14_PB_genos == '1:1:1:1')),
  size = 2
)
t14_test_vcf[tmp_inds, c('V1', 'V2', 't14_dose_pattern_1', 'tree14_PB_genos')]

# 3R-1H in Illumina, 4R in PacBio
tmp_inds <- sample(x = intersect(
  which(t14_test_vcf$t14_dose_pattern_1 == 'HRRR'),
  which(t14_test_vcf$tree14_PB_genos == '1:1:1:1')),
  size = 2
)
t14_test_vcf[tmp_inds, c('V1', 'V2', 't14_dose_pattern_1', 'tree14_PB_genos')]

# Overlap of 3H-1R between Illumina and PacBio
tmp_inds <- sample(x = intersect(match_full_inds_2, 
  which(t14_test_vcf$t14_dose_pattern_1 == 'HHHR')),
  size = 2)
t14_test_vcf[tmp_inds, c('V1', 'V2', 't14_dose_pattern_1', 'tree14_PB_genos')]



