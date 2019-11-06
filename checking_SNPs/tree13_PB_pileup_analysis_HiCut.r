# Script for analyzing the mpileup results for Tree14

## LOAD FUNCTIONS AND PACKAGES
func_file <- paste('/home/grabowsky/tools/workflows/poplar_branch_indels/', 
  'checking_SNPs/', 'pb_SNPanalysis_funcs.r', sep = '')
source(func_file)

###############

## LOAD DATA

t13_counts_in <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree13_v2SNPs_PB_full_allele_counts.txt', 
  sep = '')
t13_counts_0 <- read.table(t13_counts_in, header = F, sep = '\t',
  stringsAsFactors = F)

t13_vcf_in <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_v2SNPs_PB_full.vcf',
  sep = '')
t13_vcf_0 <- read.table(t13_vcf_in, stringsAsFactors = F)
# Need to remove INDELs from the VCF - do that below

# load Illunina SNP info based on file from Sujan
t13_Illum_snp_info_file <- paste('/home/grabowsky_scratch/', 
  'poplar_branch_files/snps_v2/sujan_092519/', 'tree13_pval_file_SNPinfo.txt',
  sep = '')
t13_Illum_snp_info <- read.table(t13_Illum_snp_info_file, header = F,
  stringsAsFactors = F)

# load Illumina-based VCF
t13_Illum_vcf_file <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree13.good_positions.v1.vcf', sep = '')
t13_Illum_vcf <- read.table(t13_Illum_vcf_file, stringsAsFactors = F)

# Load the Illumina-based genotype dosages
t13_dosage_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_good_positions_v1.dosages', sep = '')
t13_dosages <- read.table(t13_dosage_file, header = T, sep = '\t', 
  stringsAsFactors = F)

# load the sample order
## the Tree13 and Tree14 PB pileup files have the same library orders
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
indel_inds <- grep('INDEL', t13_vcf_0[,8])
scaff_inds <- grep('scaff', t13_vcf_0[,1])
remove_inds <- union(indel_inds, scaff_inds)
t13_vcf <- t13_vcf_0[-remove_inds,]
t13_vcf$snp_name <- paste(t13_vcf[ ,1], t13_vcf[ ,2], sep = '_')

# Remove INDELS and Scaffold SNPs from the count matrix
t13_counts <- t13_counts_0[-remove_inds, ]

# Add the Illumina Ref and Alt alleles to the PB mpileup VCF
t13_Illum_snp_info$snp_name <- paste(t13_Illum_snp_info[,1], 
  t13_Illum_snp_info[,2], sep = '_')

illum_inds <- which(t13_Illum_snp_info$snp_name %in% t13_vcf$snp_name)

t13_Illum_sub_1 <- t13_Illum_snp_info[illum_inds, ]
illum_ord_1 <- order(t13_Illum_sub_1[,2])
illum_ord_2 <- order(t13_Illum_sub_1[illum_ord_1,1])

t13_Illum_sub_2 <- t13_Illum_sub_1[illum_ord_1[illum_ord_2], ]

# Check that the names match up, and if so, add Illumina Ref/Alt data to VCF
if( sum(t13_Illum_sub_2$snp_name == t13_vcf$snp_name) == 
 nrow(t13_vcf) ){
  t13_vcf$Illum_Ref <- t13_Illum_sub_2[,3]
  t13_vcf$Illum_Alt <- t13_Illum_sub_2[,4]
  } else {print('Illumina and PacBio snp_names do not match')}
#########

# Generate and all allele counts to the VCF

count_split_list <- apply(t13_counts, 1, function(x) 
  strsplit(x, split = ','))
count_mat_list <- lapply(count_split_list, function(x) 
  matrix(data = as.numeric(unlist(x)), nrow = 8, byrow = T))
count_sum_list <- lapply(count_mat_list, function(x) apply(x, 2, sum))
count_sum_string <- unlist(lapply(count_sum_list, function(x) 
  paste(x, sep = ':', collapse = ':')))

t13_vcf$PB_tot_al_counts <- count_sum_string

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
t13_mat1 <- make_per_mat(counts_df = t13_counts, per_cut = 0.1)
colnames(t13_mat1) <- b_ord
rownames(t13_mat1) <- t13_vcf$snp_name

# Calculate the seq depth for each sample and SNP and 
##  find those that don't pass
t13_depth <- make_tot_depth_mat(counts_df = t13_counts)

min_depth <- 20

low_d_inds <- which(t13_depth < min_depth)

t13_mat2 <- t13_mat1
t13_mat2[low_d_inds] <- NA

na_counts <- apply(t13_mat2, 1, function(x) sum(is.na(x)))
t14_cols <- c(5:8)
t13_cols <- c(1:4)
t14_na_counts <- apply(t13_mat2[,t14_cols], 1, function(x) sum(is.na(x)))
t13_na_counts <- apply(t13_mat2[,t13_cols], 1, function(x) sum(is.na(x)))

t13_vcf$t14_PB_NAs <- t14_na_counts
t13_vcf$t13_PB_NAs <- t13_na_counts

############

# Make strings of genotypes for each tree

t13_geno_mat <- matrix(NA, nrow = nrow(t13_mat2), 
  ncol = ncol(t13_mat2))
t13_geno_mat[which(is.na(t13_mat2))] <- 'N'
for(i in seq(4)){
  tmp_inds <- which(t13_mat2 == i)
  t13_geno_mat[tmp_inds] <- as.character(i)
}

t14_geno_string <- apply(t13_geno_mat[ , c(8,7,6,5)], 1, function(x) 
  paste(x, collapse = ':'))
t13_geno_string <- apply(t13_geno_mat[ , c(4,3,2,1)], 1, function(x) 
  paste(x, collapse = ':'))

t13_vcf$tree14_PB_genos <- t14_geno_string
t13_vcf$tree13_PB_genos <- t13_geno_string

#############

# Add Illumina-based genotypes

t13_vcf$tree13_Illum_genos_1 <- t13_Illum_sub_2[,5]
tmp_genos_2 <- unlist(lapply(strsplit(t13_Illum_sub_2[,5], split = ''), 
  function(x) paste(x, collapse = ':')))
tmp_genos_2 <- gsub('R|V', '1', tmp_genos_2)
tmp_genos_2 <- gsub('H', '2', tmp_genos_2)

t13_vcf$tree13_Illum_genos_2 <- tmp_genos_2

###############

# Add Illumina-based Dosage genotypes
## need to order the SNPs so they sync up with the PB-based files
t13_dosages$snp_name <- paste(t13_dosages[,1], t13_dosages[,2], sep = '_') 

dose_ord_1 <- order(t13_dosages$pos)
dose_ord_2 <- order(t13_dosages$chr[dose_ord_1])

t13_dosages_1 <- t13_dosages[dose_ord_1[dose_ord_2], ]

## assign dosage-based discrete genotypes
t13_dosage_mat <- matrix(unlist(t13_dosages_1[ , c(3:6)]), ncol = 4, 
  nrow = nrow(t13_dosages_1), byrow = F)

t13_dosage_genos <- assign_geno_from_dosage(t13_dosage_mat, dose_dist_cut = 0.1)

# branch order needs to be reverse for the vcf
t13_dosages_1$discrete_dose_pattern_1 <- apply(t13_dosage_genos, 1, 
  function(x) paste(rev(x), collapse = ''))

t13_dosages_1$discrete_dose_pattern_2 <- gsub('R|V', '1', 
  t13_dosages_1$discrete_dose_pattern_1)
t13_dosages_1$discrete_dose_pattern_2 <- gsub('H', '2', 
  t13_dosages_1$discrete_dose_pattern_2)

t13_dosages_1$discrete_dose_pattern_2 <- unlist(lapply(
  strsplit(t13_dosages_1$discrete_dose_pattern_2, split = ''), function(x) 
  paste(x, collapse = ':')
  ))

dose_inds <- which(t13_dosages_1$snp_name %in% t13_vcf$snp_name)

t13_vcf$t13_dose_pattern_1 <- (
  t13_dosages_1$discrete_dose_pattern_1[dose_inds])
t13_vcf$t13_dose_pattern_2 <- (   
  t13_dosages_1$discrete_dose_pattern_2[dose_inds])

#####################

# save info
t13_info_file_out <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_v2SNPs_combined_info.rds',
  sep = '')

saveRDS(t13_vcf, file = t13_info_file_out)

########################

# Test the aggreement between PacBio genotypes and Illumina genotypes
#   both VarScan and discrete,dosage-based genotypes from Illumina data

sum(t13_vcf$tree13_PB_genos == t13_vcf$tree13_Illum_genos_2)
# 286

sum(t13_vcf$tree13_PB_genos == t13_vcf$t13_dose_pattern_2)
# 267

n_perms <- 10000

illum_geno_perm_vec <- rep(NA, times = n_perms)
for(i in seq(n_perms)){
  tmp_match <- sum(t13_vcf$tree13_PB_genos == 
    sample(t13_vcf$tree13_Illum_genos_2, size = nrow(t13_vcf)))
  illum_geno_perm_vec[i] <- tmp_match
}

summary(illum_geno_perm_vec)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  163.0   199.0   208.0   208.2   217.0   257.0 

sum(illum_geno_perm_vec >= 286)
# 0

# pval for getting 286 < 1/10000 = 1e-4
# ie: expect to get 286+ matches from random less than  0.01% of time
# if use mean value, then would say that ~78 matches happen more than expected

dose_geno_perm_vec <- rep(NA, times = n_perms)
for(j in seq(n_perms)){
  tmp_match <- sum(t13_vcf$tree13_PB_genos == 
    sample(t13_vcf$t13_dose_pattern_2, size = nrow(t13_vcf)))
  dose_geno_perm_vec[j] <- tmp_match
}

summary(dose_geno_perm_vec)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  134.0   174.0   183.0   183.1   192.0   231.0 

sum(dose_geno_perm_vec >= 267)
# 0

# pval for getting 267 < 1/10000 = 1e-4
# if use mean value, then get 84 matches more than expected by chance

# Compare: A) Tree14_PB_genos with no missing data vs B) Illumina/Dosage
#  genotypes with no missiong data

pb_full_inds <- which(t13_vcf$t13_PB_NAs == 0)
# 12,504

sum(t13_vcf$tree13_PB_genos[pb_full_inds] == 
  t13_vcf$tree13_Illum_genos_2[pb_full_inds])
# 286

n_perms <- 10000

sub_illum_perm_vec <- rep(NA, times = n_perms)
for(i in seq(n_perms)){
  tmp_match <- sum(t13_vcf$tree13_PB_genos[pb_full_inds] ==
    sample(t13_vcf$tree13_Illum_genos_2[pb_full_inds], 
    size = length(pb_full_inds)))
  sub_illum_perm_vec[i] <- tmp_match
}

summary(sub_illum_perm_vec)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  160.0   201.0   209.0   209.5   218.0   262.0

sum(sub_illum_perm_vec > 267)
# 0
# 209 - 267 = 58

dose_na_inds <- grep('N', t13_vcf$t13_dose_pattern_2)
# 4367
dose_full_inds <- setdiff(seq(nrow(t13_vcf)), dose_na_inds)
# 17,150

pb_dose_full_inds <- intersect(pb_full_inds, dose_full_inds)
# 10,151

sum(t13_vcf$tree13_PB_genos[pb_dose_full_inds] == 
  t13_vcf$t13_dose_pattern_2[pb_dose_full_inds])
# 262

sub_dose_perm_vec <- rep(NA, times = n_perms)
for(j in seq(n_perms)){
  tmp_match <- sum(t13_vcf$tree13_PB_genos[pb_dose_full_inds] ==
    sample(t13_vcf$t13_dose_pattern_2[pb_dose_full_inds], 
    size = length(pb_dose_full_inds)))
  sub_dose_perm_vec[j] <- tmp_match
}

summary(sub_dose_perm_vec)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 131.0   172.0   180.0   180.1   189.0   230.0

sum(sub_dose_perm_vec >= 262)
# 0
# 262-180 = 82

########################

# Look at SNPs where PB and Dosage genotypes match

# t13_info_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
#  'snps_v2/sujan_092519/', 'tree13_v2SNPs_combined_info.rds',
#  sep = '')

# t13_vcf <- readRDS(t13_info_file)

pb_full_inds <- which(t13_vcf$t13_PB_NAs == 0)

dose_na_inds <- grep('N', t13_vcf$t13_dose_pattern_2)
dose_full_inds <- setdiff(seq(nrow(t13_vcf)), dose_na_inds)

pb_dose_full_inds <- intersect(pb_full_inds, dose_full_inds) 

match_full_inds <- intersect(pb_dose_full_inds, 
  which(t13_vcf$tree13_PB_genos == t13_vcf$t13_dose_pattern_2))

# Remove SNPs that are HET in Tree14
t14_2_inds <- grep('2', t13_vcf$tree14_PB_genos)
# 9270

pb_dose_no14Het_inds <- setdiff(pb_dose_full_inds, t14_2_inds)
# 5023

#t14_N_inds <- grep('N', t13_vcf$tree14_PB_genos)
# 6338

# length(setdiff(t14_N_inds, t14_2_inds))
# 4752

match_full_inds_2 <- setdiff(match_full_inds, t14_2_inds)
# 113 with 0.1

table(t13_vcf$t13_dose_pattern_1[match_full_inds_2])
# HHHR HHRR HRHR HRRR RHHH RHHR RHRR RRHH RRHR RRRH 
#   4    1    1   21    7    2   20    7   28   22

sum(t13_vcf$tree13_PB_genos[pb_dose_no14Het_inds] == 
  t13_vcf$t13_dose_pattern_2[pb_dose_no14Het_inds])
# 113

n_perms <- 10000
filt_match_perm_vec <- rep(NA, times = n_perms)
for(j in seq(n_perms)){
  tmp_match <- sum(t13_vcf$tree13_PB_genos[pb_dose_no14Het_inds] ==
    sample(t13_vcf$t13_dose_pattern_2[pb_dose_no14Het_inds],
    size = length(pb_dose_no14Het_inds)))
  filt_match_perm_vec[j] <- tmp_match
}

summary(filt_match_perm_vec)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  16.00   37.00   41.00   41.31   45.00   64.00 
# none are close to 113

sing_het_string <- '2:1:1:1|1:2:1:1|1:1:2:1|1:1:1:2'

sing_het_inds <- grep(sing_het_string, t13_vcf$t13_dose_pattern_2)

match_full_inds_3 <- intersect(match_full_inds_2, sing_het_inds)

#length(match_full_inds_3)
#[1] 91

match_positions_out <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree13_filtered_somatic_SNPs.positions', sep = '')

write.table(t13_vcf[match_full_inds_3, c(1:2)], 
  file = match_positions_out, quote = F, sep = '\t', row.names = F, 
  col.names = F)

# matches including the branch order
good_order_strings <- '1:2:2:2|1:1:2:2'

good_ord_inds <- grep(good_order_strings, t13_vcf$t13_dose_pattern_2)

match_full_inds_4 <- intersect(match_full_inds_2, 
  c(sing_het_inds, good_ord_inds))

length(match_full_inds_4)
# 105

table(t13_vcf$t13_dose_pattern_1[match_full_inds_4])
#  HRRR RHHH RHRR RRHH RRHR RRRH 
#    21    7   20    7   28   22 

match_4_positions_out <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_filtered_somatic_SNPs_expanded.positions', 
  sep = '')

write.table(t13_vcf[match_full_inds_4, c(1:2)],
  file = match_4_positions_out, quote = F, sep = '\t', row.names = F,
  col.names = F)


####################
# Stopped adapting the Tree14 code for Tree13 at this point
####################


