# Analysis for looking at private, shared, etc SVs in the PacBio SV data

function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)

data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'
combo_file <- 'ref.ALLData.vcf'
combo_file_tot <- paste(data_dir, combo_file, sep = '')
combo_df <- make_combo_indel_df(combo_file_tot)

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')

## Remove bad branches
bad_branches <- c('13.4', '14.1')
bad_samps <- c()
for(bb in bad_branches){
  tmp_ind <- which(as.character(samp_meta$branch_name) == bb)
  bad_samps <- c(bad_samps, samp_meta$lib_name[tmp_ind])
}

bad_vcf_cols <- c()
for(bs in bad_samps){
  tmp_col_ind <- which(colnames(combo_df) == bs)
  bad_vcf_cols <- c(bad_vcf_cols, tmp_col_ind)
}

combo_df_filt <- combo_df[ , -sort(bad_vcf_cols)]

## find and remove any SVs with missing genotypes in any sample
first_samp_ind <- which(colnames(combo_df_filt) == 'FORMAT') + 1
last_samp_ind <- which(colnames(combo_df_filt) == 'full_name') - 1

miss_geno_inds <- c()
for(sn in c(first_samp_ind:last_samp_ind)){
  tmp_inds <- grep('.', combo_df_filt[,sn], fixed = T)
  miss_geno_inds <- c(miss_geno_inds, tmp_inds)
}

miss_inds_toremove <- unique(miss_geno_inds)

combo_df_2 <- combo_df_filt[-sort(miss_inds_toremove), ]

# Extract genotype information. Inclueds
# converting genos to numeric
# getting coverage
# getting number of reads supporting each allele
combo_df_2[c(1:10), first_samp_ind]

geno_info_list <- list()
for(s_ind in c(first_samp_ind:last_samp_ind)){
  samp_name <- colnames(combo_df_2)[s_ind]
  tmp_geno_list <- strsplit(combo_df_2[ , s_ind], split = ':')
  tmp_coverage <- as.numeric(unlist(lapply(tmp_geno_list, function(x) x[3])))
  tmp_num_genos <- unlist(lapply(tmp_geno_list, function(x) 
    sum(as.numeric(unlist(strsplit(x[1], split = '/'))))))
  tmp_allele_reads <- matrix(
    data = unlist(lapply(tmp_geno_list, function(x) 
      as.numeric(unlist(strsplit(x[2], split = ','))))), 
    byrow = T, ncol = 2, 
    dimnames = list(rownames = NULL, colnames = c('n_REF', 'n_ALT')))
  geno_info_list[[samp_name]] <- list()
  geno_info_list[[samp_name]][[1]] <- tmp_num_genos
  geno_info_list[[samp_name]][[2]] <- tmp_coverage
  geno_info_list[[samp_name]][[3]] <- tmp_allele_reads
}

# generate matrix of coverage
geno_cov_mat <- matrix(
  data = unlist(lapply(geno_info_list, function(x) x[[2]])),
  ncol = length(geno_info_list),
  byrow = F)

min_cov_vec <- apply(geno_cov_mat, 1, min)

min_cov_floor_vec <- c()
for(mc in seq(5, 30, 5)){
  tmp_val <- sum(min_cov_vec >= mc)
  min_cov_floor_vec <- c(min_cov_floor_vec , tmp_val)
}
names(min_cov_floor_vec) <- seq(5, 30, 5)

min_cov_floor_vec
#     5    10    15    20    25    30 
# 56642 54052 40021 10999  2149  1009 

# only 1009 SV's have a minimum coverage of 30; 54052 have min coverage of 10

# RULES I'M GOING BY AS OF 11/2 - these can change, but I'm going to use these
#  for filtering
# 1) Minimum coverage of 10 for every sample
# 2) Het calls must have at least 2 reads of each allele
# 3) Het calls must have at least a 0.15 allele ratio - may want to be more
#      conservative with this - maybe 0.2 or 0.25
# 4) Homozygous call must have at lower than 0.05 allele ratio - this may
#      need to be adjusted

# I'll try calling my own genotypes

test_ref_ratio <- apply(geno_info_list[[1]][[3]], 1, function(x) x[1]/sum(x))
test_genos <- rep(NA, times = length(test_ref_ratio))
test_call_hom_ref <- which(test_ref_ratio < 0.15)
test_call_hom_alt <- which(test_ref_ratio > 0.85)
test_call_het <- setdiff(seq(length(test_genos)), 
  union(test_call_hom_ref, test_call_hom_alt))
test_genos[test_call_hom_ref] <- 0
test_genos[test_call_hom_alt] <- 2
test_genos[test_call_het] <- 1

test_too_low_inds <- which(geno_info_list[[1]][[2]] < 10)
test_min_ratio <- apply(geno_info_list[[1]][[3]], 1, function(x) min(x)/sum(x))
test_hom_inds <- setdiff(seq(length(test_genos)), test_call_het)
test_too_skewed_homs <- test_hom_inds[
  which(test_min_ratio[test_hom_inds] > 0.05)]

test_na_inds <- union(test_too_low_inds, test_too_skewed_homs)
test_genos[test_na_inds] <- NA

# NEXT: Turn the above into a function or loop to generate genotypes for
#  each sample
#####################

too_low_inds <- which(min_cov_vec < 10)

test_het_inds <- which(geno_info_list[[1]][[1]] == 1)
test_min_allele_read <- apply(geno_info_list[[1]][[3]], 1, min)
too_low_hets <- test_het_inds[which(test_min_allele_read[test_het_inds] < 2)]

test_ratio <- apply(geno_info_list[[1]][[3]], 1, function(x) min(x)/sum(x))
too_skewed_hets <- test_het_inds[which(test_ratio[test_het_inds] < 0.15)]

test_hom_inds <- setdiff(seq(length(geno_info_list[[1]][[1]])), test_het_inds)
too_skewed_homs <- test_hom_inds[which(test_ratio[test_hom_inds] > 0.02)]

min_hom_cov <- 20

zero_inds <- which(geno_info_list[[1]][[1]] == 0)
two_inds <- which(geno_info_list[[1]][[1]] == 2)
homozyg_inds <- sort(union(zero_inds, two_inds))
homozyg_cov <- geno_info_list[[1]][[2]][homozyg_inds]
low_cov_hom_inds <- homozyg_inds[which(homozyg_cov < min_hom_cov)]

# there are a LOT of low-ish coverage SVs with homozygous genotypes. I want
#   to model the allele ratios in the heterozygotes to try to model the 
#   accuracy of the genotypes

test_ratio <- apply(geno_info_list[[1]][[3]], 1, function(x) x[2]/sum(x))
intersect(which(test_ratio[homozyg_inds] != 0), which(test_ratio[homozyg_inds] != 1))

# 1) Identify het SVs
# 2) Calculate % alt alleles in het calls
# 3) For each SV, aggregate the % alt alleles for each het call
# 4) Set desired min prob of homozygous genotype, ex: 0.001
# 5) For each SV, determine mean or min/max %alt alleles seen in het genotypes
# 5a)  depends on if is 0 (0/0) or 2 (1/1) genotype
# 6) Determine prob. of homozygous genotype given the prob of %alt alleles
#      in the het samples

het_inds <- which(geno_info_list[[1]][[1]] == 1)
test_ratio <- apply(geno_info_list[[1]][[3]], 1, function(x) x[2]/sum(x))
het_ratio <- rep(NA, times = length(geno_info_list[[1]][[1]]))
het_ratio[het_inds] <- test_ratio[het_inds]

# some rules for heterozygous genotypes: ratio bust be between 0.05

# let's see what the "bad" het inds are in each sample - these are hets called
#   using very skewed read counts


test_geno_list <- strsplit(combo_df_2[ , first_samp_ind], split = ':')
test_coverage <- as.numeric(unlist(lapply(test_geno_list, function(x) x[3])))
test_geno_calls <- unlist(lapply(test_geno_list, function(x) sum(as.numeric(unlist(strsplit(x[1], split = '/'))))))
test_allele_reads <- matrix(data = unlist(lapply(test_geno_list, function(x) as.numeric(unlist(strsplit(x[2], split = ','))))), byrow = T, ncol = 2, dimnames = list(rownames = NULL, colnames = c('n_REF', 'n_ALT')))

