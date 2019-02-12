# Analysis looking at only deletions from the pbsv pipeline

# Outline
# Do all analyses for all 4 runs to show repeatibility
# Generate trees for all 10 samples to make sure they make sense
# If 10-samp trees make sense, then do rest on only 8 "good" samples
# Tally number of genotypes by sample
# Tally number of singleton DELs by sample
# Tally clone-specific DELs
# Does pattern of sharing make sense
 
# Find overlap in DELs across runs
# Size of deletions
# Location of deletions


# LOAD LIBRARIES #
function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(ape)

# LOAD DATA #
data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'

combo_file_name_vec <- c('ref.ALLData.vcf', 'ref.ALLData_try2.vcf', 'ref.ALLData.try3.vcf', 'ref.ALLData.try4.vcf')

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')

# OUTPUT INFO #
analysis_res_dir <- paste(data_dir, 'analysis_results/', sep = '')

## NJ tree of variable Deletions
DEL_nj_file <- paste(analysis_res_dir, 'DEL_var_noNA_NJ_tree.png', sep = '')

DEL_upgma_file <- paste(analysis_res_dir, 'DEL_var_noNA_UPGMA_tree.png', 
  sep = '')
# SET CONSTANTS #
samp_order_list <- list()
samp_order_list[[1]] <- c("PAXL", "PAXN", "PAYK", "PAYZ", "PAZF", "PAZG",
  "PAZH", "PBAT", "PBAU", "PBAW")
samp_order_list[[2]] <- c('PBAW', 'PBAU', 'PBAT', 'PAZH', 'PAZG', 'PAZF',
  'PAYZ', 'PAYK', 'PAXN', 'PAXL')
samp_order_list[[3]] <- c('PBAU', 'PAYZ', 'PBAW', 'PAXN', 'PBAT', 'PAXL',
  'PAZG', 'PAYK', 'PAZF', 'PAZH')
samp_order_list[[4]] <- c('PAYZ', 'PAZF', 'PAXN', 'PAZG', 'PAXL', 'PBAT',
  'PAYK', 'PBAW', 'PAZH', 'PBAU')

# SET VARIABLES #


##################
# Load InDel data for 4 runs
precall_filt_list <- list()
for(cfnv in seq(length(combo_file_name_vec))){
  tmp_tot_name <- paste(data_dir, combo_file_name_vec[cfnv], sep = '')
  tmp_df <- make_combo_indel_df(tmp_tot_name)
  tmp_precall_df <- precalling_SV_filtering(sv_geno_df = tmp_df,
    mer_length = 8, per_mn_cutoff = 0.7, per_bn_pure_cutoff = 0.5,
    per_bn_multi_cutoff = 0.6, dist_cut = 1000, sd_cut = 3, use_sd = T)
  precall_filt_list[[cfnv]] <- tmp_precall_df
  print(cfnv)
}

filt_genotypes_list <- lapply(precall_filt_list, function(x)
  generate_filtered_genotypes(combo_df = x, min_sv_coverage = 10,
  het_ratio_cut = 0.25,
  min_minor_allele_count = 2, max_hom_ratio = 0.05, est_error = 0.032,
  est_penetrance = 0.35, good_score = 0.9, great_score = 0.98,
  same_geno_bonus = 0.67, ss_min_cov = 20))

noNA_geno_list <- lapply(filt_genotypes_list, filt_geno_mat, max_nas = 0,
  min_length = 20)

for(ngl in seq(length(noNA_geno_list))){
  tmp_samp_order <- colnames(noNA_geno_list[[ngl]])
  tmp_branch_names <- c()
  for(so in tmp_samp_order){
    tmp_meta_ind <- which(samp_meta$lib_name == so)
    tmp_branch_names <- c(tmp_branch_names,
      samp_meta$branch_name[tmp_meta_ind])
  }
  colnames(noNA_geno_list[[ngl]]) <- tmp_branch_names
}

# Get DELetions
noNA_del_list <- list()
for(ngl in seq(length(noNA_geno_list))){
  tmp_del_inds <- grep('DEL', rownames(noNA_geno_list[[ngl]]))
  tmp_del_mat <- noNA_geno_list[[ngl]][tmp_del_inds, ]
  noNA_del_list[[ngl]] <- tmp_del_mat
}

# Only those that vary...
noNA_del_var_list <- list()
for(ndvl in seq(length(noNA_del_list))){
  tmp_n_al <- apply(noNA_del_list[[ndvl]], 1, function(x) length(unique(x)))
  tmp_var_inds <- which(tmp_n_al == 2)
  noNA_del_var_list[[ndvl]] <- noNA_del_list[[ndvl]][tmp_var_inds, ]
}

n_var_dels <- unlist(lapply(noNA_del_var_list, nrow))
# [1] 35 36 38 43

# Distance-based graphs/trees
del_nj_list <- list()
del_upgma_list <- list()

for(ngl in seq(length(noNA_del_var_list))){
  tmp_geno_t <- t(noNA_del_var_list[[ngl]])
  tmp_dist_mat <- dist(tmp_geno_t, method = 'manhattan', diag = T, upper = T)
  del_nj_list[[ngl]] <- nj(tmp_dist_mat)
  del_upgma_list[[ngl]] <- hclust(tmp_dist_mat, method = 'average')
}

png(filename = DEL_nj_file, width = 1000, height = 1000)
par(mfrow = c(2,2))
for(nji in seq(length(del_nj_list))){
  plot(del_nj_list[[nji]], main = paste('Poplar Branch DEL NJ tree\nRun', nji))
}
dev.off()

png(filename = DEL_upgma_file, width = 1000, height = 1000)
par(mfrow = c(2,2))
for(upi in seq(length(del_upgma_list))){
  plot(del_upgma_list[[upi]], 
    main = paste('Poplar Branch DEL UPGMA tree\nRun', upi))
}
dev.off()

##########
# test overlap between 4 runs
del_int_1_2 <- intersect(rownames(noNA_del_var_list[[1]]), 
  rownames(noNA_del_var_list[[2]]))
length(del_int_1_2)
# 20

del_int_1_2_diff_genos <- c()
for(di12 in del_int_1_2){
  tmp_same_genos <- sum(noNA_del_var_list[[1]][di12,] == 
    noNA_del_var_list[[2]][di12, ])
  if(tmp_same_genos != 10){
    del_int_1_2_diff_genos <- c(del_int_1_2_diff_genos, di12)
  }
  #print(tmp_same_genos)
}
del_int_1_2_diff_genos
#NULL - all genotypes are the same between the two runs
setdiff(rownames(noNA_del_var_list[[1]]), del_int_1_2)
setdiff(rownames(noNA_del_var_list[[2]]), del_int_1_2)
# 8 non-identical Deletions are likely same deletion:
# List1				# List2
same_1_2_vec <- c('Chr01_35013696_DEL_1848',	'Chr01_35013673_DEL_1852',
'Chr02_19900668_DEL_12316',	'Chr02_19900665_DEL_12317',
'Chr06_24972400_DEL_477',	'Chr06_24972400_DEL_479',
'Chr07_12613882_DEL_54',	'Chr07_12613882_DEL_53',
'Chr08_17464456_DEL_661',	'Chr08_17464456_DEL_660',
'Chr10_1245178_DEL_198',	'Chr10_1245178_DEL_199',
'Chr10_16534461_DEL_472',	'Chr10_16534461_DEL_469',
'Chr17_7812838_DEL_6636',	'Chr17_7812838_DEL_6635')

same_1_2_mat <- matrix(data = same_1_2_vec, ncol = 2, byrow = T)

del_same_1_2_diff_genos <- c()
for(ds12 in seq(nrow(same_1_2_mat))){
  tmp_same_genos <- sum(noNA_del_var_list[[1]][same_1_2_mat[ds12,1], ] == 
    noNA_del_var_list[[2]][same_1_2_mat[ds12,2], ])
  if(tmp_same_genos != 10){
    del_same_1_2_diff_genos <- c(del_same_1_2_diff_genos, ds12)
  }  
#print(tmp_same_genos)
}
del_same_1_2_diff_genos
# NULL - all genotypes are the same between the DELS that aren't 
#  identical but vary similar (and prob the same) between the two runs

##

del_int_1_3 <- intersect(rownames(noNA_del_var_list[[1]]), 
  rownames(noNA_del_var_list[[3]]))
length(del_int_1_3)
# 20
del_int_1_3_diff_genos <- c()
for(di13 in del_int_1_3){
  tmp_same_genos <- sum(noNA_del_var_list[[1]][di13,] ==
    noNA_del_var_list[[3]][di13, ])
  if(tmp_same_genos != 10){
    del_int_1_3_diff_genos <- c(del_int_1_3_diff_genos, di13)
  }
  #print(tmp_same_genos)
}
del_int_1_3_diff_genos
#[1] "Chr13_5343435_DEL_20"
# This deletion is different at 4 samples - should be removed from final
#  set
setdiff(rownames(noNA_del_var_list[[1]]), del_int_1_3)
setdiff(rownames(noNA_del_var_list[[3]]), del_int_1_3)
# 8 non-identical Deletions are likely same deletion:
# List1				# List3
same_1_3_vec <- c('Chr01_35013696_DEL_1848',	'Chr01_35013673_DEL_1852',
'Chr02_15912927_DEL_55',	'Chr02_15912927_DEL_57',
'Chr03_15549628_DEL_7640',	'Chr03_15549628_DEL_7639',
'Chr06_24972400_DEL_477',	'Chr06_24972400_DEL_478',
'Chr08_17464456_DEL_661',	'Chr08_17464456_DEL_659',
'Chr10_16534461_DEL_472',	'Chr10_16534461_DEL_471',
'Chr14_8258327_DEL_142',	'Chr14_8258322_DEL_141',
'Chr17_3368868_DEL_8785',	'Chr17_3368868_DEL_8786')

same_1_3_mat <- matrix(data = same_1_3_vec, ncol = 2, byrow = T)

del_same_1_3_diff_genos <- c()
for(ds13 in seq(nrow(same_1_3_mat))){
  tmp_same_genos <- sum(noNA_del_var_list[[1]][same_1_3_mat[ds13,1], ] ==
    noNA_del_var_list[[3]][same_1_3_mat[ds13,2], ])
  if(tmp_same_genos != 10){
    del_same_1_3_diff_genos <- c(del_same_1_3_diff_genos, ds13)
  }
#print(tmp_same_genos)
}
del_same_1_3_diff_genos
# NULL - all genotypes are the same between the DELS that aren't 
#  identical but vary similar (and prob the same) between the two runs

##

del_int_1_4 <- intersect(rownames(noNA_del_var_list[[1]]),
  rownames(noNA_del_var_list[[4]]))
length(del_int_1_4)
# 23
del_int_1_4_diff_genos <- c()
for(di14 in del_int_1_4){
  tmp_same_genos <- sum(noNA_del_var_list[[1]][di14,] ==
    noNA_del_var_list[[4]][di14, ])
  if(tmp_same_genos != 10){
    del_int_1_4_diff_genos <- c(del_int_1_4_diff_genos, di14)
  }
  #print(tmp_same_genos)
}
del_int_1_4_diff_genos
# [1] "Chr13_5343435_DEL_20"
# This deletion is different at 2 samples - is same weird deletion as 1v3,
#  should be removed from final list
setdiff(rownames(noNA_del_var_list[[1]]), del_int_1_4)
setdiff(rownames(noNA_del_var_list[[4]]), del_int_1_4)
# 8 non-identical Deletions are likely same deletion:
# List1				# List4
same_1_4_vec <- c('Chr01_35013696_DEL_1848',	'Chr01_35013673_DEL_1852',
'Chr03_15549628_DEL_7640',	'Chr03_15549628_DEL_7639',
'Chr06_9510637_DEL_294',	'Chr06_9510636_DEL_295',
'Chr06_23906913_DEL_10372',	'Chr06_23906913_DEL_10373',
'Chr08_17464456_DEL_661',	'Chr08_17464456_DEL_659',
'Chr10_16534461_DEL_472',	'Chr10_16534461_DEL_471',
'Chr17_3368868_DEL_8785',	'Chr17_3368868_DEL_8784',
'Chr17_13034495_DEL_321',	'Chr17_13034495_DEL_323')

same_1_4_mat <- matrix(data = same_1_4_vec, ncol = 2, byrow = T)

del_same_1_4_diff_genos <- c()
for(ds14 in seq(nrow(same_1_4_mat))){
  tmp_same_genos <- sum(noNA_del_var_list[[1]][same_1_4_mat[ds14,1], ] ==
    noNA_del_var_list[[4]][same_1_4_mat[ds14,2], ])
  if(tmp_same_genos != 10){
    del_same_1_4_diff_genos <- c(del_same_1_4_diff_genos, ds14)
  }
#print(tmp_same_genos)
}
del_same_1_4_diff_genos
# NULL - all genotypes are the same between the DELS that aren't 
#  identical but vary similar (and prob the same) between the two runs
####

del_int_2_3 <- intersect(rownames(noNA_del_var_list[[2]]),
  rownames(noNA_del_var_list[[3]]))
length(del_int_2_3)
# 17
del_int_2_3_diff_genos <- c()
for(di23 in del_int_2_3){
  tmp_same_genos <- sum(noNA_del_var_list[[2]][di23,] ==
    noNA_del_var_list[[3]][di23, ])
  if(tmp_same_genos != 10){
    del_int_2_3_diff_genos <- c(del_int_2_3_diff_genos, di23)
  }
  #print(tmp_same_genos)
}
del_int_2_3_diff_genos
#NULL - all the same-named DELs have the same genotypes
setdiff(rownames(noNA_del_var_list[[2]]), del_int_2_3)
setdiff(rownames(noNA_del_var_list[[3]]), del_int_2_3)
# 10 non-identical Deletions are likely same deletion:
# List2				# List3
same_2_3_vec <- c('Chr02_15912927_DEL_55',	'Chr02_15912927_DEL_57',
'Chr02_19900665_DEL_12317',	'Chr02_19900668_DEL_12316',
'Chr06_24972400_DEL_479',	'Chr06_24972400_DEL_478',
'Chr07_12613882_DEL_53',	'Chr07_12613882_DEL_54',
'Chr08_17464456_DEL_660',	'Chr08_17464456_DEL_659',
'Chr10_1245178_DEL_199',	'Chr10_1245178_DEL_198',
'Chr10_16534461_DEL_469',	'Chr10_16534461_DEL_471',
'Chr14_8258327_DEL_142',	'Chr14_8258322_DEL_141',
'Chr17_3368868_DEL_8785',	'Chr17_3368868_DEL_8786',
'Chr17_7812838_DEL_6635',	'Chr17_7812838_DEL_6636')

same_2_3_mat <- matrix(data = same_2_3_vec, ncol = 2, byrow = T)

del_same_2_3_diff_genos <- c()
for(ds23 in seq(nrow(same_2_3_mat))){
  tmp_same_genos <- sum(noNA_del_var_list[[2]][same_2_3_mat[ds23,1], ] ==
    noNA_del_var_list[[3]][same_2_3_mat[ds23,2], ])
  if(tmp_same_genos != 10){
    del_same_2_3_diff_genos <- c(del_same_2_3_diff_genos, ds23)
  }
#print(tmp_same_genos)
}
del_same_2_3_diff_genos
# NULL - all genotypes are the same between the DELS that aren't 
#  identical but vary similar (and prob the same) between the two runs
####

del_int_2_4 <- intersect(rownames(noNA_del_var_list[[2]]),
  rownames(noNA_del_var_list[[4]]))
length(del_int_2_4)
# 19
del_int_2_4_diff_genos <- c()
for(di24 in del_int_2_4){
  tmp_same_genos <- sum(noNA_del_var_list[[2]][di24,] ==
    noNA_del_var_list[[4]][di24, ])
  if(tmp_same_genos != 10){
    del_int_2_4_diff_genos <- c(del_int_2_4_diff_genos, di24)
  }
  #print(tmp_same_genos)
}
del_int_2_4_diff_genos
#NULL - all the same-named DELs have the same genotypes
setdiff(rownames(noNA_del_var_list[[2]]), del_int_2_4)
setdiff(rownames(noNA_del_var_list[[4]]), del_int_2_4)
# 9 non-identical Deletions are likely the same:
# List2				# List4
same_2_4_vec <- c('Chr02_19900665_DEL_12317',	'Chr02_19900668_DEL_12316',
'Chr06_9510637_DEL_294',	'Chr06_9510636_DEL_295',
'Chr06_24972400_DEL_479',	'Chr06_24972400_DEL_477',
'Chr07_12613882_DEL_53',	'Chr07_12613882_DEL_54',
'Chr08_17464456_DEL_660',	'Chr08_17464456_DEL_659',
'Chr10_1245178_DEL_199',	'Chr10_1245178_DEL_198',
'Chr10_16534461_DEL_469',	'Chr10_16534461_DEL_471',
'Chr17_3368868_DEL_8785',	'Chr17_3368868_DEL_8784',
'Chr17_7812838_DEL_6635',	'Chr17_7812838_DEL_6636')

same_2_4_mat <- matrix(data = same_2_4_vec, ncol = 2, byrow = T)

del_same_2_4_diff_genos <- c()
for(ds24 in seq(nrow(same_2_4_mat))){
  tmp_same_genos <- sum(noNA_del_var_list[[2]][same_2_4_mat[ds24,1], ] ==
    noNA_del_var_list[[4]][same_2_4_mat[ds24,2], ])
  if(tmp_same_genos != 10){
    del_same_2_4_diff_genos <- c(del_same_2_4_diff_genos, ds24)
  }
#print(tmp_same_genos)
}
del_same_2_4_diff_genos
# NULL - all genotypes are the same between the DELS that aren't 
#  identical but vary similar (and prob the same) between the two runs
####

del_int_3_4 <- intersect(rownames(noNA_del_var_list[[3]]),
  rownames(noNA_del_var_list[[4]]))
length(del_int_3_4)
# 24
del_int_3_4_diff_genos <- c()
for(di34 in del_int_3_4){
  tmp_same_genos <- sum(noNA_del_var_list[[3]][di34,] ==
    noNA_del_var_list[[4]][di34, ])
  if(tmp_same_genos != 10){
    del_int_3_4_diff_genos <- c(del_int_3_4_diff_genos, di34)
  }
  #print(tmp_same_genos)
}
del_int_3_4_diff_genos
# [1] "Chr13_5343435_DEL_20" - different genotypes at 2 samples, should
#  be removed from final list
setdiff(rownames(noNA_del_var_list[[3]]), del_int_3_4)
setdiff(rownames(noNA_del_var_list[[4]]), del_int_3_4)
# 4 non-identical Deletions are likelye the same:
# List3				# List4
same_3_4_vec <- c('Chr06_9510637_DEL_294',	'Chr06_9510636_DEL_295',
'Chr06_24972400_DEL_478',	'Chr06_24972400_DEL_477',
'Chr14_8258322_DEL_141',	'Chr14_8258327_DEL_142',
'Chr17_3368868_DEL_8786',	'Chr17_3368868_DEL_8784')

same_3_4_mat <- matrix(data = same_3_4_vec, ncol = 2, byrow = T)

del_same_3_4_diff_genos <- c()
for(ds34 in seq(nrow(same_3_4_mat))){
  tmp_same_genos <- sum(noNA_del_var_list[[3]][same_3_4_mat[ds34,1], ] ==
    noNA_del_var_list[[4]][same_3_4_mat[ds34,2], ])
  if(tmp_same_genos != 10){
    del_same_3_4_diff_genos <- c(del_same_3_4_diff_genos, ds34)
  }
#print(tmp_same_genos)
}
del_same_3_4_diff_genos
# NULL - all genotypes are the same between the DELS that aren't 
#  identical but vary similar (and prob the same) between the two runs
####

############
# OUTLINE FOR GENERATING SAME-NAME and DIFFERENT-NAME LISTS WITH FUNCTIONS:
# - For same-name, just need that they intersect. Would be good to check
#     that the genotypes are the same between both runs - may need separate
#     function for that
# - For different names: 
#   Break up the name string to extract chromosome, position, and size, and
#     use that info to join SVs using:
#   A) Same chromosome
#   B) Similar position (within 500bp (more/less?)
#   C) Similar size:
#     - <50bp - 25%
#     - 50bp < SV <200 - 10%
#     - >200bp - 5%
############

# Generat unified list:
# 1) Start with overlap between 1 and 2
# 2) Add same-but-nonidentical-name DELs from 1v2
# 3) Generate vector of "alternate names"
# 4) For next set, look for rownames that are NOT part of the new geno matrix
#      or the "alternate names" vector
# 5) Add new genotypes to geno matrix, add to "alternate names" geno matrix
#      as need be
# 6) At end, remove any "bad" DELs with genetoypes that differ between
#      runs

##########
# Functions for combining the results fromt the different runs
get_new_int_genos <- function(geno_list, int_name_vec, 
  uni_genos = uni_del_genos, alt_names = alt_name_vec){
  # Get SV name and genotype for SVs in the same-name list of two runs and
  #  that are NOT already in the united genotype matrix 
  # INPUTS #
  # geno_list = one of the two genotype matrices that were used to make
  #               int_name_vec; ex: noNA_del_var_list[[1]]
  # int_name_vec = vector of SVs that were present, with exact same names,
  #                  in two runs
  # uni_genos = matrix of the united genotypes; ex: uni_del_genos (default)
  # alt_names = vector of alternate names for same SV; 
  #               ex: alt_name_vec (default)
  # OUTPUT #
  # list with first element the names of the new SVs and second element
  #   a vector or matrix of the genotypes for the new SVs
  ###########
  new_int_svs <- setdiff(int_name_vec, c(rownames(uni_genos), alt_names))
  new_int_genos <- geno_list[new_int_svs, ]
  new_int_list <- list()
  new_int_list[['names']] <- new_int_svs
  new_int_list[['genos']] <- new_int_genos
  return(new_int_list)
}

# test_sv <- get_new_int_genos(noNA_del_var_list[[1]], del_int_1_3)

gen_new_uni_genos <- function(new_geno_list, uni_genos = uni_del_genos, 
  alt_names = alt_name_vec){
  # Generate a new unified genotype matrix that incorporates the new
  #   genotypes
  # INPUTS #
  # new_geno_list = list of [[1]] new SV names and [[2]] new SV genotypes
  # uni_genos = starting unified genotype matrix
  # alt_names = alternate SV name vector
  # OUTPUT #
  # Genotype matrix including the original unified geno matrix and the new
  #  SVs from new_geno_list
  #####
  tmp_sv_names <- rownames(uni_genos)
  tmp_full_sv_names <- c(tmp_sv_names, new_geno_list[['names']])
#  print(tmp_full_sv_names)
  tmp_genos <- rbind(uni_genos, new_geno_list[['genos']])
  rownames(tmp_genos) <- tmp_full_sv_names
  return(tmp_genos)
}

find_new_diffname_svs <- function(diff_name_mat, uni_genos = uni_del_genos, 
  alt_names = alt_name_vec){
  # Get names of SVs in the different-name marix that is not already in 
  #   the unified genotype matrix
  # INPUTS #
  # diff_name_mat = matrix of names of the same SV that have different name
  #                   in different runs. 2 columns, each column has the names
  #                   of the SV in the different run
  # uni_genos = matrix of the united genotypes
  # alt_names = vector of alternate names
  # OUTPUT #
  # vector of the new SV names, using the names in the first column of 
  #  diff_name_mat
  ###############
  tmp_in_genos <- apply(diff_name_mat, 1, 
    function(x) sum(x %in% c(rownames(uni_genos), alt_names)))
  new_inds <- which(tmp_in_genos == 0)
  new_svnames <- diff_name_mat[new_inds, 1]  
  return(new_svnames)
}

get_new_diffname_genos <- function(geno_list, diff_name_mat, 
  uni_genos = uni_del_genos, alt_names = alt_name_vec){
  # Gets the genotypes for the new SVs that have different names that were 
  #  identified using the <find_new_diffname_svs> function
  # INPUTS #
  # geno_list = the genotype matrix of one of the runs. Should be the matrix
  #               that corresponds to the first column of diff_name_mat
  # diff_name_mat = matrix of names of the same SV that have different name
  #                   in different runs. The 1st column should contain the 
  #                   names from <geno_list>
  # uni_genos = unified genotype matrix
  # alt_names = vector of alternate names of the SVs
  # OUTPUT #
  # List. [[1]] is the name(s) of the new SV(s). [[2]] is the genotypes of the
  #  new SV(s)
  #######
  new_svs <- find_new_diffname_svs(diff_name_mat)
  new_genos <- geno_list[new_svs, ]
  new_geno_list <- list()
  new_geno_list[['names']] <- new_svs
  new_geno_list[['genos']] <- new_genos
  return(new_geno_list)
}
#################

# All the SVs shared between 2+ runs but DON'T have the same genotypes in
#   each run
bad_svs <- unique(c(del_int_1_2_diff_genos, del_same_1_2_diff_genos,
del_int_1_3_diff_genos, del_same_1_3_diff_genos,
del_int_1_4_diff_genos, del_same_1_4_diff_genos,
del_int_2_3_diff_genos, del_same_2_3_diff_genos,
del_int_2_4_diff_genos, del_same_2_4_diff_genos,
del_int_3_4_diff_genos, del_same_3_4_diff_genos))

# List containing all the vectors containing the names of the SVs that are the
#  same in pairwise comparisons of runs
int_vec_list <- list()
int_vec_list[[1]] <- del_int_1_2
int_vec_list[[2]] <- del_int_1_3
int_vec_list[[3]] <- del_int_1_4
int_vec_list[[4]] <- del_int_2_3
int_vec_list[[5]] <- del_int_2_4
int_vec_list[[6]] <- del_int_3_4

# shows which run is the first run in each comparison for int_vec_list and 
#   diff_mat_list
run_ind_vec <- c(1,1,1,2,2,3)

# List containing the matrices of the names of SVs that are the same but have
#   different names in the different pairwise comparisons
diff_mat_list <- list()
diff_mat_list[[1]] <- same_1_2_mat
diff_mat_list[[2]] <- same_1_3_mat
diff_mat_list[[3]] <- same_1_4_mat
diff_mat_list[[4]] <- same_2_3_mat
diff_mat_list[[5]] <- same_2_4_mat
diff_mat_list[[6]] <- same_3_4_mat

# generate united genotype matrix from same- and different-named SVs from
#  runs 1 and 2
uni_del_genos <- noNA_del_var_list[[1]][del_int_1_2, ]
tmp_same_genos <- noNA_del_var_list[[1]][same_1_2_mat[,1], ]
uni_del_genos <- rbind(uni_del_genos, tmp_same_genos)
# generate alternate name vector
alt_name_vec <- setdiff(as.vector(same_1_2_mat), rownames(uni_del_genos))

setdiff(del_int_1_3, c(rownames(uni_del_genos), alt_name_vec))

# Add SVs with the same name in 2 runs and aren't included in the first 2 runs
for(j in c(2:length(int_vec_list))){
  tmp_run_mat <- noNA_del_var_list[[run_ind_vec[j]]]
  comp_names <- int_vec_list[[j]]
  tmp_sv_geno_list <- get_new_int_genos(tmp_run_mat, comp_names, 
    uni_genos = uni_del_genos)
  uni_del_genos <- gen_new_uni_genos(tmp_sv_geno_list, 
    uni_genos = uni_del_genos)
}

# Add SVs with different names in 2 runs and aren't included in 1st 2 runs
for(k in c(2:length(diff_mat_list))){
  tmp_run_mat <- noNA_del_var_list[[run_ind_vec[k]]]
  tmp_name_mat <- diff_mat_list[[k]]
  tmp_diffname_genos <- get_new_diffname_genos(geno_list = tmp_run_mat, 
    diff_name_mat = tmp_name_mat, uni_genos = uni_del_genos, 
    alt_names = alt_name_vec)
  uni_del_genos <- gen_new_uni_genos(tmp_diffname_genos, 
    uni_genos = uni_del_genos, alt_names = alt_name_vec)
  tmp_new_altnames <- setdiff(as.vector(tmp_name_mat), 
    c(rownames(uni_del_genos), alt_name_vec))
  alt_name_vec <- c(alt_name_vec, tmp_new_altnames)
}

good_inds <- setdiff(seq(nrow(uni_del_genos)), 
  which(rownames(uni_del_genos) %in% bad_svs))
uni_del_genos <- uni_del_genos[good_inds, ]

uni_2_inds <- which(apply(uni_del_genos, 1, function(x) sum(x == 2) > 0))

uni_del_genos_no2s <- uni_del_genos[-uni_2_inds, ]

table(apply(uni_del_genos, 1, function(x) min(table(x))))
#  1  2  4 
# 32  3  1 

table(apply(uni_del_genos_no2s, 1, function(x) min(table(x))))
# 1 
# 20

save.image(file = paste(analysis_res_dir, 'pb_DEL_analysis_workspace.RData', 
  sep = ''))
# /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/analysis_results/pb_DEL_analysis_workspace.RData
# load('/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/analysis_results/pb_DEL_analysis_workspace.RData')

uni_del_geno_file <- paste(analysis_res_dir, 'pbsv_4run_unified_DELs.txt',
  sep = '')
write.table(uni_del_genos, file = uni_del_geno_file, quote = F, sep = '\t',
  row.names = T, col.names = T)

# generate multifasta file for BLAST querries

which((rownames(uni_del_genos) %in% precall_filt_list[[1]]$full_name_long) 
  == F)
# [1] 31 32 33 34
which((rownames(uni_del_genos) %in% precall_filt_list[[3]]$full_name_long) 
  == F)
# [1]  4  5 17 18 21 23 25 27 35 36

set_1_inds <- which(precall_filt_list[[1]]$full_name_long 
  %in% rownames(uni_del_genos))
name_and_seq <- precall_filt_list[[1]][set_1_inds, c('full_name_long', 'REF')]
set_3_inds <- which(precall_filt_list[[3]]$full_name_long 
  %in% rownames(uni_del_genos))
name_and_seq <- rbind(name_and_seq, 
  precall_filt_list[[3]][set_3_inds, c('full_name_long', 'REF')])

name_and_seq <- name_and_seq[-which(duplicated(name_and_seq$full_name_long)), 
  ]
name_and_seq$full_name_long <- paste('>', name_and_seq$full_name_long, sep = '')

name_and_seq_vec <- c(rbind(name_and_seq[,1], name_and_seq[,2]))

name_and_seq_file <- paste(analysis_res_dir, 'uni_del_seqs.fasta', sep = '')
write.table(name_and_seq_vec, file = name_and_seq_file, quote = F, sep = '\t',
  row.names = F, col.names = F)

# import BLAST results

uni_del_blast_res_file <- paste(analysis_res_dir, 'blast55728.fmt_6.out', 
  sep = '')

del_blast_res <- read.table(file = uni_del_blast_res_file, sep = '\t', 
  header = F, stringsAsFactors = F)

colnames(del_blast_res) <- c('query_id', 'subject_id', 'per_identity', 
  'alignment_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 
  's_start', 's_end', 'evalue', 'bit_score')

del_blast_list <- list()
for(dn in unique(del_blast_res$query_id)){
  tmp_inds <- which(del_blast_res$query_id == dn)
  del_blast_list[[dn]] <- del_blast_res[tmp_inds, ]
}

# have checked through index 18

get_hit_chroms <- function(blast_hit_df, per_score_cut = 0.9){
  tmp_df_sort <- blast_hit_df[order(blast_hit_df$bit_score, decreasing = T), ]
  ok_scores <- which(tmp_df_sort$bit_score > 
    (per_score_cut*tmp_df_sort$bit_score[1]))
  quer_chrom <- unlist(strsplit(tmp_df_sort$query_id[1], split = '_'))[1]
  hit_chr_split <- strsplit(tmp_df_sort$subject_id[ok_scores], split = ' ')
  hit_chrs <- unlist(lapply(hit_chr_split, function(x) x[1]))
  hit_chrs_list <- list()
  hit_chrs_list[['quer_chrom']] <- quer_chrom
  hit_chrs_list[['hit_chrs']] <- hit_chrs
  return(hit_chrs_list)
}

hit_chroms_list <- lapply(del_blast_list, get_hit_chroms, per_score_cut = 0.9)

n_diff_chroms <- lapply(hit_chroms_list, 
  function(x) sum(x$quer_chrom != x$hit_chrs))

diff_chrom_del_names <- names(n_diff_chroms)

n_diff_chrom_df <- data.frame(del_name = diff_chrom_del_names, 
  n_diff_chroms = unlist(n_diff_chroms))

# DELs that vary at presence of REF allele: only have REF/ALT and ALT/ALT
# Number of DELs that vary at REF allele with no hits on other chromosomes = 5
# Number of DELS that vary at REF with hits on 1 or more other chromoeomes = 11
#  1 = 4; 2 = 3; 3 = 1; 21 = 1; 28 = 1; 36 = 1

# DELs that vary at presence of ALT allele: only have REF/REF and REF/ALT
# Number of DELS that vary at ALT allele with no hits on other chromosomes = 10
# Number of DELs that vary at ALT with hits on 1 + other chromosomes = 8
# 1 = 1; 2 = 4; 12 = 2; 43 = 1; NA = 2

quit(save = 'no')

