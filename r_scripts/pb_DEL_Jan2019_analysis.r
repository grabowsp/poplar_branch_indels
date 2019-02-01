# Analysis looking at only deletions from the pbsv pipeline

# Outline
# Do all analyses for all 4 runs to show repeatibility
# Generate trees for all 10 samples to make sure they make sense
# If 10-samp trees make sense, then do rest on only 8 "good" samples
# Tally number of genotypes by sample
# Tally number of singleton DELs by sample
# Tally clone-specific DELs
# Tally 2-, 3-, and 4-sample DELs

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

noNA_del_list <- list()
for(ngl in seq(length(noNA_geno_list))){
  tmp_del_inds <- grep('DEL', rownames(noNA_geno_list[[ngl]]))
  tmp_del_mat <- noNA_geno_list[[ngl]][tmp_del_inds, ]
  noNA_del_list[[ngl]] <- tmp_del_mat
}

noNA_del_var_list <- list()
for(ndvl in seq(length(noNA_del_list))){
  tmp_n_al <- apply(noNA_del_list[[ndvl]], 1, function(x) length(unique(x)))
  tmp_var_inds <- which(tmp_n_al == 2)
  noNA_del_var_list[[ndvl]] <- noNA_del_list[[ndvl]][tmp_var_inds, ]
}

n_var_dels <- unlist(lapply(noNA_del_var_list, nrow))
# [1] 35 36 38 43

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


quit(save = 'no')

