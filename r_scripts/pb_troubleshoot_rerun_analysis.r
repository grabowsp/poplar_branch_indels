# Script for analysis of VCFs from re-running the PacBio pipeline
#   to troubleshoot the sample-order heterozygosity bias I saw in
#   the results of the first run

function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(ape)

# LOAD DATA #
data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'

analysis_res_dir <- paste(data_dir, 'analysis_results/', sep = '')


combo_file_name_vec <- c('ref.ALLData.vcf', 'ref.ALLData_try2.vcf', 'ref.ALLData.try3.vcf', 'ref.ALLData.try4.vcf')

samp_order_list <- list()
samp_order_list[[1]] <- c("PAXL", "PAXN", "PAYK", "PAYZ", "PAZF", "PAZG", 
  "PAZH", "PBAT", "PBAU", "PBAW")
samp_order_list[[2]] <- c('PBAW', 'PBAU', 'PBAT', 'PAZH', 'PAZG', 'PAZF', 
  'PAYZ', 'PAYK', 'PAXN', 'PAXL')
samp_order_list[[3]] <- c('PBAU', 'PAYZ', 'PBAW', 'PAXN', 'PBAT', 'PAXL', 
  'PAZG', 'PAYK', 'PAZF', 'PAZH')
samp_order_list[[4]] <- c('PAYZ', 'PAZF', 'PAXN', 'PAZG', 'PAXL', 'PBAT', 
  'PAYK', 'PBAW', 'PAZH', 'PBAU')

geno_tab_list <- list()
for(cfnv in seq(length(combo_file_name_vec))){
  tmp_tot_name <- paste(data_dir, combo_file_name_vec[cfnv], sep = '')
  tmp_df <- make_combo_indel_df(tmp_tot_name)
  tmp_geno_info <- make_allsamp_geno_info_list(combo_df = tmp_df)
  tmp_lib_names <- names(tmp_geno_info)
  tmp_lib_ord <- samp_order_list[[cfnv]]
  tmp_samp_ord <- match(tmp_lib_ord, tmp_lib_names)
  tmp_geno_df <- data.frame(lib = tmp_lib_names, samp_order = tmp_samp_ord,
    homRef = NA, het = NA, homAlt = NA, stringsAsFactors = F)
  for(i in names(tmp_geno_info)){
    tab_ind <- which(tmp_geno_df$lib == i)
    tmp_geno_tab <- table(tmp_geno_info[[i]][1])
    tmp_geno_df$homRef[tab_ind] <- tmp_geno_tab[1]
    tmp_geno_df$het[tab_ind] <- tmp_geno_tab[2]
    tmp_geno_df$homAlt[tab_ind] <- tmp_geno_tab[3]
    tmp_geno_df$per_het <- (tmp_geno_df$het / 
      apply(tmp_geno_df[,c('homRef', 'het', 'homAlt')], 1, sum))
    tmp_geno_df$het_ord <- order(tmp_geno_df$per_het, decreasing = T)
  }
  colnames(tmp_geno_df) <- paste(colnames(tmp_geno_df), cfnv, sep = '_')
  geno_tab_list[[cfnv]] <- tmp_geno_df
  print(cfnv)
}

tot_geno_df <- cbind(cbind(geno_tab_list[[1]], geno_tab_list[[2]]), 
  cbind(geno_tab_list[[3]], geno_tab_list[[4]]))

geno_tab_file_out <- paste(analysis_res_dir, 
  'reruns_combo_raw_genotype_table.txt', sep = '')

write.table(tot_geno_df, file = geno_tab_file_out, quote = F, sep = '\t', 
  row.names = F, col.names = T)


##########

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

var_geno_list <- list()
for(j in seq(length(noNA_geno_list))){
  tmp_genos <- noNA_geno_list[[j]]
  tmp_tab <- apply(tmp_genos, 1, table)
  nonvar_inds <- which(unlist(lapply(tmp_tab, length)) == 1)
  var_genos <- tmp_genos[-nonvar_inds, ]
  var_geno_list[[j]] <- var_genos
}

var_geno_tally_list <- list()
for(vgl in seq(length(var_geno_list))){
  tmp_genos <- var_geno_list[[vgl]]
  tmp_lib_names <- colnames(tmp_genos)
  tmp_lib_ord <- samp_order_list[[vgl]]
  tmp_samp_ord <- match(tmp_lib_ord, tmp_lib_names)
  tmp_geno_df <- data.frame(lib = tmp_lib_names, samp_order = tmp_samp_ord,
    homRef = NA, het = NA, homAlt = NA, stringsAsFactors = F)
  tmp_geno_df$homRef <- apply(tmp_genos, 2, function(x) sum(x == 0))
  tmp_geno_df$het <- apply(tmp_genos, 2, function(x) sum(x == 1))
  tmp_geno_df$homAlt <- apply(tmp_genos, 2, function(x) sum(x == 2))
  tmp_geno_df$het_ord <- order(tmp_geno_df$het, decreasing = T)
  colnames(tmp_geno_df) <- paste(colnames(tmp_geno_df), vgl, sep = '_')
  var_geno_tally_list[[vgl]] <- tmp_geno_df
  print(vgl)
}

tot_filt_df <- cbind(cbind(var_geno_tally_list[[1]], var_geno_tally_list[[2]]),
  cbind(var_geno_tally_list[[3]], var_geno_tally_list[[4]]))

filt_tab_file_out <- paste(analysis_res_dir,
  'reruns_combo_filtered_genotype_table.txt', sep = '')

write.table(tot_filt_df, file = filt_tab_file_out, quote = F, sep = '\t',
  row.names = F, col.names = T)


########
# Continue from here
# Next would be generating trees



combo_file_2 <- 'ref.ALLData_try2.vcf'
combo_file_2_tot <- paste(data_dir, combo_file_2, sep = '')
combo_df_2 <- make_combo_indel_df(combo_file_2_tot)

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')
# samp_meta[,c('lib_name', 'branch_name')]

geno_info_2 <- make_allsamp_geno_info_list(combo_df = combo_df_2)
# make basic dataframe to tally genotypes

samp_order_2 <- c('PBAW', 'PBAU', 'PBAT', 'PAZH', 'PAZG', 'PAZF', 'PAYZ', 
  'PAYK', 'PAXN', 'PAXL')

basic_samp_tab_df <- data.frame(lib =names(geno_info_2), samp_order = c(1:10),
  homRef = NA, het = NA, homAlt = NA,
  stringsAsFactors = F)

for(i in names(geno_info)){
  tab_ind <- which(basic_samp_tab_df$lib == i)
  tmp_geno_tab <- table(geno_info[[i]][1])
  basic_samp_tab_df$homRef[tab_ind] <- tmp_geno_tab[1]
  basic_samp_tab_df$het[tab_ind] <- tmp_geno_tab[2]
  basic_samp_tab_df$homAlt[tab_ind] <- tmp_geno_tab[3]
}

basic_samp_tab_df$per_het <- (basic_samp_tab_df$het / 
  apply(basic_samp_tab_df[,c('homRef', 'het', 'homAlt')], 1, sum))

# NEXT: 
# What I want to do:
# For each run, get sample order, percent heterozygosity, and per_het 
#  rank/order




# Change order of samples and make name vector for SV_df

branch_name_ord <- c(13.1, 13.2, 13.3, 13.4, 13.5, 14.5, 14.1, 14.4, 14.3, 14.2)
lib_name_ord <- c()
for(bn in branch_name_ord){
  tmp_ind <- which(samp_meta$branch_name == bn)
  lib_name_ord <- c(lib_name_ord, samp_meta$lib_name[tmp_ind])
}
branch_name_vec <- paste('branch', branch_name_ord, sep = '_')

first_samp_col <- which(colnames(combo_df_filt) == 'FORMAT') + 1
last_samp_col <- which(colnames(combo_df_filt) == 'full_name') - 1

combo_df_filt[, c(first_samp_col:last_samp_col)] <- combo_df_filt[,lib_name_ord]
colnames(combo_df_filt)[c(first_samp_col:last_samp_col)] <- lib_name_ord

# make filtered metadata

fs_inds <- c()
for(bn in branch_name_ord){
  tmp_m_ind <- which(samp_meta$branch_name == bn)
  fs_inds <- c(fs_inds, tmp_m_ind)
}

samp_meta_filt <- samp_meta[fs_inds, ]



basic_singleton_df <- data.frame(lib = lib_name_ord,
  branch = factor(branch_name_vec, levels = branch_name_vec),
  tree = c(rep('13', times = 5), rep('14', times = 5)),
  homRef = 0, het = 0, homAlt = 0, stringsAsFactors = F)

# required for functions tallying tree-specific SVs
branch_13_lab <- c(13.1, 13.2, 13.3, 14.5)
branch_14_lab <- c(14.2, 14.3, 14.4, 13.5)

# Things to show for SVs:
#  1) trees of samples
#  2) Total number of SVs
#  3) Total number of variant SVs
#  4) Number of homozygous/heterozygous SVs in each sample
#  4) Numbers of tree-specific SVs
#  5) Numbers of branch-specific SVs
#  6) Discordant SVs - present AND absent in both trees

# Filtering criteria to tinker with: sv_size, dist_cut, 
#  het_ratio_cut, max_hom_ratio,
#  good_score, ss_min_cov

##########################################################
# First Genotypes: strict, 100bp min
precall_filt_df_1 <- precalling_SV_filtering(sv_geno_df = combo_df_filt,
  mer_length = 8, per_mn_cutoff = 0.7, per_bn_pure_cutoff = 0.5, 
  per_bn_multi_cutoff = 0.6, dist_cut = 1000, sd_cut = 3, use_sd = T)

combo_genotypes_1 <- generate_filtered_genotypes(combo_df = precall_filt_df_1,
  min_sv_coverage = 10, het_ratio_cut = 0.25, min_minor_allele_count = 2,
  max_hom_ratio = 0.05, est_error = 0.032, est_penetrance = 0.35,
  good_score = 0.9, great_score = 0.98, same_geno_bonus = 0.67,
  ss_min_cov = 20)

combo_genos_noNA_100 <- filt_geno_mat(geno_mat = combo_genotypes_1,
  max_nas = 0, min_length = 100)

library(ape)
# generate distance matrix and trees
geno1_t <- t(combo_genos_noNA_100)

rownames(geno1_t) <- branch_name_vec
geno1_dist_mat <- dist(geno1_t, method = 'manhattan', diag = T,
  upper = T)

geno1_nj <- nj(geno1_dist_mat)
geno1_upgma <- hclust(geno1_dist_mat, method = 'average')

geno1_tree_file <- paste(data_dir,
  'analysis_results/figs/filt1_sv_dist_trees.png', sep = '')

png(filename = geno1_tree_file, width = 1000, height = 500)
par(mfrow = c(1,2))
plot(geno1_nj, main = 'Poplar SV Branch NJ Tree\nFiltered SV set 1')
plot(geno1_upgma, main = 'Poplar SV Branch UPGMA Tree\nFiltered SV set 1')
dev.off()

# Total number of SVs
geno1_totSV <- nrow(combo_genos_noNA_100)
# 8392

n_geno_vec <- apply(combo_genos_noNA_100, 1, function(x) length(table(x)))

geno1_var <- combo_genos_noNA_100[which(n_geno_vec > 1), ]
geno1_varSV <- nrow(geno1_var)
# 59

# Number of each genotype in each sample
geno1_samp_tab_df <- gen_samp_geno_tab_df(geno_mat = geno1_var, 
  basic_df = basic_samp_tab_df)

geno1_barplot_info <- gen_samp_geno_barplot_ggs(
  geno_tab_df = geno1_samp_tab_df, 
  data_lab = 'Filtered SV Set 1')

geno1_combo_bar_file <- paste(data_dir,
  'analysis_results/figs/filt1_sv_genotype_bars.png', sep = '')

png(geno1_combo_bar_file, width = 1200, height = 400)
do.call(grid.arrange, c(geno1_barplot_info, ncol = 3))
dev.off()

# Number of unique/singleton genotypes in each sample

geno1_singletons <- gen_singleton_df(geno_mat = geno1_var, 
  basic_df = basic_singleton_df)

geno1_combo_singleton_bar_file <- paste(data_dir,
  'analysis_results/figs/filt1_sv_singleton_bars.png', sep = '')

geno1_singleton_bar_info <- gen_samp_geno_barplot_ggs(
  geno_tab_df = geno1_singletons,
  data_lab = 'Filtered SV Set 1 Singleton Genotypes')

png(geno1_combo_singleton_bar_file, width = 1200, height = 400)
do.call(grid.arrange, c(geno1_singleton_bar_info, ncol = 3))
dev.off()

# Calculate number of tree-specific SVs
# Start with clone 14 specific SVs
geno1_14_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno1_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_14_lab, meta = samp_meta)

length(geno1_14_spec_hets)
# 1

geno1_13_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno1_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_13_lab, meta = samp_meta)

length(geno1_13_spec_hets)
# 1

#########################################################
# Strict filter, shorter length cutoff
combo_genos_noNA_50 <- filt_geno_mat(geno_mat = combo_genotypes_1,
  max_nas = 0, min_length = 50)

library(ape)
# generate distance matrix and trees
geno2_t <- t(combo_genos_noNA_50)

rownames(geno2_t) <- branch_name_vec
geno2_dist_mat <- dist(geno2_t, method = 'manhattan', diag = T,
  upper = T)

geno2_nj <- nj(geno2_dist_mat)
geno2_upgma <- hclust(geno2_dist_mat, method = 'average')

geno2_tree_file <- paste(data_dir,
  'analysis_results/figs/filt2_sv_dist_trees.png', sep = '')

png(filename = geno2_tree_file, width = 1000, height = 500)
par(mfrow = c(1,2))
plot(geno2_nj, main = 'Poplar SV Branch NJ Tree\nFiltered SV set 2')
plot(geno2_upgma, main = 'Poplar SV Branch UPGMA Tree\nFiltered SV set 2')
dev.off()

# Total number of SVs
geno2_totSV <- nrow(combo_genos_noNA_50)
# 12690

geno2_n_geno_vec <- apply(combo_genos_noNA_50, 1, function(x) length(table(x)))

geno2_var <- combo_genos_noNA_50[which(geno2_n_geno_vec > 1), ]
geno2_varSV <- nrow(geno2_var)
# 239

# Number of each genotype in each sample
geno2_samp_tab_df <- gen_samp_geno_tab_df(geno_mat = geno2_var,
  basic_df = basic_samp_tab_df)

geno2_barplot_info <- gen_samp_geno_barplot_ggs(
  geno_tab_df = geno2_samp_tab_df,
  data_lab = 'Filtered SV Set 2')

geno2_combo_bar_file <- paste(data_dir,
  'analysis_results/figs/filt2_sv_genotype_bars.png', sep = '')

png(geno2_combo_bar_file, width = 1200, height = 400)
do.call(grid.arrange, c(geno2_barplot_info, ncol = 3))
dev.off()

geno2_singletons <- gen_singleton_df(geno_mat = geno2_var,
  basic_df = basic_singleton_df)

geno2_combo_singleton_bar_file <- paste(data_dir,
  'analysis_results/figs/filt2_sv_singleton_bars.png', sep = '')

geno2_singleton_bar_info <- gen_samp_geno_barplot_ggs(
  geno_tab_df = geno2_singletons,
  data_lab = 'Filtered SV Set 2 Singleton Genotypes')

png(geno2_combo_singleton_bar_file, width = 1200, height = 400)
do.call(grid.arrange, c(geno2_singleton_bar_info, ncol = 3))
dev.off()

# Calculate number of tree-specific SVs
# Start with clone 14 specific SVs
geno2_14_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno2_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_14_lab, meta = samp_meta)

length(geno2_14_spec_hets)
# 33

geno2_13_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno2_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_13_lab, meta = samp_meta)

length(geno2_13_spec_hets)
# 3

##########################
# Finally, try cutting off at SV length of 32
# Strict filter, shorter length cutoff
combo_genos_noNA_32 <- filt_geno_mat(geno_mat = combo_genotypes_1,
  max_nas = 0, min_length = 32)

library(ape)
# generate distance matrix and trees
geno3_t <- t(combo_genos_noNA_32)

rownames(geno3_t) <- branch_name_vec
geno3_dist_mat <- dist(geno3_t, method = 'manhattan', diag = T,
  upper = T)

geno3_nj <- nj(geno3_dist_mat)
geno3_upgma <- hclust(geno3_dist_mat, method = 'average')

geno3_tree_file <- paste(data_dir,
  'analysis_results/figs/filt3_sv_dist_trees.png', sep = '')

png(filename = geno3_tree_file, width = 1000, height = 500)
par(mfrow = c(1,2))
plot(geno3_nj, main = 'Poplar SV Branch NJ Tree\nFiltered SV set 3')
plot(geno3_upgma, main = 'Poplar SV Branch UPGMA Tree\nFiltered SV set 3')
dev.off()

# Total number of SVs
geno3_totSV <- nrow(combo_genos_noNA_32)
# 15538

geno3_n_geno_vec <- apply(combo_genos_noNA_32, 1, function(x) length(table(x)))

geno3_var <- combo_genos_noNA_32[which(geno3_n_geno_vec > 1), ]
geno3_varSV <- nrow(geno3_var)
# 821

# Number of each genotype in each sample
geno3_samp_tab_df <- gen_samp_geno_tab_df(geno_mat = geno3_var,
  basic_df = basic_samp_tab_df)

geno3_barplot_info <- gen_samp_geno_barplot_ggs(
  geno_tab_df = geno3_samp_tab_df,
  data_lab = 'Filtered SV Set 3')

geno3_combo_bar_file <- paste(data_dir,
  'analysis_results/figs/filt3_sv_genotype_bars.png', sep = '')

png(geno3_combo_bar_file, width = 1200, height = 400)
do.call(grid.arrange, c(geno3_barplot_info, ncol = 3))
dev.off()

geno3_singletons <- gen_singleton_df(geno_mat = geno3_var,
  basic_df = basic_singleton_df)

geno3_combo_singleton_bar_file <- paste(data_dir,
  'analysis_results/figs/filt3_sv_singleton_bars.png', sep = '')

geno3_singleton_bar_info <- gen_samp_geno_barplot_ggs(
  geno_tab_df = geno3_singletons,
  data_lab = 'Filtered SV Set 3 Singleton Genotypes')

png(geno3_combo_singleton_bar_file, width = 1200, height = 400)
do.call(grid.arrange, c(geno3_singleton_bar_info, ncol = 3))
dev.off()

# Calculate number of tree-specific SVs
# Start with clone 14 specific SVs
geno3_14_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno3_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_14_lab, meta = samp_meta)

length(geno3_14_spec_hets)
# 114

geno3_13_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno3_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_13_lab, meta = samp_meta)

length(geno3_13_spec_hets)
# 4

#################################
# Try not using minimum coverage for at least one sample

combo_genotypes_2 <- generate_filtered_genotypes(combo_df = precall_filt_df_1,
  min_sv_coverage = 10, het_ratio_cut = 0.25, min_minor_allele_count = 2,
  max_hom_ratio = 0.05, est_error = 0.032, est_penetrance = 0.35,
  good_score = 0.9, great_score = 0.98, same_geno_bonus = 0.67,
  ss_min_cov = 10)

combo_genos_4_32 <- filt_geno_mat(geno_mat = combo_genotypes_2,
  max_nas = 0, min_length = 32)

# generate distance matrix and trees
geno4_t <- t(combo_genos_4_32)

rownames(geno4_t) <- branch_name_vec
geno4_dist_mat <- dist(geno4_t, method = 'manhattan', diag = T,
  upper = T)

geno4_nj <- nj(geno4_dist_mat)
geno4_upgma <- hclust(geno4_dist_mat, method = 'average')

geno4_tree_file <- paste(data_dir,
  'analysis_results/figs/filt4_sv_dist_trees.png', sep = '')

png(filename = geno4_tree_file, width = 1000, height = 500)
par(mfrow = c(1,2))
plot(geno4_nj, main = 'Poplar SV Branch NJ Tree\nFiltered SV set 4')
plot(geno4_upgma, main = 'Poplar SV Branch UPGMA Tree\nFiltered SV set 4')
dev.off()

# Total number of SVs
geno4_totSV <- nrow(combo_genos_4_32)
# 15551

geno4_n_geno_vec <- apply(combo_genos_4_32, 1, function(x) length(table(x)))

geno4_var <- combo_genos_4_32[which(geno4_n_geno_vec > 1), ]
geno4_varSV <- nrow(geno4_var)
# 826

# using the miminum single-sample coverage hardly affects the total number of 
#   SVs

########################################
# Try using more permissive het_ratio_cut and estimated penetrance

combo_genotypes_3 <- generate_filtered_genotypes(combo_df = precall_filt_df_1,
  min_sv_coverage = 10, het_ratio_cut = 0.15, min_minor_allele_count = 2,
  max_hom_ratio = 0.05, est_error = 0.032, est_penetrance = 0.3,
  good_score = 0.9, great_score = 0.98, same_geno_bonus = 0.67,
  ss_min_cov = 10)

combo_genos_5_32 <- filt_geno_mat(geno_mat = combo_genotypes_3,
  max_nas = 0, min_length = 32)

# generate distance matrix and trees
geno5_t <- t(combo_genos_5_32)

rownames(geno5_t) <- branch_name_vec
geno5_dist_mat <- dist(geno5_t, method = 'manhattan', diag = T,
  upper = T)

geno5_nj <- nj(geno5_dist_mat)
geno5_upgma <- hclust(geno5_dist_mat, method = 'average')

geno5_tree_file <- paste(data_dir,
  'analysis_results/figs/filt5_sv_dist_trees.png', sep = '')

png(filename = geno5_tree_file, width = 1000, height = 500)
par(mfrow = c(1,2))
plot(geno5_nj, main = 'Poplar SV Branch NJ Tree\nFiltered SV set 5')
plot(geno5_upgma, main = 'Poplar SV Branch UPGMA Tree\nFiltered SV set 5')
dev.off()

# Total number of SVs
geno5_totSV <- nrow(combo_genos_5_32)
# 21797

geno5_n_geno_vec <- apply(combo_genos_5_32, 1, function(x) length(table(x)))

geno5_var <- combo_genos_5_32[which(geno5_n_geno_vec > 1), ]
geno5_varSV <- nrow(geno5_var)
# 1330

# Number of each genotype in each sample
geno5_samp_tab_df <- gen_samp_geno_tab_df(geno_mat = geno5_var,
  basic_df = basic_samp_tab_df)

geno5_barplot_info <- gen_samp_geno_barplot_ggs(
  geno_tab_df = geno5_samp_tab_df,
  data_lab = 'Filtered SV Set 5')

geno5_combo_bar_file <- paste(data_dir,
  'analysis_results/figs/filt5_sv_genotype_bars.png', sep = '')

png(geno5_combo_bar_file, width = 1200, height = 400)
do.call(grid.arrange, c(geno5_barplot_info, ncol = 3))
dev.off()

geno5_singletons <- gen_singleton_df(geno_mat = geno5_var,
  basic_df = basic_singleton_df)

geno5_combo_singleton_bar_file <- paste(data_dir,
  'analysis_results/figs/filt5_sv_singleton_bars.png', sep = '')

geno5_singleton_bar_info <- gen_samp_geno_barplot_ggs(
  geno_tab_df = geno5_singletons,
  data_lab = 'Filtered SV Set 5 Singleton Genotypes')

png(geno5_combo_singleton_bar_file, width = 1200, height = 400)
do.call(grid.arrange, c(geno5_singleton_bar_info, ncol = 3))
dev.off()

# Calculate number of tree-specific SVs
# Start with clone 14 specific SVs
geno5_14_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno5_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_14_lab, meta = samp_meta)

length(geno5_14_spec_hets)
# 156

geno5_13_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno5_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_13_lab, meta = samp_meta)

length(geno5_13_spec_hets)
# 6


quit(save = 'n')







