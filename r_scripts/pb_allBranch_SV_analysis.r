# Script for analysis of filtered SVs generated by PacBio SV caller pipeline

function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(ape)

# LOAD DATA #
data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'
combo_file <- 'ref.ALLData.vcf'
combo_file_tot <- paste(data_dir, combo_file, sep = '')
combo_df <- make_combo_indel_df(combo_file_tot)

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')
# samp_meta[,c('lib_name', 'branch_name')]

## Remove bad branches
#bad_branches <- c('13.4', '14.1')
#bad_samps <- c()
#for(bb in bad_branches){
#  tmp_ind <- which(as.character(samp_meta$branch_name) == bb)
#  bad_samps <- c(bad_samps, samp_meta$lib_name[tmp_ind])
#}

#bad_vcf_cols <- c()
#for(bs in bad_samps){
#  tmp_col_ind <- which(colnames(combo_df) == bs)
#  bad_vcf_cols <- c(bad_vcf_cols, tmp_col_ind)
#}

#combo_df_filt <- combo_df[ , -sort(bad_vcf_cols)]

combo_df_filt <- combo_df

# Change order of samples and make name vector for SV_df

branch_name_ord <- c(13.1, 13.2, 13.3, 14.5, 13.4, 14.2, 14.3, 14.4, 14.1, 13.5)
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

# make basic dataframe to tally genotypes
basic_samp_tab_df <- data.frame(lib = lib_name_ord,
  branch = factor(branch_name_vec, levels = branch_name_vec),
  tree = c(rep('13', times = 5), rep('14', times = 5)),
  homRef = NA, het = NA, homAlt = NA,
  stringsAsFactors = F)

basic_singleton_df <- data.frame(lib = lib_name_ord,
  branch = factor(branch_name_vec, levels = branch_name_vec),
  tree = c(rep('13', times = 5), rep('14', times = 5)),
  homRef = 0, het = 0, homAlt = 0, stringsAsFactors = F)

# required for functions tallying tree-specific SVs
branch_13_lab <- c(13.1, 13.2, 13.3, 14.5, 13.4)
branch_14_lab <- c(14.2, 14.3, 14.4, 14.1, 13.5)

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
# First Genotypes: strict, 32bp min
precall_filt_df_1 <- precalling_SV_filtering(sv_geno_df = combo_df_filt,
  mer_length = 8, per_mn_cutoff = 0.7, per_bn_pure_cutoff = 0.5, 
  per_bn_multi_cutoff = 0.6, dist_cut = 1000, sd_cut = 3, use_sd = T)

combo_genotypes_1 <- generate_filtered_genotypes(combo_df = precall_filt_df_1,
  min_sv_coverage = 10, het_ratio_cut = 0.25, min_minor_allele_count = 2,
  max_hom_ratio = 0.05, est_error = 0.032, est_penetrance = 0.35,
  good_score = 0.9, great_score = 0.98, same_geno_bonus = 0.67,
  ss_min_cov = 20)

combo_genos_noNA_32 <- filt_geno_mat(geno_mat = combo_genotypes_1,
  max_nas = 0, min_length = 32)

# generate distance matrix and trees
geno1_t <- t(combo_genos_noNA_32)

rownames(geno1_t) <- branch_name_vec
geno1_dist_mat <- dist(geno1_t, method = 'manhattan', diag = T,
  upper = T)

geno1_nj <- nj(geno1_dist_mat)
geno1_upgma <- hclust(geno1_dist_mat, method = 'average')

geno1_tree_file <- paste(data_dir,
  'analysis_results/figs/tenBranch_sv_dist_trees.png', sep = '')

png(filename = geno1_tree_file, width = 1000, height = 500)
par(mfrow = c(1,2))
plot(geno1_nj, main = 'Ten Branch SV NJ Tree')
plot(geno1_upgma, main = 'Ten Branch SV UPGMA Tree')
dev.off()

# Total number of SVs
geno1_totSV <- nrow(combo_genos_noNA_32)
# 12892

geno1_n_geno_vec <- apply(combo_genos_noNA_32, 1, function(x) length(table(x)))

geno1_var <- combo_genos_noNA_32[which(geno1_n_geno_vec > 1), ]
geno1_varSV <- nrow(geno1_var)
# 583

# Number of each genotype in each sample
geno1_samp_tab_df <- gen_samp_geno_tab_df(geno_mat = geno1_var,
  basic_df = basic_samp_tab_df)

geno1_barplot_info <- gen_samp_geno_barplot_ggs(
  geno_tab_df = geno1_samp_tab_df,
  data_lab = 'Ten Branch SVs')

geno1_combo_bar_file <- paste(data_dir,
  'analysis_results/figs/tenBranch_sv_genotype_bars.png', sep = '')

png(geno1_combo_bar_file, width = 1200, height = 400)
do.call(grid.arrange, c(geno1_barplot_info, ncol = 3))
dev.off()

geno1_singletons <- gen_singleton_df(geno_mat = geno1_var,
  basic_df = basic_singleton_df)

geno1_combo_singleton_bar_file <- paste(data_dir,
  'analysis_results/figs/tenBranch_sv_singleton_bars.png', sep = '')

geno1_singleton_bar_info <- gen_samp_geno_barplot_ggs(
  geno_tab_df = geno1_singletons,
  data_lab = 'Ten Branch SV Singleton Genotypes')

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
# 19

geno1_13_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno1_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_13_lab, meta = samp_meta)

length(geno1_13_spec_hets)
# 1

quit(save = 'n')







