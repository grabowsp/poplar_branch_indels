# Analysis of Samples to Try to find what sample is closest to 
#  reference used and if there are sample switches

# LOAD LIBRARIES #
function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)

# LOAD DATA #
data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'
combo_file <- 'ref.ALLData.vcf'
combo_file_tot <- paste(data_dir, combo_file, sep = '')
combo_df <- make_combo_indel_df(combo_file_tot)

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')

# SET OUTPUTS #
## distance matrix analysis files
dist_mat_file <- paste(data_dir,
  'analysis_results/poplar_lib_manhattan_distance_50bp.txt', sep = '')
nj_file <- paste(data_dir,
  'analysis_results/figs/sv_dist_50bp_nj_tree.png', sep = '')
upgma_file <- paste(data_dir,
  'analysis_results/figs/sv_dist_50bp_upgma_tree.png', sep = '')

## percent different genotypes
num_SV_genotype_file <- paste(data_dir,
  'analysis_results/percent_SV_genotypes_by_samp.txt', sep = '')

##############
# Filter samples, SVs, and genotypes
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

# indices to remove that are either 70% mononucleotide repeats, 50% the same
#   binucleotide repeat, or 60% of two binucleotide repeats
remove_8mer_SV_inds <- id_prob_SV_seqs(sv_geno_df = combo_df_filt,
  mer_length = 8, per_mn_cutoff = 0.7, per_bn_pure_cutoff = 0.5,
  per_bn_multi_cutoff = 0.6)
combo_filt_1 <- combo_df_filt[-remove_8mer_SV_inds, ]

# remove duplicated positions - the smaller SVs are removed
rm_dup_inds <- dup_inds_to_remove(sv_geno_df = combo_filt_1)
combo_filt_2 <- combo_filt_1[-rm_dup_inds, ]

# remove SV's that are within 30bp of eachother - keep the larger SV
too_close_inds_30 <- overlap_inds_to_remove(sv_geno_df = combo_filt_2,
  dist_cut = 30)
combo_filt_3 <- combo_filt_2[-too_close_inds_30, ]

## find and remove any SVs with missing genotypes in any sample
combo_df_2 <- remove_missing_combo_SVs(combo_filt_3)

###
# Extract genotype information.
combo_geno_info_list <- make_allsamp_geno_info_list(combo_df_2)

# Calling my own genotypes
combo_filtered_genotypes <- call_allsamp_genotypes(combo_geno_info_list,
  min_sv_coverage = 10, het_ratio_cut = 0.15, min_minor_allele_count = 2,
  max_hom_ratio = 0.05)

# remove any SVs with missing data
combo_filt_genos_noNAs <- filt_geno_mat(combo_filtered_genotypes, max_nas = 0)

# remove SV's smaller than 50bp
combo_filt_genos_50bp <- filt_geno_mat(combo_filtered_genotypes, max_nas = 0,
  min_length = 50)

####
# Distance Matrix Analysis to ID possible sample switches
## re-format genotype matrix
filt_geno_t <- t(combo_filtered_genotypes)
## re-order samples
branch_ord <- c()
for(i in seq(nrow(filt_geno_t))){
  tmp_ind <- which(samp_meta$lib_name == rownames(filt_geno_t)[i])
  branch_ord <- c(branch_ord, samp_meta$branch_name[tmp_ind])
}
filt_geno_t_ord <- filt_geno_t[order(branch_ord),]
tmp_names <- paste('branch', sort(branch_ord), sep = '_')
rownames(filt_geno_t_ord) <- tmp_names

# euclidean distance
filt_geno_euc_dist_mat <- dist(filt_geno_t_ord, method = 'euclidean',
  diag = T, upper = T)
# manhattan distance
filt_geno_manhat_dist_mat <- dist(filt_geno_t_ord, method = 'manhattan',
  diag = T, upper = T)

# calculate distance for genotypes filtered for length and missing data
long_geno_t <- t(combo_filt_genos_50bp)
branch_ord <- c()
for(j in seq(nrow(long_geno_t))){
  tmp_ind <- which(samp_meta$lib_name == rownames(long_geno_t)[j])
  branch_ord <- c(branch_ord, samp_meta$branch_name[tmp_ind])
}
long_geno_t_ord <- long_geno_t[order(branch_ord), ]
tmp_names <- paste('branch', sort(branch_ord), sep = '_')
rownames(long_geno_t_ord) <- tmp_names

long_geno_euc_dist_mat <- dist(long_geno_t_ord, method = 'euclidean',
  diag = T, upper = T)
long_geno_manhat_dist_mat <- dist(long_geno_t_ord, method = 'manhattan',
  diag = T, upper = T)

# save the manhattan distance matrix for the filtered SVs
write.table(as.matrix(long_geno_manhat_dist_mat), file = dist_mat_file,
  quote = F, sep = '\t')

# Generate trees
library(ape)

long_nj <- nj(long_geno_manhat_dist_mat)
long_upgma <- hclust(long_geno_manhat_dist_mat, method = 'average')

png(filename = nj_file)
plot(long_nj, main = 'Poplar SV Libs Neighbor Joining')
dev.off()

png(filename = upgma_file)
plot(long_upgma, main = 'Poplar SV Libs UPGMA')
dev.off()

# TakeHome: 14.5 and 13.5 are switched. Should consider 14.5 part of tree 13
#  and 13.5 as part of tree 14
######

# Tally number of SVs to examine relationship with reference; reference
#   sample should have the fewest number of SVs
## Unfiltered SVs and raw PacBio genotype calls
raw_gt_tabs <- lapply(combo_geno_info_list, function(x)
  round(table(x[[1]])/length(x[[1]]), digits = 4))
raw_gt_df <- data.frame(lib = names(raw_gt_tabs),
  raw_per_0 = unlist(lapply(raw_gt_tabs, function(x) x[1])),
  raw_per_1 = unlist(lapply(raw_gt_tabs, function(x) x[2])),
  raw_per_2 = unlist(lapply(raw_gt_tabs, function(x) x[3])),
  stringsAsFactors = F)

raw_gt_df$branch = NA
for(bn in seq(nrow(raw_gt_df))){
  tmp_ind <- which(samp_meta$lib_name == raw_gt_df$lib[bn])
  raw_gt_df$branch[bn] <- samp_meta$branch_name[tmp_ind]
}
raw_gt_df <- raw_gt_df[order(raw_gt_df$branch),c(1,5,2:4)]

# use filtered genotypes
filt_gt_tab <- apply(combo_filtered_genotypes, 2,
  function(x) round(table(x)/sum(table(x)), digits = 4))
raw_gt_df[, c('filt_per_0', 'filt_per_1', 'filt_per_2')] <- NA
for(i in seq(nrow(raw_gt_df))){
  tmp_ind <- which(colnames(filt_gt_tab) == raw_gt_df$lib[i])
  raw_gt_df[i,
  c('filt_per_0', 'filt_per_1', 'filt_per_2')] <- filt_gt_tab[ , tmp_ind]
}

# use SVs filtered for missing data
noNA_gt_df <- apply(combo_filt_genos_noNAs, 2,
  function(x) round(table(x)/sum(table(x)), digits = 4))
raw_gt_df[, c('noNA_per_0', 'noNA_per_1', 'noNA_per_2')] <- NA
for(i in seq(nrow(raw_gt_df))){
  tmp_ind <- which(colnames(noNA_gt_df) == raw_gt_df$lib[i])
  raw_gt_df[i,
  c('noNA_per_0', 'noNA_per_1', 'noNA_per_2')] <- noNA_gt_df[ , tmp_ind]
}

# try SVs larger than 49bp
long_gt_df <- apply(combo_filt_genos_50bp, 2,
  function(x) round(table(x)/sum(table(x)), digits = 4))
raw_gt_df[, c('long_per_0', 'long_per_1', 'long_per_2')] <- NA
for(i in seq(nrow(raw_gt_df))){
  tmp_ind <- which(colnames(long_gt_df) == raw_gt_df$lib[i])
  raw_gt_df[i,
  c('long_per_0', 'long_per_1', 'long_per_2')] <- long_gt_df[ , tmp_ind]
}

write.table(raw_gt_df, file = num_SV_genotype_file, quote = F, sep = '\t',
  row.names = F, col.names = T)

raw_gt_df

#         lib branch raw_per_0 raw_per_1 raw_per_2 filt_per_0 filt_per_1
# PBAU.0 PBAU   13.1    0.0685    0.9124    0.0191     0.1156     0.8749
# PBAW.0 PBAW   13.2    0.0973    0.8880    0.0147     0.1359     0.8552
# PBAT.0 PBAT   13.3    0.0631    0.9182    0.0186     0.1083     0.8823
# PAZF.0 PAZF   13.5    0.0615    0.9210    0.0175     0.0867     0.9047
# PAXN.0 PAXN   14.2    0.0332    0.9470    0.0198     0.0476     0.9414
# PAXL.0 PAXL   14.3    0.0204    0.9592    0.0204     0.0405     0.9508
# PAYK.0 PAYK   14.4    0.0500    0.9299    0.0202     0.0656     0.9263
# PAZH.0 PAZH   14.5    0.0671    0.9146    0.0183     0.1011     0.8911
#        filt_per_2 noNA_per_0 noNA_per_1 noNA_per_2 long_per_0 long_per_1
# PBAU.0     0.0095     0.0620     0.9359     0.0021     0.0139     0.9831
# PBAW.0     0.0088     0.0653     0.9332     0.0015     0.0165     0.9814
# PBAT.0     0.0094     0.0583     0.9396     0.0020     0.0135     0.9838
# PAZF.0     0.0087     0.0420     0.9563     0.0017     0.0096     0.9884
# PAXN.0     0.0109     0.0218     0.9757     0.0025     0.0048     0.9923
# PAXL.0     0.0087     0.0128     0.9851     0.0021     0.0043     0.9932
# PAYK.0     0.0081     0.0300     0.9682     0.0018     0.0084     0.9894
# PAZH.0     0.0078     0.0533     0.9450     0.0018     0.0128     0.9849
#        long_per_2
# PBAU.0     0.0030
# PBAW.0     0.0020
# PBAT.0     0.0027
# PAZF.0     0.0020
# PAXN.0     0.0029
# PAXL.0     0.0025
# PAYK.0     0.0022
# PAZH.0     0.0023

# TakeHome: It looks like Tree 13 has the fewest SVs and is probably the tree
#  that supplied the refrence sample

# OVERALL TAKEHOME: Sample 14.5 and 13.5 were switched. Reference was sample
#  14.5 but actually represents branch 13.5.

quit(save = 'no')



