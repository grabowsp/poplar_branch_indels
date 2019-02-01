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

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')
# samp_meta[,c('lib_name', 'branch_name')]

# Generate trees
# Change names to branch names

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

# generate distance matrix and trees
nj_list <- list()
upgma_list <- list()

for(ngl in seq(length(noNA_geno_list))){
  tmp_geno_t <- t(noNA_geno_list[[ngl]])
  tmp_dist_mat <- dist(tmp_geno_t, method = 'manhattan', diag = T, upper = T)
  nj_list[[ngl]] <- nj(tmp_dist_mat)
  upgma_list[[ngl]] <- hclust(tmp_dist_mat, method = 'average')
}

data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'

analysis_res_dir <- paste(data_dir, 'analysis_results/', sep = '')

ts_nj_tree_file <- paste(analysis_res_dir, 
  'PB_pipe_troubleshooting_NJ_tree.png', sep = '')

png(filename = ts_nj_tree_file, width = 1000, height = 1000)
par(mfrow = c(2,2))
for(nji in seq(length(nj_list))){
  plot(nj_list[[nji]], main = paste('Poplar Branch SV NJ tree\nRun', nji))
}
dev.off()

ts_upgma_tree_file <- paste(analysis_res_dir,
  'PB_pipe_troubleshooting_UPGMA_tree.png', sep = '')

png(filename = ts_upgma_tree_file, width = 1000, height = 1000)
par(mfrow = c(2,2))
for(upi in seq(length(upgma_list))){
  plot(upgma_list[[upi]], main = paste('Poplar Branch SV UPGMA tree\nRun', upi))
}
dev.off()

################

# Code to answer some of pbsv developer Aaron's questions

# load data from first round
combo_file <- paste(data_dir, combo_file_name_vec[1], sep = '')

combo_vcf_0 <- read.table(combo_file, sep = '\t', stringsAsFactors = F)
  # add header to data.frame
combo_vcf_header <- scan(combo_file, nlines = 100,
  what = 'character', sep = '\n', quiet = T)
col_info_ind <- grep('#CHROM', combo_vcf_header, fixed = T)
col_info <- gsub('#', '', unlist(strsplit(combo_vcf_header[col_info_ind],
  split = '\t')), fixed = T)
colnames(combo_vcf_0) <- col_info
  #
combo_vcf_0$full_name <- paste(combo_vcf_0[,1], combo_vcf_0[,2], sep = '_')
#  combo_dup_inds <- which(duplicated(combo_vcf_0$full_name))
#  combo_vcf <- combo_vcf_0[-combo_dup_inds,]
combo_vcf <- combo_vcf_0
combo_info_list <- strsplit(combo_vcf[ ,8], split = ';')
combo_sv_vec <- gsub('SVTYPE=', '', unlist(lapply(combo_info_list,
  function(x) x[1])))
combo_vcf$type <- combo_sv_vec

bnd_inds <- which(combo_vcf$type == 'BND')
del_inds <- which(combo_vcf$type == 'DEL')
ins_inds <- which(combo_vcf$type == 'INS')
inv_inds <- which(combo_vcf$type == 'INV')

file_1_geno_list <- make_allsamp_geno_info_list(combo_vcf)

lapply(file_1_geno_list, function(x) sum(x[[1]][bnd_inds] == 1, na.rm = T))
lapply(file_1_geno_list, function(x) sum(table(x[[1]][bnd_inds])))

lapply(file_1_geno_list, function(x) sum(x[[1]][bnd_inds] == 1, na.rm = T)/
  sum(table(x[[1]][bnd_inds])))

svtype_ind_list <- list()
for(sv_type in unique(combo_vcf$type)){
  svtype_ind_list[[sv_type]] <- which(combo_vcf$type == sv_type)
}

sv_het_list <- list()
for(sv_type in names(svtype_ind_list)){
  sv_het_list[[sv_type]][[1]] <- unlist(lapply(file_1_geno_list, function(x) 
    sum(x[[1]][svtype_ind_list[[sv_type]]] == 1, na.rm = T)))
  sv_het_list[[sv_type]][[2]] <- unlist(lapply(file_1_geno_list, function (x)
    sum(table(x[[1]][svtype_ind_list[[sv_type]]]))))
  sv_het_list[[sv_type]][[3]] <- sv_het_list[[sv_type]][[1]] / sv_het_list[[
    sv_type]][[2]]
}

sv_het_mat <- matrix(data = unlist(sv_het_list), nrow = 10, byrow = F)

colnames(sv_het_mat) <- paste(rep(names(sv_het_list), each = 3), 
  c('n_het_genos', 'tot_genos', 'percent_het'), sep = '_')

sv_het_df <- data.frame(lib_name = names(file_1_geno_list), sv_het_mat, 
  stringsAsFactors = F)

analysis_res_dir <- paste(data_dir, 'analysis_results/', sep = '')

run1_sv_het_out <- paste(analysis_res_dir, 'run1_het_by_SVtype.txt', sep = '')

write.table(sv_het_df, file = run1_sv_het_out, quote = F, sep = '\t', 
  row.names = F, col.names = T)

#########
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

noNA_ins_inds <- grep('INS', rownames(noNA_geno_list[[1]]))
noNA_del_inds <- grep('DEL', rownames(noNA_geno_list[[1]]))

no_NA_del_hets <- apply(noNA_geno_list[[1]], 2, 
  function(x) sum(x[noNA_del_inds] == 1))

no_NA_ins_hets <- apply(noNA_geno_list[[1]], 2, 
  function(x) sum(x[noNA_ins_inds] == 1))

noNA_het_tab <- data.frame(lib_name = names(no_NA_del_hets), 
  n_DEL_hets = no_NA_del_hets, n_INS_hets = no_NA_ins_hets, 
  stringsAsFactors = F)

run1_filt_sv_het_out <- paste(analysis_res_dir, 'run1_het_filt_InsDel.txt', 
  sep = '')

write.table(noNA_het_tab, file = run1_filt_sv_het_out, quote = F, sep = '\t',
  row.names = F, col.names = T)

quit(save = 'n')





