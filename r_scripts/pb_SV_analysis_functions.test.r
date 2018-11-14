# Script for testing pb_SV_analysis_functions.r

data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'
paxl_file <- 'ref.PAXL.vcf'
paxl_file_tot <- paste(data_dir, paxl_file, sep = '')

test_paxl <- load_ind_pbnew_vcf(paxl_file_tot)

test_paxl_raw <- gen_indel_raw_df(test_paxl)

combo_file <- 'ref.ALLData.vcf'
combo_file_tot <- paste(data_dir, combo_file, sep = '')

test_combo <- make_combo_indel_df(combo_file_tot)

test_paxl_2 <- gen_ind_uni_df(indiv_df = test_paxl, combo_df = test_combo)

test_NO_paxl <- gen_nonoverlap_df(uni_ind_df = test_paxl_2, use_svsize = T)

paxl_full_process <- load_to_nonoverlap_df(in_file = paxl_file_tot, 
  combo_df = test_combo)

paxl_indiv_miss <- gen_ind_missing_df(indiv_df = test_paxl, 
  combo_df = test_combo)

test_combo_filt <- remove_missing_combo_SVs(test_combo)

test_ind_geno_info <- make_ind_geno_info_list(test_combo_filt[,12], 
  sv_names = test_combo_filt$full_name_long)

test_full_geno_list <- make_allsamp_geno_info_list(test_combo_filt)

test_ind_num_geno <- call_ind_genotype(test_ind_geno_info, 
  min_sv_coverage = 10, het_ratio_cut = 0.15, min_minor_allele_count = 2, 
  max_hom_ratio = 0.05)

test_allsamp_num_genos <- call_allsamp_genotypes(test_full_geno_list)

test_filt_allsamp_genos <- filt_geno_mat(test_allsamp_num_genos, max_nas = 0)

test_filt_genos_50 <- filt_geno_mat(test_allsamp_num_genos, max_nas = 0, 
  min_length = 50)

test_filt_genos_del <- filt_geno_mat(test_allsamp_num_genos, max_nas = 0,
  sv_type = 'DEL')

test_genos_sing_tab <- tally_singleton_svs(test_filt_allsamp_genos)
test_genos_sing_tab_2 <- tally_singleton_svs(test_filt_genos_50)

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')

test_b13_labels <- c(13.1, 13.2, 13.3, 13.5)
test_b13_libs <- get_branch_lib_names(branch_name_vec = test_b13_labels, 
  meta = samp_meta)

test_b13_cols <- get_lib_col_inds(geno_mat = test_filt_allsamp_genos, 
  lib_name_vec = test_b13_libs)

test_b13_allHet <- get_samegeno_inds(geno_mat = test_filt_allsamp_genos, 
  branch_name_vec = test_b13_labels, genotype = 1, meta = samp_meta)
test_b13_allHomRef <- get_samegeno_inds(geno_mat = test_filt_allsamp_genos, 
  branch_name_vec = test_b13_labels, genotype = 0, meta = samp_meta)


