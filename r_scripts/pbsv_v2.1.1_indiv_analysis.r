# Analysis of individual sample VCF outputs from the pbsv v2.1.1 run

# LOAD LIBRARIES AND PACKAGES
function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)

# LOAD DATA
data_dir <- '/home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL/'
rep_1_vcf_short <- 'Ptr145v1.PAXL_v2_r1.vcf'
rep_1_vcf_file <- paste(data_dir, rep_1_vcf_short, sep = '')

old_vcf_short <- 'ref.PAXL.vcf'
old_vcf_file <- paste(data_dir, old_vcf_short, sep = '')

# SET OUTPUT


# SET VARIABLES


# SET CONSTANTS


############
rep1_in <- load_ind_pbnew_vcf(rep_1_vcf_file)
rep1_indel_vcf <- make_combo_indel_df(rep_1_vcf_file)

rep1_precall_filt <- precalling_SV_filtering(sv_geno_df = rep1_indel_vcf,
  mer_length = 8, per_mn_cutoff = 0.7, per_bn_pure_cutoff = 0.5, 
  per_bn_multi_cutoff = 0.6, dist_cut = 1000, sd_cut = 3, use_sd = T)

rep1_genotypes_1 <- generate_filtered_genotypes(combo_df = rep1_precall_filt,
  min_sv_coverage = 10, het_ratio_cut = 0.25, min_minor_allele_count = 2,
  max_hom_ratio = 0.05, est_error = 0.032, est_penetrance = 0.35,
  good_score = 0.9, great_score = 0.98, same_geno_bonus = 0.67,
  ss_min_cov = 10)

rep1_genos_2 <- rep1_genotypes_1[-grep('scaffold', rownames(rep1_genotypes_1)),]

size_vec <- as.numeric(unlist(lapply(strsplit(names(rep1_genos_2),
  split = '_'), function(x) x[[4]])))

length(grep('INS', names(rep1_genos_2)))
# 17919

length(intersect(grep('INS', names(rep1_genos_2)), which(size_vec >= 100)))
# 7789

length(grep('DEL', names(rep1_genos_2)))
# [1] 17938

length(intersect(grep('DEL', names(rep1_genos_2)), which(size_vec >= 100)))
# [1] 7880

###########
old_indel_vcf <- make_combo_indel_df(old_vcf_file)

old_precall_filt <- precalling_SV_filtering(sv_geno_df = old_indel_vcf,
  mer_length = 8, per_mn_cutoff = 0.7, per_bn_pure_cutoff = 0.5, 
  per_bn_multi_cutoff = 0.6, dist_cut = 1000, sd_cut = 3, use_sd = T)

old_genotypes_1 <- generate_filtered_genotypes(combo_df = old_precall_filt,
  min_sv_coverage = 10, het_ratio_cut = 0.25, min_minor_allele_count = 2,
  max_hom_ratio = 0.05, est_error = 0.032, est_penetrance = 0.35,
  good_score = 0.9, great_score = 0.98, same_geno_bonus = 0.67,
  ss_min_cov = 10)

mn_mer_seqs <- gen_mer_seqs(mer_length = 8, n_bp = 1)
  bn_mer_seqs <- gen_mer_seqs(mer_length = 8, n_bp = 2)

mn_per_mer <- per_mer_length(sv_geno_df = old_indel_vcf, nuc_seqs = mn_mer_seqs)
mn_per_inds <- unlist(per_mer_inds(sv_geno_df = old_indel_vcf,
    per_mer_list = mn_per_mer, per_cutoff = 0.7))

bn_per_mer <- per_mer_length(sv_geno_df = old_indel_vcf, 
   nuc_seqs = bn_mer_seqs)
bn_per_pure_inds <- unlist(per_mer_inds(sv_geno_df = old_indel_vcf,
    per_mer_list = bn_per_mer, per_cutoff = 0.5))

mult_bn_cut <- 0.6 / 2
  tmp_multi_inds <- unlist(per_mer_inds(sv_geno_df = old_indel_vcf,
    per_mer_list = bn_per_mer, per_cutoff = mult_bn_cut))
  tmp_mult_table <- table(unlist(tmp_multi_inds))
  bn_per_multi_inds <- as.numeric(names(tmp_mult_table))[
                         which(tmp_mult_table > 1)]


quit(save = 'n')

