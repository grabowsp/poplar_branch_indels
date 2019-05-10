# Analysis specifically of the BND entries from the joint/combined genotype 
#  calling from pbsv - written for v2.1.1 results but SHOULD work for v2.2 
#  results

# LOAD LIBRARIES AND PACKAGES
function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)
# library(Hmisc)

# LOAD DATA
args <- commandArgs(trailingOnly = T)

info_file <- args[1]
info_list <- readRDS(info_file)

data_dir <- info_list[['data_dir']]
combo_1_vcf_short <- info_list[['vcf_short']]

# data_dir <- '/home/f1p1/tmp/PBSV/Poplar14.5/LIBS/'
# combo_1_vcf_short <- 'Ptr145v1.ALLData.v2.try1.vcf'

combo_1_vcf_file <- paste(data_dir, combo_1_vcf_short, sep = '')

meta_in <- info_list[['meta_in']]
#meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')

# SET OUTPUT
# General stats about the BND results
inv_stat_list <- list()
inv_stat_list_out <- paste(combo_1_vcf_file, '_INVstats.txt', sep = '')

# SET VARIABLES
# for now, don't need this info
#lib_order <- c('PAXL', 'PAXN', 'PAYK', 'PAYZ', 'PAZF', 'PAZG', 'PAZH', 'PBAT',
#  'PBAU', 'PBAW')

# for now, don't need this info
#branch_13_lab <- c(13.4, 13.5, 13.3, 13.2, 13.1)
#branch_14_lab <- c(14.5, 14.1, 14.4, 14.3, 14.2)

# SET CONSTANTS


############
# Raw Counts
raw_vcf <- read.table(combo_1_vcf_file, header = F, stringsAsFactors = F,
  sep = '\t')

raw_info <- strsplit(gsub('IMPRECISE;', '', raw_vcf[,8]), split = ';')
type_info <- unlist(lapply(raw_info, 
  function(x) unlist(strsplit(x[[1]], split = '='))[2]))

raw_type_tab <- table(type_info)
#   BND   DEL   INS   INV 
# 28432 29089 21618    31

inv_stat_list[['Raw INV count']] <- sum(type_info == 'INV')

# add column names to VCF and "full_name" that is required for generating
#   genotypes
test_vcf_header <- scan(combo_1_vcf_file, nlines = 100, what = 'character',
  sep = '\n', quiet = T)
test_vcf_info_ind <- grep('#CHROM', test_vcf_header, fixed = T)
test_vcf_col_info <-  gsub('#', '',
  unlist(strsplit(test_vcf_header[test_vcf_info_ind],
  split = '\t')), fixed = T)
colnames(raw_vcf) <- test_vcf_col_info
raw_vcf$full_name <- paste(raw_vcf[,1], raw_vcf[,2], sep = '_')

#######

# bnd_inds <- which(type_info == 'BND')

# Find likely INSERTIONS that are represented by BND entries

inv_inds <- which(type_info == 'INV')

inv_end_vec <- as.numeric(gsub('SVTYPE=INV;END=', '', raw_vcf$INFO[inv_inds]))

inv_length_vec <- inv_end_vec - raw_vcf$POS[inv_inds]

raw_inv_size_tab <- summary(inv_length_vec)

for(i in seq(length(raw_inv_size_tab))){
  inv_stat_list[[paste('raw INV size', 
    names(raw_inv_size_tab), sep = ' ')[i]]] <- as.vector(raw_inv_size_tab)[i]
}

inv_precall_geno_info <- make_allsamp_geno_info_list(raw_vcf[inv_inds, ])

tmp_inv_genotypes <- call_allsamp_genotypes(
  geno_info_list = inv_precall_geno_info, min_sv_coverage = 10,
  het_ratio_cut = 0.25, min_minor_allele_count = 2, max_hom_ratio = 0.05)

rownames(tmp_inv_genotypes) <- paste(raw_vcf$CHROM[inv_inds], 
  raw_vcf$POS[inv_inds], 'INV', sep = '_')

tmp_inv_genos_noNA <- filt_geno_mat(geno_mat = tmp_inv_genotypes, max_nas = 0)

inv_stat_list[['Number Filtered INV']] <- nrow(tmp_inv_genos_noNA)

n_geno_vec <- apply(tmp_inv_genos_noNA, 1, function(x) length(table(x)))

geno1_var <- tmp_inv_genos_noNA[which(n_geno_vec > 1), ]

inv_stat_list[['Number variable Filtered INV']] <- nrow(geno1_var)

########

inv_stat_df <- data.frame(label = names(inv_stat_list), 
  value = unlist(inv_stat_list),
  stringsAsFactors = F)

write.table(inv_stat_df, file = inv_stat_list_out, quote = F, sep = '\t',
  row.names = F, col.names = T)

quit(save = 'no')

