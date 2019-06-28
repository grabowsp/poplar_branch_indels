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
res_info_df <- read.table(info_file, header = F, stringsAsFactors = F,
  sep = '\t', row.names = 1)

combo_1_vcf_file <- trimws(res_info_df['vcf_full', 1])

res_file_parts <- unlist(strsplit(combo_1_vcf_file, split = '/'))
data_dir <- paste(
  paste(res_file_parts[-length(res_file_parts)], collapse = '/'),
  '/', sep = '')
vcf_short <- res_file_parts[length(res_file_parts)]
#data_dir <- '/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/'
#vcf_short <- 'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r01.vcf'

meta_in <- trimws(res_info_df['meta_in', 1])
#meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')

# SET OUTPUT
# General stats about the BND results
inv_stat_list <- list()
inv_stat_list_out <- paste(combo_1_vcf_file, '_INVstats.txt', sep = '')

filt_inv_names_out <- paste(combo_1_vcf_file, '_filtered_INV_names.rds',
  sep = '')

# SET VARIABLES
lib_order <- trimws(unlist(strsplit(res_info_df['lib_order',1], split = ',')))
#lib_order <- c('PAXL', 'PAXN', 'PAYK', 'PAYZ', 'PAZF', 'PAZG', 'PAZH', 'PBAT',
#  'PBAU', 'PBAW')

branch_14_lab <- as.numeric(
  trimws(unlist(strsplit(res_info_df['branch_13_lab',1], split = ','))))
branch_13_lab <- as.numeric(
  trimws(unlist(strsplit(res_info_df['branch_14_lab',1], split = ','))))

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
  raw_vcf$POS[inv_inds], inv_length_vec, 'INV', sep = '_')

tmp_inv_genos_noNA <- filt_geno_mat(geno_mat = tmp_inv_genotypes, max_nas = 0)

filt_inv_names <- rownames(tmp_inv_genos_noNA)

inv_stat_list[['Number Filtered INV']] <- nrow(tmp_inv_genos_noNA)

filt_size_vec <- as.numeric(unlist(lapply(
  strsplit(rownames(tmp_inv_genos_noNA), split = '_'), function(x) x[[3]])))

inv_stat_list[['N Filtered INV > 100bp']] <- sum(filt_size_vec > 100)
inv_stat_list[['N Filtered INV > 1kbp']] <- sum(filt_size_vec > 1000)
inv_stat_list[['N Filtered INV > 5kbp']] <- sum(filt_size_vec > 5000)
inv_stat_list[['N Filtered INV > 10kbp']] <- sum(filt_size_vec > 10000)
inv_stat_list[['N Filtered INV > 25kbp']] <- sum(filt_size_vec > 25000)
inv_stat_list[['N Filtered INV > 50kbp']] <- sum(filt_size_vec > 50000)

n_geno_vec <- apply(tmp_inv_genos_noNA, 1, function(x) length(table(x)))

geno1_var <- tmp_inv_genos_noNA[which(n_geno_vec > 1), ]

inv_stat_list[['Number variable Filtered INV']] <- nrow(geno1_var)

###################
# Added stuff below from InDel script to try to get "good" variable INVs

if(nrow(geno1_var) == 0){
  inv_stat_list[['N variable INVs with loss-of-heterozygosity']] <- 0
  inv_stat_list[['N variable INVs with gain-of-heterozygosity']] <- 0
  inv_stat_list[['N "decent" variable INVs']] <- 0
} else{

# Find variable DUP inds variable in either or both clones
tmp_13_libs <- c()
for(bl in branch_13_lab){
  tmp_ml_ind <- which(samp_meta$branch_name == bl)
  tmp_13_libs <- c(tmp_13_libs, samp_meta$lib_name[tmp_ml_ind])
}

t13_cols <- which(colnames(geno1_var) %in% tmp_13_libs)
t13_ngenos <- apply(geno1_var[, t13_cols], 1, function(x) length(unique(x)))
t13_var_inds <- which(t13_ngenos > 1)

tmp_14_libs <- c()
for(bl in branch_14_lab){
  tmp_ml_ind <- which(samp_meta$branch_name == bl)
  tmp_14_libs <- c(tmp_14_libs, samp_meta$lib_name[tmp_ml_ind])
}

t14_cols <- which(colnames(geno1_var) %in% tmp_14_libs)
t14_ngenos <- apply(geno1_var[, t14_cols], 1, function(x) length(unique(x)))
t14_var_inds <- which(t14_ngenos > 1)

# Look for loss-of-heterozygosity - find where homozygosity is found in
#  branch(es) above branches that are heterozygous
low_13_het <- apply(geno1_var[ , t13_cols], 1, function(x) min(which(x == 1)))
hi_13_hom <- apply(geno1_var[ , t13_cols], 1, function(x)
  max(union(which(x == 0), which(x == 2))))

soloHet_13 <- which(apply(geno1_var[, t13_cols], 1, function(x) sum(x == 1))
  == 1)

low_13_change <- low_13_het - hi_13_hom
lossOfHet_13 <- setdiff(which(low_13_change < 0 & low_13_change > -Inf),
  soloHet_13)

low_14_het <- apply(geno1_var[ , t14_cols], 1, function(x) min(which(x == 1)))
hi_14_hom <- apply(geno1_var[ , t14_cols], 1, function(x)
  max(union(which(x == 0), which(x == 2))))

soloHet_14 <- which(apply(geno1_var[, t14_cols], 1, function(x) sum(x == 1))
  == 1)

low_14_change <- low_14_het - hi_14_hom
lossOfHet_14 <- setdiff(which(low_14_change < 0 & low_14_change > -Inf),
  soloHet_14)

lossOfHet_tot <- union(lossOfHet_13, lossOfHet_14)

inv_stat_list[['N variable INVs with loss-of-heterozygosity']] <- length(
  lossOfHet_tot)
inv_stat_list[['N variable INVs with gain-of-heterozygosity']] <- nrow(
  geno1_var) - length(lossOfHet_tot)

# find SVs that are fixed for 1 in one clone and variable in the other clone
t13_fixed_het <- which(apply(geno1_var[, t13_cols], 1,
  function(x) sum(x == 1)) == length(t13_cols))

t13fixHet_t14Var <- intersect(t13_fixed_het, t14_var_inds)

t14_fixed_het <- which(apply(geno1_var[, t14_cols], 1,
  function(x) sum(x == 1)) == length(t14_cols))

t14fixHet_t13Var <- intersect(t14_fixed_het, t13_var_inds)

bad_var_SVs <- union(lossOfHet_tot,
  union(union(t13fixHet_t14Var, t14fixHet_t13Var),
  intersect(t13_var_inds, t14_var_inds)))

decent_var_SVs <- setdiff(seq(nrow(geno1_var)), bad_var_SVs)

inv_stat_list[['N "decent" variable INVs']] <- length(decent_var_SVs)

}
########

inv_stat_df <- data.frame(label = names(inv_stat_list), 
  value = unlist(inv_stat_list),
  stringsAsFactors = F)

write.table(inv_stat_df, file = inv_stat_list_out, quote = F, sep = '\t',
  row.names = F, col.names = T)

saveRDS(filt_inv_names, file = filt_inv_names_out)

quit(save = 'no')

