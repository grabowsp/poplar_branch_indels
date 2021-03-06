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
dup_stat_list <- list()
dup_stat_list_out <- paste(combo_1_vcf_file, '_DUPstats.txt', sep = '')

filt_dup_names_out <- paste(combo_1_vcf_file, '_filtered_DUP_names.rds',
  sep = '')

dup_allvar_out <- paste(combo_1_vcf_file, '_allVarDUP_genos.rds', sep = '')

dup_var_geno_out <- paste(combo_1_vcf_file, '_goodDUP_genos.rds', sep = '')

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

dup_stat_list[['Raw DUP count']] <- sum(type_info == 'DUP')

r_dup_inds <- which(type_info == 'DUP')
r_dup_size_ls <- strsplit(raw_vcf[r_dup_inds, 8], split = ';')
tmp_r_dup_size <- unlist(lapply(r_dup_size_ls, function(x) x[grep('SVLEN', x)]))
r_dup_size <- abs(as.numeric(gsub('SVLEN=', '', tmp_r_dup_size)))

dup_stat_list[['Raw DUP > 100bp']] <- sum(r_dup_size > 100)

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

dup_inds <- which(type_info == 'DUP')

dup_size_list <- strsplit(raw_vcf$INFO[dup_inds], split = ';')
dup_size_vec <- as.numeric(gsub('SVLEN=', '', 
  unlist(lapply(dup_size_list, function(x) x[[length(x)]]))))

dup_size_tab <- summary(dup_size_vec)

for(i in seq(length(dup_size_tab))){
  dup_stat_list[[paste('raw DUP size', 
    names(dup_size_tab), sep = ' ')[i]]] <- as.vector(dup_size_tab)[i]
}

dup_precall_geno_info <- make_allsamp_geno_info_list(raw_vcf[dup_inds, ])

tmp_dup_genotypes <- call_allsamp_genotypes(
  geno_info_list = dup_precall_geno_info, min_sv_coverage = 10,
  het_ratio_cut = 0.25, min_minor_allele_count = 2, max_hom_ratio = 0.05)

rownames(tmp_dup_genotypes) <- paste(raw_vcf$CHROM[dup_inds], 
  raw_vcf$POS[dup_inds], 'DUP', dup_size_vec, sep = '_')

tmp_dup_genos_noNA <- filt_geno_mat(geno_mat = tmp_dup_genotypes, max_nas = 0)

dup_stat_list[['Number Filtered DUP']] <- nrow(tmp_dup_genos_noNA)

filt_dup_size_list <- strsplit(rownames(tmp_dup_genos_noNA), split = '_')
filt_dup_size_vec <- as.numeric(
  unlist(lapply(filt_dup_size_list, function(x) x[[length(x)]])))

filt_dup_names <- rownames(tmp_dup_genos_noNA)

dup_stat_list[['Number Filtered DUP > 100bp']] <- sum(filt_dup_size_vec > 100)
dup_stat_list[['Number Filtered DUP > 1kbp']] <- sum(filt_dup_size_vec > 1000)
dup_stat_list[['Number Filtered DUP > 5kbp']] <- sum(filt_dup_size_vec > 5000)
dup_stat_list[['Number Filtered DUP > 10kbp']] <- sum(filt_dup_size_vec > 10000)
dup_stat_list[['Number Filtered DUP > 25kbp']] <- sum(filt_dup_size_vec > 25000)
dup_stat_list[['Number Filtered DUP > 50kbp']] <- sum(filt_dup_size_vec > 50000)

filt_dup_size_tab <- summary(filt_dup_size_vec)

for(i in seq(length(filt_dup_size_tab))){
  dup_stat_list[[paste('Filtered DUP size',
    names(filt_dup_size_tab), sep = ' ')[i]]] <- as.vector(
    filt_dup_size_tab)[i]
}

n_geno_vec <- apply(tmp_dup_genos_noNA, 1, function(x) length(table(x)))

geno1_var <- tmp_dup_genos_noNA[which(n_geno_vec > 1), ]

if(length(geno1_var > 0)){
  saveRDS(geno1_var, file = dup_allvar_out)
}
# pineapple

dup_stat_list[['Number variable Filtered DUP']] <- nrow(geno1_var)

###################
# Added stuff below from InDel script to try to get "good" variable DUPs

if(nrow(geno1_var) == 0){
  dup_stat_list[['N variable DUPs with loss-of-heterozygosity']] <- 0
  dup_stat_list[['N variable DUPs with gain-of-heterozygosity']] <- 0
  dup_stat_list[['N "decent" variable DUPs']] <- 0
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

dup_stat_list[['N variable DUPs with loss-of-heterozygosity']] <- length(
  lossOfHet_tot)
dup_stat_list[['N variable DUPs with gain-of-heterozygosity']] <- nrow(
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

decent_names <- rownames(geno1_var)[decent_var_SVs]
d_name_split <- strsplit(decent_names, split = '_')
decent_sizes <- as.numeric(unlist(
  lapply(d_name_split, function(x) x[length(x)])))

dup_stat_list[['N "decent" variable DUPs']] <- length(decent_var_SVs)
dup_stat_list[['N "decent" DUPs > 100bp']] <- sum(decent_sizes > 100)

if(length(decent_var_SVs) > 0){
  saveRDS(data.frame(geno1_var)[decent_var_SVs, ], file = dup_var_geno_out)
}

}
########

dup_stat_df <- data.frame(label = names(dup_stat_list), 
  value = unlist(dup_stat_list),
  stringsAsFactors = F)

write.table(dup_stat_df, file = dup_stat_list_out, quote = F, sep = '\t',
  row.names = F, col.names = T)

saveRDS(filt_dup_names, file = filt_dup_names_out)

quit(save = 'no')

