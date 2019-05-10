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
bnd_stat_list <- list()
bnd_stat_list_out <- paste(combo_1_vcf_file, '_BNDstats.txt', sep = '')

# SET VARIABLES
# Maximum distance allowed between the 5' and 3' locations of the receiving
#  portion of the potential INS
bnd_receiv_dist <- info_list[['bnd_receiv_dist']]
# Maximum size of potential INS to be included
insert_max_size <- info_list[['insert_max_size']]

# cutoff for CIPOS length for BND-INS to be considered "good"
cipos_cut <- info_list[['cipos_cut']]

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

bnd_stat_list[['Raw BND count']] <- sum(type_info == 'BND')
bnd_stat_list[['Raw BND count / 4']] <- (sum(type_info == 'BND') / 4)

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

test_bnd_ins <- potential_BND_insertions(vcf_df = raw_vcf, 
  bnd_receiv_dist = bnd_receiv_dist, 
  insert_max_size = insert_max_size, show_progress = T)

bnd_stat_list[['Initial BND-INS count']] <- nrow(test_bnd_ins)

raw_bnd_size_tab <- summary(test_bnd_ins$ins_size)

for(i in seq(length(raw_bnd_size_tab))){
  bnd_stat_list[[paste('Initial BND-INS size', 
    names(raw_bnd_size_tab), sep = ' ')[i]]] <- as.vector(raw_bnd_size_tab)[i]
}

# nrow(test_bnd_ins)
# 1441
# There seem to be 1441 INSERTIONS that are represented by BND entries

## Check/filter by POS confidence interval

ci_5prime_length <- get_cipos_length(raw_vcf, 
  bnd_inds = test_bnd_ins$rec_5prime_ind)

ci_3prime_length <- get_cipos_length(raw_vcf,
  bnd_inds = test_bnd_ins$rec_3prime_ind)

cipos_bad <- union(which(ci_5prime_length > cipos_cut), 
  which(ci_3prime_length > cipos_cut))

# remove BND-INS with excessive CIPOS
bnd_ins_filt1 <- test_bnd_ins[-cipos_bad, ]

bnd_stat_list[['Count of BND-INS with good CIPOS']] <- nrow(bnd_ins_filt1)

good_bnd_size_tab <- summary(bnd_ins_filt1$ins_size)

for(i in seq(length(good_bnd_size_tab))){
  bnd_stat_list[[paste('Good BND-INS size',
    names(good_bnd_size_tab), sep = ' ')[i]]] <- as.vector(good_bnd_size_tab)[i]
}

# 5' genotypes

test_5prime_geno_info <- make_allsamp_geno_info_list(
  raw_vcf[bnd_ins_filt1$rec_5prime_ind, ])

tmp_5prime_genotypes <- call_allsamp_genotypes(
  geno_info_list = test_5prime_geno_info, min_sv_coverage = 10, 
  het_ratio_cut = 0.25, min_minor_allele_count = 2, max_hom_ratio = 0.05)

# need rownames so can easily compare the 5' and 3' genotype matrices
rownames(tmp_5prime_genotypes) <- paste( paste(
  paste('BND_INS', seq(nrow(bnd_ins_filt1)), sep = '_'), 
  'BND', sep = '_'),
  bnd_ins_filt1$ins_size, sep = '_')

tmp_5prime_noNA <- filt_geno_mat(geno_mat = tmp_5prime_genotypes, max_nas = 0)

##
# 3' genotypes

test_3prime_geno_info <- make_allsamp_geno_info_list(
  raw_vcf[bnd_ins_filt1$rec_3prime_ind, ])

tmp_3prime_genotypes <- call_allsamp_genotypes(
  geno_info_list = test_3prime_geno_info, min_sv_coverage = 10, 
  het_ratio_cut = 0.25, min_minor_allele_count = 2, max_hom_ratio = 0.05)

rownames(tmp_3prime_genotypes) <- paste( paste(
  paste('BND_INS', seq(nrow(bnd_ins_filt1)), sep = '_'), 
  'BND', sep = '_'),
  bnd_ins_filt1$ins_size, sep = '_')

tmp_3prime_noNA <- filt_geno_mat(geno_mat = tmp_3prime_genotypes, max_nas = 0)

# make sure have BND-INS that pass for both 5' and 3'
both_side_names <- intersect(rownames(tmp_5prime_noNA), 
  rownames(tmp_3prime_noNA))

bnd_stat_list[['N BND-INS with 5 and 3prime passing genotype filters']] <- (
  length(both_side_names))

unmatch_rownames <- both_side_names[which(
  apply(tmp_5prime_noNA[both_side_names, ], 1, function(x) 
  paste(x, collapse = '')) != 
  apply(tmp_3prime_noNA[both_side_names, ], 1,
  function(x) paste(x, collapse = '')) )]

good_geno_names <- setdiff(both_side_names, unmatch_rownames)

bnd_good_genos <- tmp_5prime_noNA[good_geno_names, ]

bnd_stat_list[['N BND-INS with matching 5 and 3prime genotypes']] <- (
  nrow(bnd_good_genos))
# nrow(bnd_good_genos)
# 43
# 43 BND-INS pass my filtering criteria: good coverage, good quality genotypes,
#     no missing genotypes, 5' and 3' genotypes match

good_geno_size <- as.numeric(unlist(
  lapply(strsplit(good_geno_names, split = '_'), function(x) x[[5]])))

good_geno_size_tab <- summary(good_geno_size)

for(i in seq(length(good_geno_size_tab))){
  bnd_stat_list[[paste('BND-INS with good genotypes  size',
    names(good_geno_size_tab), sep = ' ')[i]]] <- (
    as.vector(good_geno_size_tab)[i])
}

bnd_n_genotypes <- apply(bnd_good_genos, 1, function(x) length(unique(x)))

bnd_stat_list[['N Variable Good BND-INS']] <- sum(bnd_n_genotypes > 1)
# sum(bnd_n_genotypes > 1)
# [1] 0
# no variable filtered BND-INSERTIONS

######
# cnv analysis

cnv_inds <- grep('cnv', type_info)

bnd_stat_list[['Raw cnv count']] <- sum(type_info == 'cnv')
bnd_stat_list[['Shadowed raw cnv count']] <- length(grep('SHADOWED',
  raw_vcf[cnv_inds, 8]))
bnd_stat_list[['Un-shadowed raw cnv count']] <- (sum(type_info == 'cnv') -
  length(grep('SHADOWED', raw_vcf[cnv_inds, 8])))

cnv_shad_inds <- intersect(which(type_info == 'cnv'),
  grep('SHADOWED', raw_vcf[, 8]))
cnv_nonshad_inds <- setdiff(which(type_info == 'cnv'),
  grep('SHADOWED', raw_vcf[, 8]))

tmp_cnv_nonshad_bnd_inds <- cnv_nonshad_inds + 1
cnv_nonshad_bnd_inds <- setdiff(tmp_cnv_nonshad_bnd_inds, cnv_nonshad_inds)

cnv_nonshad_geno_info <- make_allsamp_geno_info_list(
  raw_vcf[cnv_nonshad_bnd_inds, ])

cnv_nonshad_bnd_genos <- call_allsamp_genotypes(
  geno_info_list = cnv_nonshad_geno_info, min_sv_coverage = 10)

cnv_na_count <-apply(cnv_nonshad_bnd_genos, 1, function(x) sum(is.na(x)))

cnv_bnd_noNA_genos <- cnv_nonshad_bnd_genos[which(cnv_na_count == 0), ]

bnd_stat_list[['N cnv POS with adequate coverage']] <- nrow(cnv_bnd_noNA_genos)

########

bnd_stat_df <- data.frame(label = names(bnd_stat_list), 
  value = unlist(bnd_stat_list),
  stringsAsFactors = F)

write.table(bnd_stat_df, file = bnd_stat_list_out, quote = F, sep = '\t',
  row.names = F, col.names = T)

quit(save = 'no')

