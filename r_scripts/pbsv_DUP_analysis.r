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
dup_stat_list <- list()
dup_stat_list_out <- paste(combo_1_vcf_file, '_DUPstats.txt', sep = '')

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

dup_stat_list[['Raw DUP count']] <- sum(type_info == 'DUP')

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

filt_dup_size_tab <- summary(filt_dup_size_vec)

for(i in seq(length(filt_dup_size_tab))){
  dup_stat_list[[paste('Filtered DUP size',
    names(filt_dup_size_tab), sep = ' ')[i]]] <- as.vector(
    filt_dup_size_tab)[i]
}

n_geno_vec <- apply(tmp_dup_genos_noNA, 1, function(x) length(table(x)))

geno1_var <- tmp_dup_genos_noNA[which(n_geno_vec > 1), ]

dup_stat_list[['Number variable Filtered DUP']] <- nrow(geno1_var)

########

dup_stat_df <- data.frame(label = names(dup_stat_list), 
  value = unlist(dup_stat_list),
  stringsAsFactors = F)

write.table(dup_stat_df, file = dup_stat_list_out, quote = F, sep = '\t',
  row.names = F, col.names = T)

quit(save = 'no')

