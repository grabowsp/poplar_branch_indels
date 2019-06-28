# Analysis of the joint/combined genotype calling from pbsv v2.1.1

# LOAD LIBRARIES AND PACKAGES
function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)
library(Hmisc)

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
## General Stat list
stat_list <- list()
stat_list_out <- paste(combo_1_vcf_file, '_SVstats.txt', sep = '')

# dataframe containing tallys of variable INDELs
# samp_info_df
samp_info_out <- paste(combo_1_vcf_file, '_sampSVInfo.txt', sep = '')

# dataframe containing the number of singletons of each SV type and
#   genotype class for each sample
# samp_sing_df
samp_sing_out <- paste(combo_1_vcf_file, '_sampSingletonInfo.txt', sep = '')

samp_ord_corr_out <- paste(combo_1_vcf_file, '_varSVsampOrderCorr.txt', 
  sep = '')

sing_ord_corr_out <- paste(combo_1_vcf_file, '_singletonSampOrderCorr.txt', 
  sep = '')

filt_indel_names_out <- paste(combo_1_vcf_file, '_filtered_INDEL_names.rds',
  sep = '')

all_var_genos_out <- paste(combo_1_vcf_file, '_allVarINDEL_genos.rds', sep = '')

good_var_genos_out <- paste(combo_1_vcf_file, '_goodINDEL_genos.rds', sep = '')

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

for(i in seq(length(raw_type_tab))){
  stat_list[[paste('raw', names(raw_type_tab), sep = ' ')[i]]] <- as.vector(
    raw_type_tab)[i]
}

bnd_ins_equiv <- raw_type_tab[names(raw_type_tab) == 'BND'] / 4

stat_list[['BND count / 4 (possible INS)']] <- bnd_ins_equiv

r_del_inds <- which(type_info == 'DEL')
r_del_size_ls <- strsplit(raw_vcf[r_del_inds, 8], split = ';')
tmp_r_del_size <- unlist(lapply(r_del_size_ls, function(x) x[grep('SVLEN', x)]))
r_del_size <- abs(as.numeric(gsub('SVLEN=', '', tmp_r_del_size)))

stat_list[['Raw DEL > 100bp']] <- sum(r_del_size > 100)

r_ins_inds <- which(type_info == 'INS')
r_ins_size_ls <- strsplit(raw_vcf[r_ins_inds, 8], split = ';')
tmp_r_ins_size <- unlist(lapply(r_ins_size_ls, function(x) x[grep('SVLEN', x)]))
r_ins_size <- abs(as.numeric(gsub('SVLEN=', '', tmp_r_ins_size)))

stat_list[['Raw INS > 100bp']] <- sum(r_ins_size > 100)

# Filtered InDel Analysis
combo1_vcf <- make_combo_indel_df(combo_1_vcf_file)

combo1_precall_df <- precalling_SV_filtering(sv_geno_df = combo1_vcf,
  mer_length = 8, per_mn_cutoff = 0.7, per_bn_pure_cutoff = 0.5,
  per_bn_multi_cutoff = 0.6, dist_cut = 1000, sd_cut = 3, use_sd = T)

combo1_genotypes_1 <- generate_filtered_genotypes(combo_df = combo1_precall_df,
  min_sv_coverage = 10, het_ratio_cut = 0.25, min_minor_allele_count = 2,
  max_hom_ratio = 0.05, est_error = 0.032, est_penetrance = 0.35,
  good_score = 0.9, great_score = 0.98, same_geno_bonus = 0.67,
  ss_min_cov = 20)

combo1_genos_noNA <- filt_geno_mat(geno_mat = combo1_genotypes_1,
  max_nas = 0, min_length = 15)

#nrow(combo1_genos_noNA)
# [1] 12646

stat_list[['Total number of filtered INDELs']] <- nrow(combo1_genos_noNA)

scaf_inds <- grep('scaffold', rownames(combo1_genos_noNA))
# length(scaf_inds)
# [1] 3

stat_list[['Filtered INDELs in scaffolds']] <- length(scaf_inds)

chr_inds <- grep('Chr', rownames(combo1_genos_noNA))
# length(chr_inds)
# [1] 12643

stat_list[['Filtered INDELs in chromosomes']] <- length(chr_inds)

combo1_genos_2 <- combo1_genos_noNA[-scaf_inds, ]
combo1_genos_info <- strsplit(rownames(combo1_genos_2), split = '_')
combo1_genos_type <- unlist(lapply(combo1_genos_info, function(x) x[[3]]))

filt_indel_names <- rownames(combo1_genos_2)
# table(combo1_genos_type)
#  DEL  INS 
# 7647 4996 

combo1_g2_del_inds <- grep('DEL', combo1_genos_type)
combo1_g2_ins_inds <- grep('INS', combo1_genos_type)

stat_list[['Number filtered Deletions']] <- length(combo1_g2_del_inds)

stat_list[['Number filtered Insertions']] <- length(combo1_g2_ins_inds)

# Start generating chromosome-level data
chrom_info_df <- data.frame(
  chrom = paste('Chr', sprintf('%02d', c(1:19)), sep = ''), 
  stringsAsFactors = F)

add_chr_tab_info <- function(chr_table, chr_df, new_lab){
  # Function to add chromosome-tally data from tables to a combined dataframe
  # INPUTS
  # chr_table = table for values tallied up by chromosome
  # chr_df = the dataframe that want to add the table info to
  # new_lab = column label for the table values added to the dataframe
  # OUTPUT
  # dataframe containing new column with values from chr_table
  ###################
  chrs_in_tab <- which(chr_df$chrom %in% names(chr_table))
  chr_df[, new_lab] <- 0
  chr_df[chrs_in_tab, new_lab] <- as.vector(chr_table) 
  return(chr_df)
}

####

combo1_genos_chr <- unlist(lapply(combo1_genos_info, function(x) x[[1]]))

c1_del_chr_tab <- table(combo1_genos_chr[combo1_g2_del_inds])
# Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 
#   867   529   417   522   532   531   334   366   254   462   399   313   380 
# Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 
#   337   235   308   308   290   263

chrom_info_df <- add_chr_tab_info(chr_table = c1_del_chr_tab, 
  chr_df = chrom_info_df, 
  new_lab = 'n_DEL')

c1_ins_chr_tab <- table(combo1_genos_chr[combo1_g2_ins_inds])
# Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 
#   593   358   300   328   367   375   209   227   213   321   204   162   209 
# Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 
#   225   165   208   188   213   131

chrom_info_df <- add_chr_tab_info(chr_table = c1_ins_chr_tab, 
  chr_df = chrom_info_df, 
  new_lab = 'n_INS')

combo1_genos_size <- unlist(lapply(combo1_genos_info, 
  function(x) as.numeric(x[[4]])))
combo1_g2_100bp_inds <- which(combo1_genos_size >= 100)

stat_list[['Number INDELs > 100bp']] <- length(combo1_g2_100bp_inds)

stat_list[['Number DELs > 100bp']] <- length(
  intersect(combo1_g2_100bp_inds, combo1_g2_del_inds))
# [1] 2562

stat_list[['Number INSs > 100bp']] <- length(
  intersect(combo1_g2_100bp_inds, combo1_g2_ins_inds))
# [1] 2593

combo1_g2_1kbp_inds <- which(combo1_genos_size >= 1000)
combo1_g2_5kbp_inds <- which(combo1_genos_size >= 5000)
combo1_g2_10kbp_inds <- which(combo1_genos_size >= 10000)
combo1_g2_25kbp_inds <- which(combo1_genos_size >= 25000)
combo1_g2_50kbp_inds <- which(combo1_genos_size >= 50000)

stat_list[['Number DELs > 1kbp']] <- length(
  intersect(combo1_g2_1kbp_inds, combo1_g2_del_inds))

stat_list[['Number INSs > 1kbp']] <- length(
  intersect(combo1_g2_1kbp_inds, combo1_g2_ins_inds))
#
stat_list[['Number DELs > 5kbp']] <- length(
  intersect(combo1_g2_5kbp_inds, combo1_g2_del_inds))

stat_list[['Number INSs > 5kbp']] <- length(
  intersect(combo1_g2_5kbp_inds, combo1_g2_ins_inds))
#
stat_list[['Number DELs > 10kbp']] <- length(
  intersect(combo1_g2_10kbp_inds, combo1_g2_del_inds))

stat_list[['Number INSs > 10kbp']] <- length(
  intersect(combo1_g2_10kbp_inds, combo1_g2_ins_inds))
#
stat_list[['Number DELs > 25kbp']] <- length(
  intersect(combo1_g2_25kbp_inds, combo1_g2_del_inds))

stat_list[['Number INSs > 25kbp']] <- length(
  intersect(combo1_g2_25kbp_inds, combo1_g2_ins_inds))
#
stat_list[['Number DELs > 50kbp']] <- length(
  intersect(combo1_g2_50kbp_inds, combo1_g2_del_inds))

stat_list[['Number INSs > 50kbp']] <- length(
  intersect(combo1_g2_50kbp_inds, combo1_g2_ins_inds))

#

c1_del_100_chr_tab <- table(
  combo1_genos_chr[intersect(combo1_g2_100bp_inds, combo1_g2_del_inds)])
# Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 
#   296   191   136   177   175   192   116   128    95   155   128   101   114 
# Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 
#   110    75   104   101    90    78 

chrom_info_df <- add_chr_tab_info(chr_table = c1_del_100_chr_tab,
  chr_df = chrom_info_df,
  new_lab = 'n_DEL_100bp')

c1_ins_100_chr_tab <- table(
  combo1_genos_chr[intersect(combo1_g2_100bp_inds, combo1_g2_ins_inds)])
# Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 
#   303   186   161   168   191   175   109   133   121   176   107    74   113 
# Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 
#   103    86    96   100   118    73

chrom_info_df <- add_chr_tab_info(chr_table = c1_ins_100_chr_tab,
  chr_df = chrom_info_df,
  new_lab = 'n_INS_100bp')

# Variable INDEL analysis

n_geno_vec <- apply(combo1_genos_2, 1, function(x) length(table(x)))

geno1_var <- combo1_genos_2[which(n_geno_vec > 1), ]

# pineapple

if(length(geno1_var > 0)){
  saveRDS(geno1_var, file = all_var_genos_out)
}

stat_list[['Number variable INDELs']] <- nrow(geno1_var)
# 160

geno1_var_info <- strsplit(rownames(geno1_var), split = '_')
geno1_var_type <- unlist(lapply(geno1_var_info, function(x) x[[3]]))

stat_list[['Number variable DELs']] <- sum(geno1_var_type == 'DEL')
# 77

stat_list[['Number variable INSs']] <- sum(geno1_var_type == 'INS')
# 83

geno1_var_size <- unlist(lapply(geno1_var_info, function(x) as.numeric(x[[4]])))

stat_list[['Number variable DELs > 100bp']] <- length(
  intersect(which(geno1_var_type == 'DEL'), which(geno1_var_size >= 100)))
# 62

stat_list[['Number variable INSs > 100bp']] <-  length(
  intersect(which(geno1_var_type == 'INS'), which(geno1_var_size >= 100)))
# 28

# number of invariant Deletions with each genotype class
stat_list[['Number homREF invariant DELs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1), 
  combo1_g2_del_inds), 1] == 0)
stat_list[['Percent homREF invariant DELs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1), 
  combo1_g2_del_inds), 1] == 0)/length(intersect(which(n_geno_vec == 1),
  combo1_g2_del_inds))

stat_list[['Number homALT invariant DELs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1), 
  combo1_g2_del_inds), 1] == 2)
stat_list[['Percent homALT invariant DELs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1), 
  combo1_g2_del_inds), 1] == 2)/length(intersect(which(n_geno_vec == 1),
  combo1_g2_del_inds))

stat_list[['Number HET invariant DELs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1),
  combo1_g2_del_inds), 1] == 1)
stat_list[['Percent HET invariant DELs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1),
  combo1_g2_del_inds), 1] == 1)/length(intersect(which(n_geno_vec == 1),
  combo1_g2_del_inds))

# invariant Insertions with each genotype class
stat_list[['Number homREF invariant INSs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1),
  combo1_g2_ins_inds), 1] == 0)
stat_list[['Percent homREF invariant INSs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1),
  combo1_g2_ins_inds), 1] == 0)/length(intersect(which(n_geno_vec == 1),
  combo1_g2_ins_inds))

stat_list[['Number homALT invariant INSs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1),
  combo1_g2_ins_inds), 1] == 2)
stat_list[['Percent homALT invariant INSs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1),
  combo1_g2_ins_inds), 1] == 2)/length(intersect(which(n_geno_vec == 1),
  combo1_g2_ins_inds))

stat_list[['Number HET invariant INSs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1),
  combo1_g2_ins_inds), 1] == 1)
stat_list[['Percent HET invariant INSs']] <- sum(
  combo1_genos_2[intersect(which(n_geno_vec == 1),
  combo1_g2_ins_inds), 1] == 1)/length(intersect(which(n_geno_vec == 1),
  combo1_g2_ins_inds))

# Number of variable INDELs with each genotype class
del1_var <- geno1_var[which(geno1_var_type == 'DEL'), ]
del1_var_100 <- geno1_var[intersect(which(geno1_var_type == 'DEL'), 
  which(geno1_var_size >= 100)), ]

ins1_var <- geno1_var[which(geno1_var_type == 'INS'), ]
ins1_var_100 <- geno1_var[intersect(which(geno1_var_type == 'INS'),
  which(geno1_var_size >= 100)), ]

## Deletions
stat_list[['Number of variant DELs with homREF genotypes']] <- sum(
  apply(del1_var, 1, function(x) 0 %in% x))
stat_list[['Percent variant DELs with homREF genotypes']] <- sum(
  apply(del1_var, 1, function(x) 0 %in% x))/nrow(del1_var)

stat_list[['Number of variant DELs with homALT genotypes']] <- sum(
  apply(del1_var, 1, function(x) 2 %in% x))
stat_list[['Percent variant DELs with homALT genotypes']] <- sum(
  apply(del1_var, 1, function(x) 2 %in% x))/nrow(del1_var)

stat_list[['Number of variant DELs with HET genotypes']] <- sum(
  apply(del1_var, 1, function(x) 1 %in% x))
stat_list[['Percent variant DELs with HET genotypes']] <- sum(
  apply(del1_var, 1, function(x) 1 %in% x))/nrow(del1_var)

# variable Insertions with each genotype class
stat_list[['Number of variant INSs with homREF genotypes']] <- sum(
  apply(ins1_var, 1, function(x) 0 %in% x))
stat_list[['Percent variant INSs with homREF genotypes']] <- sum(
  apply(ins1_var, 1, function(x) 0 %in% x))/nrow(ins1_var)

stat_list[['Number of variant INSs with homALT genotypes']] <- sum(
  apply(ins1_var, 1, function(x) 2 %in% x))
stat_list[['Percent variant INSs with homALT genotypes']] <- sum(
  apply(ins1_var, 1, function(x) 2 %in% x))/nrow(ins1_var)

stat_list[['Number of variant INSs with HET genotypes']] <- sum(
  apply(ins1_var, 1, function(x) 1 %in% x))
stat_list[['Percent variant INSs with HET genotypes']] <- sum(
  apply(ins1_var, 1, function(x) 1 %in% x))/nrow(ins1_var)

# Start sample-level counts

samp_info_df <- data.frame(samp = colnames(geno1_var), stringsAsFactors = F)

samp_info_df$samp_order <- NA
#samp_order_ind <- c()
for(i in seq(length(lib_order))){
  tmp_match <- which(samp_info_df$samp == lib_order[i])
  samp_info_df$samp_order[tmp_match] <- i
}

samp_info_df$n_var_INDEL_het <- apply(geno1_var, 2, function(x) sum(x == 1))
# PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW 
#   92  105  101   96   99  102  108  108  119  136

samp_info_df$n_var_DEL_het <- apply(del1_var, 2, function(x) sum(x == 1))
# PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW 
#   65   72   70   65   65   68   69   64   64   67

samp_info_df$n_var_DEL_100bp_het <- apply(del1_var_100, 2, 
  function(x) sum(x == 1))
# PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW 
#   51   57   56   52   52   55   55   49   51   56

samp_info_df$n_var_INS_het <- apply(ins1_var, 2, function(x) sum(x == 1))
# PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW 
#   27   33   31   31   34   34   39   44   55   69

samp_info_df$n_var_INS_100bp_het <- apply(ins1_var_100, 2, 
  function(x) sum(x == 1))
# PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW 
#   11   18   17   14   18   15   17   18   20   21

#   look at corrleation between order and number of variable INDELS

corr_info <- rcorr(as.matrix(samp_info_df[, c(2:ncol(samp_info_df))]))

samp_order_r <- corr_info[[1]]['samp_order', -1]
samp_order_pval <- corr_info[[3]]['samp_order', -1]

samp_ord_df <- data.frame('sv_count' = names(samp_order_r), 
  'r_samp_order' = samp_order_r, 'cor_pval' = samp_order_pval, 
  stringsAsFactors = F)

# banana

### Singleton Analysis
basic_singleton_df <- data.frame(lib = colnames(del1_var), homRef = 0, het = 0,
 homAlt = 0, stringsAsFactors = F)

del1_singletons <- gen_singleton_df(geno_mat = del1_var,
  basic_df = basic_singleton_df)

# apply(del1_singletons[,c(2:4)], 2, sum)
# homRef    het homAlt 
#     30      0     34

stat_list[['Number singleton DELs']] <- sum(apply(
  del1_singletons[,c(2:4)], 2, sum))
# [1] 64
# 64 of the 77 variable deletions are singletons - troubling

stat_list[['Perecent singleton variable  DELs']] <- sum(apply(
  del1_singletons[,c(2:4)], 2, sum))/nrow(del1_var)
# [1] 0.8311688
 
stat_list[['Number singleton HET DELs']] <- sum(del1_singletons$het)

del1_100_singletons <- gen_singleton_df(geno_mat = del1_var_100,
  basic_df = basic_singleton_df)

# apply(del1_100_singletons[,c(2:4)], 2, sum)
# homRef    het homAlt 
#     21      0     28

stat_list[['Number singleton DELs > 100bp']] <- sum(apply(
  del1_100_singletons[,c(2:4)], 2, sum))
# [1] 49

stat_list[['Percent singleton variable DELs > 100bp']] <- sum(apply(
  del1_100_singletons[,c(2:4)], 2, sum))/nrow(del1_var_100)
# [1] 0.7903226

stat_list[['Number singleton HET DELs > 100bp']] <- sum(del1_100_singletons$het)

# INSERTIONS
ins1_singletons <- gen_singleton_df(geno_mat = ins1_var,
  basic_df = basic_singleton_df)

# apply(ins1_singletons[,c(2:4)], 2, sum)
# homRef    het homAlt 
#    21     27     12

stat_list[['Number singleton INSs']] <- sum(apply(
  ins1_singletons[,c(2:4)], 2, sum))
# [1] 60
# 60 of 83 variable insertions are singletons
stat_list[['Percent singleton variable INSs']] <- sum(apply(
  ins1_singletons[,c(2:4)], 2, sum))/nrow(ins1_var)
# [1] 0.7228916
stat_list[['Number singleton HET INSs']] <- sum(ins1_singletons$het)

ins1_100_singletons <- gen_singleton_df(geno_mat = ins1_var_100,
  basic_df = basic_singleton_df)

# apply(ins1_100_singletons[,c(2:4)], 2, sum)
# homRef    het homAlt 
#    13     8     4

stat_list[['Number singleton INSs > 100bp']] <- sum(apply(
  ins1_100_singletons[,c(2:4)], 2, sum))
# [1] 25
stat_list[['Percent singleton variable INSs > 100bp']] <- sum(apply(
  ins1_100_singletons[,c(2:4)], 2, sum))/nrow(ins1_var_100)
# 0.8928571
stat_list[['Number of singleton HET INSs > 100bp']] <- sum(
  ins1_100_singletons$het)

# make dataframe that contains the number of singleton genotypes for
#  each sample
sing_list <- list()
sing_list[['DEL_20bp']] <- del1_singletons
sing_list[['DEL_100bp']] <- del1_100_singletons
sing_list[['INS_20bp']] <- ins1_singletons
sing_list[['INS_100bp']] <- ins1_100_singletons

samp_sing_df <- data.frame(lib = del1_singletons$lib, stringsAsFactors = F)
# artichoke

samp_sing_df$samp_order <- NA
for(i in seq(length(lib_order))){
  tmp_match <- which(samp_sing_df$lib == lib_order[i])
  samp_sing_df$samp_order[tmp_match] <- i
}

for(sl in seq(length(sing_list))){
  type_name <- names(sing_list)[sl]
  tmp_colnames <- paste(type_name, colnames(sing_list[[sl]])[2:4], sep = '_')
  samp_sing_df[, tmp_colnames] <- sing_list[[sl]][, c(2:4)]
}

# look at correlation between sample order and number of singletons
sing_corr_info <- rcorr(as.matrix(samp_sing_df[, c(2:ncol(samp_sing_df))]))

sing_order_r <- sing_corr_info[[1]]['samp_order', -1]
sing_order_pval <- sing_corr_info[[3]]['samp_order', -1]

sing_ord_df <- data.frame('sv_count' = names(sing_order_r), 
  'r_samp_order' = sing_order_r, 'cor_pval' = sing_order_pval, 
  stringsAsFactors = F)

# Clone specific SVs
geno1_14_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno1_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_14_lab, meta = samp_meta)
# integer(0) - none

stat_list[['Number Tree 14 fixed INDELs']] <- length(geno1_14_spec_hets)

geno1_13_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno1_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_13_lab, meta = samp_meta)
# integer(0) - none
# I'm a bit worried about these results...

stat_list[['Number Tree 13 fixed INDELs']] <- length(geno1_13_spec_hets)

# Based on these preliminary results, I don't think there is the sensitivity
#  to detect branch-specific deletions.

# INFO ABOUT QUALITY OF SVs

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

stat_list[['Number of SVs variable in both trees']] <- length(
  intersect(t13_var_inds, t14_var_inds))
# 16
stat_list[['Number of SVs variable in Tree13']] <- length(t13_var_inds)
# 112
stat_list[['Number of SVs variable in Tree14']] <- length(t14_var_inds)
# 64
stat_list[['Number of singleton variable SVs']] <- sum(apply(
  del1_singletons[,c(2:4)], 2, sum)) + sum(apply(
  ins1_singletons[,c(2:4)], 2, sum))
# 124

# Look for loss-of-heterozygosity - find where homozygosity is found in
#  branch(es) above branches that are heterozygous
# Possible approach: set branch names in correct order,
#  get indices of branches with het and hom genotypes
#  subtract max index of each - ex max(het) - max(hom); negative number means
#  loss of heterozygosity
# Subtract highest homozygous branch number from lowest heterozygous branch
#  number; if number is negative, then there is a loss-of-heterozygosity on
#  that branch
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

stat_list[['N variable SVs with loss-of-heterozygosity']] <- length(
  lossOfHet_tot)
stat_list[['N variable SVs with gain-of-heterozygosity']] <- nrow(
  geno1_var) - length(lossOfHet_tot)

###
# find SVs that are fixed for 1 in one clone and variable in the other clone;
#  these are suspicious because should not be getting heterozygosity arising
#  twice independently
t13_fixed_het <- which(apply(geno1_var[, t13_cols], 1, 
  function(x) sum(x == 1)) == length(t13_cols))

t13fixHet_t14Var <- intersect(t13_fixed_het, t14_var_inds)

t14_fixed_het <- which(apply(geno1_var[, t14_cols], 1, 
  function(x) sum(x == 1)) == length(t14_cols))

t14fixHet_t13Var <- intersect(t14_fixed_het, t13_var_inds)

stat_list[['N SVs fixedHet in 1 tree and var in other']] <- length(
  union(t13fixHet_t14Var, t14_fixed_het))

# Final tally of good SVs
bad_var_SVs <- union(lossOfHet_tot, 
  union(union(t13fixHet_t14Var, t14fixHet_t13Var),
  intersect(t13_var_inds, t14_var_inds)))
decent_var_SVs <- setdiff(seq(nrow(geno1_var)), bad_var_SVs)

singleton_inds <- which(
  apply(geno1_var, 1, function(x) sum(table(x) == 1)) == 1)

stat_list[['N gain-of-het SVs variable in 1 tree']] <- length(decent_var_SVs)
stat_list[['N gain-of-het SVs var in 1 tree and var in 2+ branches']] <- (
 length(setdiff(decent_var_SVs, singleton_inds)) )

# geno1_var[decent_var_SVs, c(t13_cols, t14_cols)]

stat_list[['N gain-of-het DELs var in 1 tree']] <- length(intersect(
  decent_var_SVs, which(geno1_var_type == 'DEL')))
stat_list[['N gain-of-het DELs > 100bp var in 1 tree']] <- length(intersect(
  intersect(decent_var_SVs, which(geno1_var_type == 'DEL')), 
  which(geno1_var_size >= 100)))

stat_list[['N gain-of-het INSs var in 1 tree']] <- length(intersect(
  decent_var_SVs, which(geno1_var_type == 'INS')))
stat_list[['N gain-of-het INSs > 100bp var in 1 tree']] <- length(intersect(
  intersect(decent_var_SVs, which(geno1_var_type == 'INS')),
  which(geno1_var_size >= 100)))

# save genotypes if there are any "decent" SVs
if(length(decent_var_SVs) > 0){
  saveRDS(data.frame(geno1_var)[decent_var_SVs, ], file = good_var_genos_out)
}

# Write final files
stat_df <- data.frame(label = names(stat_list), value = unlist(stat_list),
  stringsAsFactors = F)

write.table(stat_df, file = stat_list_out, quote = F, sep = '\t', 
  row.names = F, col.names = T)

write.table(samp_info_df, file = samp_info_out, quote = F, sep = '\t', 
  row.names = F, col.names = T)

write.table(samp_ord_df, file = samp_ord_corr_out, quote = F, sep = '\t',
  row.names = F, col.names = T)

write.table(samp_sing_df, file = samp_sing_out, quote = F, sep = '\t', 
  row.names = F, col.names = T)

write.table(sing_ord_df, file = sing_ord_corr_out, quote = F, sep = '\t',
  row.names = F, col.names = T)

saveRDS(filt_indel_names, file = filt_indel_names_out)

quit(save = 'no')

