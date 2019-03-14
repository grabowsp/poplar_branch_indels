# Analysis of the joint/combined genotype calling from pbsv v2.1.1

# LOAD LIBRARIES AND PACKAGES
function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)

# LOAD DATA
data_dir <- '/home/f1p1/tmp/PBSV/Poplar14.5/LIBS/'
combo_1_vcf_short <- 'Ptr145v1.ALLData.v2.try1.vcf'
combo_1_vcf_file <- paste(data_dir, combo_1_vcf_short, sep = '')

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')
# SET OUTPUT
## General Stat list
stat_list <- list()

# SET VARIABLES
branch_13_lab <- c(13.1, 13.2, 13.3, 13.4, 13.5)
branch_14_lab <- c(14.1, 14.2, 14.3, 14.4, 14.5)

# SET CONSTANTS


############
# Raw Counts
raw_vcf <- read.table(combo_1_vcf_file, header = F, stringsAsFactors = F,
  sep = '\t')

raw_info <- strsplit(raw_vcf[,8], split = ';')
type_info <- unlist(lapply(raw_info, 
  function(x) unlist(strsplit(x[[1]], split = '='))[2]))

raw_type_tab <- table(type_info)
#   BND   DEL   INS   INV 
# 28432 29089 21618    31

stat_list[['raw_type_count']][['label']] <- paste('raw', names(raw_type_tab), 
  sep = ' ')
stat_list[['raw_type_count']][['value']] <- as.vector(raw_type_tab)

bnd_ins_equiv <- raw_type_tab[names(raw_type_tab) == 'BND'] / 4

stat_list[['BND_INS_equiv']][['label']] <- 'BND count / 4 (possible INS)'
stat_list[['BND_INS_equiv']][['value']] <- bnd_ins_equiv
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

nrow(combo1_genos_noNA)
# [1] 12646

stat_list[['tot_filt_INDELs']][['label']] <- 'Total number of filtered INDELs'
stat_list[['tot_filt_INDELs']][['value']] <- nrow(combo1_genos_noNA)

# CONTINUE GENERATING STANDARD OBJECTS FOR SAVING INFO FROM THIS POINT
scaf_inds <- grep('scaffold', rownames(combo1_genos_noNA))
length(scaf_inds)
# [1] 3

combo1_genos_2 <- combo1_genos_noNA[-scaf_inds, ]
combo1_genos_info <- strsplit(rownames(combo1_genos_2), split = '_')
combo1_genos_type <- unlist(lapply(combo1_genos_info, function(x) x[[3]]))

table(combo1_genos_type)
#  DEL  INS 
# 7647 4996 

combo1_g2_del_inds <- grep('DEL', combo1_genos_type)
combo1_g2_ins_inds <- grep('INS', combo1_genos_type)

combo1_genos_chr <- unlist(lapply(combo1_genos_info, function(x) x[[1]]))

table(combo1_genos_chr[combo1_g2_del_inds])
# Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 
#   867   529   417   522   532   531   334   366   254   462   399   313   380 
# Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 
#   337   235   308   308   290   263

table(combo1_genos_chr[combo1_g2_ins_inds])
# Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 
#   593   358   300   328   367   375   209   227   213   321   204   162   209 
# Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 
#   225   165   208   188   213   131

combo1_genos_size <- unlist(lapply(combo1_genos_info, 
  function(x) as.numeric(x[[4]])))
combo1_g2_100bp_inds <- which(combo1_genos_size >= 100)

length(combo1_g2_100bp_inds)

length(intersect(combo1_g2_100bp_inds, combo1_g2_del_inds))
# [1] 2562

length(intersect(combo1_g2_100bp_inds, combo1_g2_ins_inds))
# [1] 2593

table(combo1_genos_chr[intersect(combo1_g2_100bp_inds, combo1_g2_del_inds)])
# Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 
#   296   191   136   177   175   192   116   128    95   155   128   101   114 
# Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 
#   110    75   104   101    90    78 

table(combo1_genos_chr[intersect(combo1_g2_100bp_inds, combo1_g2_ins_inds)])
# Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 
#   303   186   161   168   191   175   109   133   121   176   107    74   113 
# Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 
#   103    86    96   100   118    73

# Variable INDEL analysis

n_geno_vec <- apply(combo1_genos_2, 1, function(x) length(table(x)))

geno1_var <- combo1_genos_2[which(n_geno_vec > 1), ]

nrow(geno1_var)
# 160

geno1_var_info <- strsplit(rownames(geno1_var), split = '_')
geno1_var_type <- unlist(lapply(geno1_var_info, function(x) x[[3]]))

sum(geno1_var_type == 'DEL')
# 77
sum(geno1_var_type == 'INS')
# 83

geno1_var_size <- unlist(lapply(geno1_var_info, function(x) as.numeric(x[[4]])))

length(intersect(which(geno1_var_type == 'DEL'), which(geno1_var_size >= 100)))
# 62

length(intersect(which(geno1_var_type == 'INS'), which(geno1_var_size >= 100)))
# 28

apply(geno1_var, 2, function(x) sum(x == 1))
# PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW 
#   92  105  101   96   99  102  108  108  119  136

del1_var <- geno1_var[which(geno1_var_type == 'DEL'), ]
apply(del1_var, 2, function(x) sum(x == 1))
# PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW 
#   65   72   70   65   65   68   69   64   64   67

del1_var_100 <- geno1_var[intersect(which(geno1_var_type == 'DEL'), 
  which(geno1_var_size >= 100)), ]
apply(del1_var_100, 2, function(x) sum(x == 1))
# PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW 
#   51   57   56   52   52   55   55   49   51   56

ins1_var <- geno1_var[which(geno1_var_type == 'INS'), ]
apply(ins1_var, 2, function(x) sum(x == 1))
# PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW 
#   27   33   31   31   34   34   39   44   55   69

ins1_var_100 <- geno1_var[intersect(which(geno1_var_type == 'INS'),
  which(geno1_var_size >= 100)), ]
apply(ins1_var_100, 2, function(x) sum(x == 1))
# PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW 
#   11   18   17   14   18   15   17   18   20   21

### Singleton Analysis
basic_singleton_df <- data.frame(lib = colnames(del1_var), homRef = 0, het = 0,
 homAlt = 0, stringsAsFactors = F)

del1_singletons <- gen_singleton_df(geno_mat = del1_var,
  basic_df = basic_singleton_df)

apply(del1_singletons[,c(2:4)], 2, sum)
# homRef    het homAlt 
#     30      0     34

sum(apply(del1_singletons[,c(2:4)], 2, sum))
# [1] 64
# 64 of the 77 variable deletions are singletons - troubling
sum(apply(del1_singletons[,c(2:4)], 2, sum))/nrow(del1_var)
# [1] 0.8311688
 
del1_100_singletons <- gen_singleton_df(geno_mat = del1_var_100,
  basic_df = basic_singleton_df)

apply(del1_100_singletons[,c(2:4)], 2, sum)
# homRef    het homAlt 
#     21      0     28

sum(apply(del1_100_singletons[,c(2:4)], 2, sum))
# [1] 49
sum(apply(del1_100_singletons[,c(2:4)], 2, sum))/nrow(del1_var_100)
# [1] 0.7903226

ins1_singletons <- gen_singleton_df(geno_mat = ins1_var,
  basic_df = basic_singleton_df)

apply(ins1_singletons[,c(2:4)], 2, sum)
# homRef    het homAlt 
#    21     27     12

sum(apply(ins1_singletons[,c(2:4)], 2, sum))
# [1] 60
# 60 of 83 variable insertions are singletons
sum(apply(ins1_singletons[,c(2:4)], 2, sum))/nrow(ins1_var)
# [1] 0.7228916

ins1_100_singletons <- gen_singleton_df(geno_mat = ins1_var_100,
  basic_df = basic_singleton_df)

apply(ins1_100_singletons[,c(2:4)], 2, sum)
# homRef    het homAlt 
#    13     8     4

sum(apply(ins1_100_singletons[,c(2:4)], 2, sum))
# [1] 25
sum(apply(ins1_100_singletons[,c(2:4)], 2, sum))/nrow(ins1_var_100)
# 0.8928571

sing_list <- list()
sing_list[['DEL_20bp']] <- del1_singletons
sing_list[['DEL_100bp']] <- del1_100_singletons
sing_list[['INS_20bp']] <- ins1_singletons
sing_list[['INS_100bp']] <- ins1_100_singletons

samp_sing_df <- data.frame(lib = del1_singletons$lib, stringsAsFactors = F)

for(sl in seq(length(sing_list))){
  type_name <- names(sing_list)[sl]
  tmp_colnames <- paste(type_name, colnames(sing_list[[sl]])[2:4], sep = '_')
  samp_sing_df[, tmp_colnames] <- sing_list[[sl]][, c(2:4)]
}


# NEXT STEPS:
# 2) Generate Dataframe for counts of singletons by type for each sv type
#      and size
# 3) CLONE SVs
# 4) Error estimates



geno1_14_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno1_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_14_lab, meta = samp_meta)
# integer(0) - none

geno1_13_spec_hets <- get_shared_unique_het_inds(
  geno_mat = geno1_var,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_13_lab, meta = samp_meta)
# integer(0) - none
# I'm a bit worried about these results...

# Based on these preliminary results, I don't think there is the sensitivity
#  to detect branch-specific deletions.


quit(save = 'no')

