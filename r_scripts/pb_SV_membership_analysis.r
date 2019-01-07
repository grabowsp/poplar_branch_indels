# Analysis for looking at private, shared, etc SVs in the PacBio SV data

# LOAD LIBRARIES #
function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)

function_test_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.test.r'
# source(function_test_file)

# LOAD DATA #
data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'
combo_file <- 'ref.ALLData.vcf'
combo_file_tot <- paste(data_dir, combo_file, sep = '')
combo_df <- make_combo_indel_df(combo_file_tot)

meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')
# samp_meta[,c('lib_name', 'branch_name')]

####################

# Filter samples, SVs, and genotypes
## Remove bad branches
bad_branches <- c('13.4', '14.1')
bad_samps <- c()
for(bb in bad_branches){
  tmp_ind <- which(as.character(samp_meta$branch_name) == bb)
  bad_samps <- c(bad_samps, samp_meta$lib_name[tmp_ind])
}

bad_vcf_cols <- c()
for(bs in bad_samps){
  tmp_col_ind <- which(colnames(combo_df) == bs)
  bad_vcf_cols <- c(bad_vcf_cols, tmp_col_ind)
}

combo_df_filt <- combo_df[ , -sort(bad_vcf_cols)]

# Filtering example statistics
# how many positions are removed by the <make_combo_indel_df> function that
#   selects only INS and DEL
# I got these values by manipulating some of the code in <make_combo_indel_df>

# table(combo_vcf$type)
#   BND   DEL   INS   INV 
#  2154 29852 26942    12


# REF is one letter for INSertions; ALT is one letter for DELetions

# indices to remove that are either 70% mononucleotide repeats, 50% the same
#   binucleotide repeat, or 60% of two binucleotide repeats
remove_8mer_SV_inds <- id_prob_SV_seqs(sv_geno_df = combo_df_filt, 
  mer_length = 8, per_mn_cutoff = 0.7, per_bn_pure_cutoff = 0.5, 
  per_bn_multi_cutoff = 0.6)

length(remove_8mer_SV_inds)
# 6073
length(remove_8mer_SV_inds)/nrow(combo_df_filt)
# 10.7%

combo_filt_1 <- combo_df_filt[-remove_8mer_SV_inds, ]

# remove duplicated positions - the smaller SVs are removed
rm_dup_inds <- dup_inds_to_remove(sv_geno_df = combo_filt_1)

length(rm_dup_inds)
# 185
length(rm_dup_inds)/nrow(combo_filt_1)
# 0.003648196

combo_filt_2 <- combo_filt_1[-rm_dup_inds, ]

# remove SV's that are within 30bp of eachother - keep the larger SV
too_close_inds_30 <- overlap_inds_to_remove(sv_geno_df = combo_filt_2,
  dist_cut = 30) 

length(too_close_inds_30)
# 543
length(too_close_inds_30)/nrow(combo_filt_2)
#  0.01074715

combo_filt_3 <- combo_filt_2[-too_close_inds_30, ]

## find and remove any SVs with missing genotypes in any sample
combo_filt_4 <- remove_missing_combo_SVs(combo_filt_3)

nrow(combo_filt_3) - nrow(combo_filt_4)
# 31
(nrow(combo_filt_3) - nrow(combo_filt_4))/nrow(combo_filt_3)
# 0.0006202233

# find and remove SVs with excess coverage
excess_cov_inds <- id_excess_coverage(combo_filt_4, sd_cut = 3)

length(excess_cov_inds)
# 26
length(excess_cov_inds)/nrow(combo_filt_4)
# 0.0005205101

combo_df_2 <- combo_filt_4[-excess_cov_inds, ]

combo_geno_info_list <- make_allsamp_geno_info_list(combo_df_2)

combo_filtered_genotypes <- call_allsamp_genotypes(combo_geno_info_list,
  min_sv_coverage = 10, het_ratio_cut = 0.15, min_minor_allele_count = 2,
  max_hom_ratio = 0.05)

sum(is.na(combo_filtered_genotypes))
# 36541

sum(is.na(combo_filtered_genotypes)) / length(combo_filtered_genotypes)
#  0.09148973

geno_SV_NAs <- apply(combo_filtered_genotypes, 1, function(x) sum(is.na(x)))
sum(geno_SV_NAs > 0)
# 16105
sum(geno_SV_NAs > 0)/length(geno_SV_NAs)
# 0.3225839

# assign scores for genotypes based on penetrance/allele count ratios
#   and presnece of other sampes with same, high-score genotypes
combo_filt_adj_scores <- get_adjusted_geno_scores(
  geno_info_list = combo_geno_info_list, 
  genotype_mat = combo_filtered_genotypes, 
  est_error = 0.032, est_penetrance = 0.35, good_score = 0.9, 
  great_score = 0.98, same_geno_bonus = 0.67)

# find SVs that do not have at least one sample above a coverage threshold
#   ex: minimum coverage may be 10, but want at least one sample with coverage
#   above a higher number, like 20
low_ss_svs <- find_SVs_below_ss_cov(combo_geno_info_list, ss_min_cov = 20)

combo_filt_genos_2 <- combo_filtered_genotypes
combo_filt_genos_2[which(combo_filt_adj_scores < 0.9)] <- NA

sum(is.na(combo_filt_genos_2)) - sum(is.na(combo_filtered_genotypes))
# 2630

((sum(is.na(combo_filt_genos_2)) - sum(is.na(combo_filtered_genotypes)))/
  length(combo_filt_genos_2))
# 0.006584877

geno2_SV_NAs <- apply(combo_filt_genos_2, 1, function(x) sum(is.na(x)))
length(setdiff(which(geno2_SV_NAs > 0), which(geno_SV_NAs > 0)))
# 672

(length(setdiff(which(geno2_SV_NAs > 0), which(geno_SV_NAs > 0))) / 
  length(geno2_SV_NAs))
# 0.01346019

length(setdiff(low_ss_svs, which(geno2_SV_NAs > 0)))
# 26

length(setdiff(low_ss_svs, which(geno2_SV_NAs > 0))) / length(geno2_SV_NAs)
# 0.0005207812

combo_filt_genos_2[low_ss_svs, ] <- NA

combo_filt_genos_test_noNAs <- filt_geno_mat(combo_filt_genos_2, max_nas = 0)
nrow(combo_filt_genos_test_noNAs)
# [1] 33122

orig_n_SVs <- sum(c(2154, 29852, 26942, 12))
orig_n_SVs - nrow(combo_filt_genos_test_noNAs)
# 25838
(orig_n_SVs - nrow(combo_filt_genos_test_noNAs))/orig_n_SVs
# 0.4382293

combo_filt_genos_noNAs <- filt_geno_mat(combo_filtered_genotypes, max_nas = 0)
nrow(combo_filt_genos_noNAs)
# [1] 33820

filt_genos_test_noNAs_names 

df_in_genos_inds <- which(combo_df_2$full_name_long %in% rownames(
  combo_filt_genos_test_noNAs))

combo_filt_genos_noNAs_df <- combo_df_2[df_in_genos_inds, ]

length_vec <- c(15, 50, 100, 500, 1000, 10000)
sapply(length_vec, function(x) sum(combo_filt_genos_noNAs_df$sv_length > x))


sum(combo_df_2$full_name_long[df_in_genos_inds] != rownames(combo_filt_genos_test_noNAs))
# 0

##########

# Extract genotype information. Includes
## converting genos to numeric
## getting coverage
## getting number of reads supporting each allele
combo_geno_info_list <- make_allsamp_geno_info_list(combo_df_2)

# generate matrix of coverage
geno_cov_mat <- matrix(
  data = unlist(lapply(combo_geno_info_list, function(x) x[[2]])),
  ncol = length(combo_geno_info_list),
  byrow = F)

min_cov_vec <- apply(geno_cov_mat, 1, min)

min_cov_floor_vec <- c()
for(mc in seq(5, 30, 5)){
  tmp_val <- sum(min_cov_vec >= mc)
  min_cov_floor_vec <- c(min_cov_floor_vec , tmp_val)
}
names(min_cov_floor_vec) <- seq(5, 30, 5)

min_cov_floor_vec
#     5    10    15    20    25    30 
# 49830 47436 34821  9707  1958   900 

# only 900 SV's have a minimum coverage of 30; 47436 have min coverage of 10

# RULES I'M GOING BY AS OF 11/2 - these can change, but I'm going to use these
#  for filtering
# 1) Minimum coverage of 10 for every sample
# 2) Het calls must have at least 2 reads of each allele
# 3) Het calls must have at least a 0.15 allele ratio - may want to be more
#      conservative with this - maybe 0.2 or 0.25
# 4) Homozygous call must have at lower than 0.05 allele ratio - this may
#      need to be adjusted

# Calling my own genotypes

combo_filtered_genotypes <- call_allsamp_genotypes(combo_geno_info_list, 
  min_sv_coverage = 10, het_ratio_cut = 0.15, min_minor_allele_count = 2, 
  max_hom_ratio = 0.05)

num_nas <- apply(combo_filtered_genotypes, 1, function(x) sum(is.na(x)))
table(num_nas)
#     0     1     2     3     4     5     6     7     8 
# 33820  6921  3636  2342  1581   976   458   158    33
## 33.8k SVs have no missing data

######
# Estimate (biased) error rates using homozygous genotype calls
filt_error_obj <- estimate_error_objs(
  genotype_mat = combo_filtered_genotypes, 
  geno_info_list = combo_geno_info_list)

# see what happens if allow for higher higher allele ratios in homozygotes
alt_filtered_genotypes <- call_allsamp_genotypes(combo_geno_info_list,
  min_sv_coverage = 10, het_ratio_cut = 0.15, min_minor_allele_count = 2,
  max_hom_ratio = 0.1)

alt_error_obj <- estimate_error_objs(
  genotype_mat = alt_filtered_genotypes,
  geno_info_list = combo_geno_info_list)

# mean error rates based on the number of homozyotes at a SV
round(filt_error_obj[[2]],5)
#      0       1       2       3       4       5       6       7       8 
#     NA      NA 0.01655 0.01481 0.01202 0.01000 0.00761 0.00539 0.00554

round(alt_error_obj[[2]],5)
#      0       1       2       3       4       5       6       7       8 
#     NA      NA 0.04932 0.04151 0.03349 0.02862 0.02370 0.01638 0.02055

table(filt_error_obj[[1]][,2])
#     0     1     2     3     4     5     6     7     8 
# 39058  3271  2131  1716  1452  1117   680   476    24 

table(alt_error_obj[[1]][,2])
#     0     1     2     3     4     5     6     7     8 
# 35696  3421  2389  2162  2211  1813  1208   937    88

# Toughts:
# 1) When there are fewer homozygous genotypes at a SV, the allele ratio
#     gets higher and higher, meaning that the singletons are probably
#     fairly suspicious in that many have a high allele ratio
# 2) If we use a higher allele ratio cutoff for homozygotes, then,
#      predictably, the error rate for homozygotes goes up
# 3) Related to (2), I think that our estimate for error rates can be tied
#      to the cutoff we use for calling homozygotes 
# 4) To assigne score/probability to homozygous calls, can choose
#      an error rate of mean(x) + 1*sd(x) where x is the error rates when
#      there are 2 or more homozygous genotypes.

round(filt_error_obj[[2]] + (1*filt_error_obj[[3]]),5)
#      0       1       2       3       4       5       6       7       8 
#     NA      NA 0.03181 0.02727 0.02210 0.01904 0.01552 0.01188 0.01204

round(alt_error_obj[[2]] + (1*alt_error_obj[[3]]),5)
#      0       1       2       3       4       5       6       7       8 
#     NA      NA 0.07578 0.06515 0.05545 0.04941 0.04238 0.03225 0.03570
 
# assign scores for genotypes based on penetrance/allele count ratios
#   and presnece of other sampes with same, high-score genotypes
combo_filt_adj_scores <- get_adjusted_geno_scores(
  geno_info_list = combo_geno_info_list, 
  genotype_mat = combo_filtered_genotypes, 
  est_error = 0.032, est_penetrance = 0.35, good_score = 0.9, 
  great_score = 0.98, same_geno_bonus = 0.67)

# find SVs that do not have at least one sample above a coverage threshold
#   ex: minimum coverage may be 10, but want at least one sample with coverage
#   above a higher number, like 20
low_ss_svs <- find_SVs_below_ss_cov(combo_geno_info_list, ss_min_cov = 20)

# continue from here...
combo_filt_genos_2 <- combo_filtered_genotypes
combo_filt_genos_2[which(combo_filt_adj_scores < 0.9)] <- NA

combo_filt_genos_2[low_ss_svs, ] <- NA

combo_filt_genos_test_noNAs <- filt_geno_mat(combo_filt_genos_2, max_nas = 0)
nrow(combo_filt_genos_test_noNAs)
# [1] 33122

combo_filt_genos_noNAs <- filt_geno_mat(combo_filtered_genotypes, max_nas = 0)
nrow(combo_filt_genos_noNAs)
# [1] 33820

# filtering based on genotype scores and minimum single-sample coverage
#   removes 698 (2%) of SVs - worth doing to be conservative



#####################

# SIDE TASK: COVERAGE VS PENETRANCE (allele ratio)
samp_1_cov <- combo_geno_info_list[[1]][[2]]
samp_1_al_df <- combo_geno_info_list[[1]][[3]]
samp_1_alt_count <- samp_1_al_df[,2]
samp_1_min_count <- apply(samp_1_al_df, 1, min)

seq_cov <- lapply(combo_geno_info_list, function(x) x[[2]])
seq_alt_count <- lapply(combo_geno_info_list, function(x) x[[3]][,2])
seq_min_count <- lapply(combo_geno_info_list, function(x) apply(x[[3]], 1, min))

combo_filt_cov_v_alt <- coverage_vs_penetrance(combo_geno_info_list, 
  use_alt = T)
combo_filt_cov_v_min_0 <- coverage_vs_penetrance(combo_geno_info_list, 
  use_alt = F)

# want to remove sites with suspiciously high coverage and homozygous
#   genotype calls because only interested in allele ratio in heterozygotes
too_hi_cov <- which(combo_filt_cov_v_min_0$coverage > 400)

combo_filt_cov_v_min_filt <- combo_filt_cov_v_min_0[-too_hi_cov, ]

min_per_vec <- c(0.01, 0.05, 0.1, 0.2)
library(ggplot2)
pen_list <- list()
for(mp in seq(length(min_per_vec))){
  tmp_min_per <- min_per_vec[mp]
  tmp_no_var_inds <- which(combo_filt_cov_v_min_filt$per_min < tmp_min_per)
  tmp_filt_df <- combo_filt_cov_v_min_filt[-tmp_no_var_inds, ]
  tmp_cov_v_pen <- tapply(tmp_filt_df$per_min, tmp_filt_df$coverage, mean)
  tmp_cov_v_pen_df <- data.frame(
    coverage = as.numeric(names(tmp_cov_v_pen)), mean_per_min = tmp_cov_v_pen, 
    stringsAsFactors = F)
  tmp_g_cov_mean_min <- ggplot(tmp_cov_v_pen_df, aes(x = coverage, 
    y = mean_per_min)) +
    geom_point() +
    labs(title = paste('Mean penetrance v coverage;\nmin het ratio = ', 
      tmp_min_per, sep = ''))
  pen_list[[mp]] <- tmp_g_cov_mean_min
}

combo_filt_cov_v_pen_png <- paste('/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/analysis_results/figs/',
  'filtered_SVs_coverage_v_penetrance.png', sep = '')

library(gridExtra)

png(combo_filt_cov_v_pen_png, width = 600, height = 600)
#g_combo_filt_cov_v_min
#for(i in seq(length(pen_list))){
#  do.call(grid.arrange, c(pen_list[[i]], ncol = 2))
#}
do.call(grid.arrange, c(pen_list, ncol = 2))
dev.off()

###############

# CONTINUE FROM HERE
# need to figure out a way to id het postions and look at penetrance to make
#  sure is a good position 


combo_filt_genos_noNAs <- filt_geno_mat(combo_filtered_genotypes, max_nas = 0)

##########
# Streamline filtering and genotype calling

precall_filt_df <- precalling_SV_filtering(combo_df_filt)
# sampe as combo_df_2

combo_genotypes_1 <- generate_filtered_genotypes(combo_df = precall_filt_df,
  min_sv_coverage = 10, het_ratio_cut = 0.15, min_minor_allele_count = 2, 
  max_hom_ratio = 0.05, est_error = 0.032, est_penetrance = 0.35,
  good_score = 0.9, great_score = 0.98, same_geno_bonus = 0.67, 
  ss_min_cov = 20)

combo_genos_noNA_50 <- filt_geno_mat(geno_mat = combo_genotypes_1, 
  max_nas = 0, min_length = 50)

###############
###############




tally_singleton_svs(combo_filt_genos_noNAs)
# $homRef_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#  62   25   75   82   74   47   53  247 
# $het_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#  413   26   31   40   27   14    9   26 
# $homAlt_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#    6   13    4    3    4    9   16    4 

# note: allowing 1 NA generate exact same tally (more SV in genoype matrix
#  but same numbers of singletons

genos_noNAs_50 <- filt_geno_mat(combo_filtered_genotypes, max_nas = 0, 
  min_length = 50)
tally_singleton_svs(genos_noNAs_50)
# $homRef_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#   30    5   39   26   31   10    7   81 
# $het_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#  46    6    1    8    1    1    2    4 
# $homAlt_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#    5   11    4    2    4    9   15    3 

# note: allowing 1 NA generates exact same tally

genos_noNAs_100 <- filt_geno_mat(combo_filtered_genotypes, max_nas = 0, 
  min_length = 100)
tally_singleton_svs(genos_noNAs_100)
# $homRef_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#   19    2   30   15   21    1    3   35 
# $het_singletons
# PAXL PAXN PAYK PAZF PBAT PBAU PBAW 
#   12    3    1    4    1    1    1 
# $homAlt_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#    5   10    4    2    4    9   15    3 

genos_noNAs_1000 <- filt_geno_mat(combo_filtered_genotypes, max_nas = 0,
  min_length = 1000)
tally_singleton_svs(genos_noNAs_1000)
# $homRef_singletons
# PAXL PAXN PAYK PAZF PAZH PBAU PBAW 
#    6    2   11    3   13    1   15 
# $het_singletons
# PAXL 
#    2 
# $homAlt_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#    3    6    3    2    4    8   14    3

# Look at DELetions vs INSertions
del_genos_noNAs <- filt_geno_mat(combo_filtered_genotypes, max_nas = 0,
  sv_type = 'DEL')
tally_singleton_svs(del_genos_noNAs)
# $homRef_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#   18    7   26   13   17    4    1   37 
# $het_singletons
# PAXL PAYK PAZF PAZH PBAT PBAW 
#    3    1    1    1    1    3 
# $homAlt_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#    5   13    4    3    4    9   15    4

ins_genos_noNAs <- filt_geno_mat(combo_filtered_genotypes, max_nas = 0,
  sv_type = 'INS')
tally_singleton_svs(ins_genos_noNAs)
# $homRef_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#   44   18   49   69   57   43   52  210 
# $het_singletons
# PAXL PAXN PAYK PAZF PAZH PBAT PBAU PBAW 
#  410   26   30   39   26   13    9   23 
# $homAlt_singletons
# PAXL PBAU 
#    1    1 

ins_noNAs_numHets <- apply(ins_genos_noNAs, 1, function(x) sum(x == 1))
table(apply(ins_genos_noNAs[which(ins_noNAs_numHets == 1), ], 1, sum))
#   1   5  11  13 
# 576   1   1   1
# most of the singleton het INSertions are new instance of that SV 


# Assemble into a table
filt_geno_list_for_singletons <- list(combo_filt_genos_noNAs, genos_noNAs_50, 
  genos_noNAs_100, genos_noNAs_1000, del_genos_noNAs, ins_genos_noNAs)

singleton_res_df <- data.frame(lib = colnames(combo_filtered_genotypes), 
  stringsAsFactors = F)
singleton_res_df$branch <- NA
for(sl in seq(nrow(singleton_res_df))){
  meta_ind <- which(samp_meta$lib_name == singleton_res_df$lib[sl])
  singleton_res_df$branch[sl] <- samp_meta$branch_name[meta_ind]
}

type_labs <- c('homRef', 'het', 'homAlt')
test_tab <- singleton_res_df

top_lab_vec <- c('InDel_20bp', 'InDel_50bp', 'InDel_100bp', 'InDel_1000bp',
  'DEL_20bp', 'INS_20bp')

for(k in seq(length(filt_geno_list_for_singletons))){
  # 
  tmp_sing_tal <- tally_singleton_svs(filt_geno_list_for_singletons[[k]])
  top_lab <- top_lab_vec[k]
  #
  for(j in seq(length(tmp_sing_tal))){
    tmp_lab <- paste(top_lab, type_labs[j], sep = '_')
    tmp_tal <- tmp_sing_tal[[j]]
    for(i in seq(length(tmp_tal))){
      tab_ind <- which(singleton_res_df$lib == names(tmp_tal[i]))
      test_tab[tab_ind, tmp_lab] <- tmp_tal[i]
    }
  }
}

singleton_res_df <- test_tab[order(test_tab$branch), ]

sing_tab_out <- paste(data_dir, 'analysis_results/singleton_SV_tables.tab', 
  sep = '')

write.table(singleton_res_df, file = sing_tab_out, quote = F, sep = '\t', 
  row.names = F, col.names = T)


# Some notes:
# PBAW, 13.2, always has more singleton homozygous Ref SVs, meaning the 
#  Alt allele is missing ONLY from that sample; PBAW has lower seq depth, but
#  very close to PAZF/13.5 and PAZF does not show the same pattern
# PBAU has more large homozygous Alt SVs, meaning the rest of the samples
#  are heterozygous
# PAXL(14.3) has noticably more singleton heterozygous genotypes
# Most singleton homozygous Ref and heterozygous SVs are INSertions; most
#  singleton homozygous Alt SVs are DELetions. From a technical standpoint,
#  the singleton homozygous refs are missing the allele so could be that they
#  didn't have enough coverage to call the SV - should check if that changes
#  with higher coverage cutoff. Singleton hets are mainly a new/unique
#  occurance of that INSertion.
# I'm having trouble making sense of these patterns... maybe looking at other
#   patterns will help

branch_13_lab <- c(13.1, 13.2, 13.3, 13.5)
branch_14_lab <- c(14.2, 14.3, 14.4, 14.5)

b13_all_hets <- get_samegeno_inds(geno_mat = combo_filt_genos_noNAs, 
  branch_name_vec = branch_13_lab, genotype = 1, meta = samp_meta)
b13_all_homRef <- get_samegeno_inds(geno_mat = combo_filt_genos_noNAs, 
  branch_name_vec = branch_13_lab, genotype = 0, meta = samp_meta)
b13_all_homAlt <- get_samegeno_inds(geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = branch_13_lab, genotype = 2, meta = samp_meta)

b14_all_hets <- get_samegeno_inds(geno_mat = combo_filt_genos_noNAs, 
  branch_name_vec = branch_14_lab, genotype = 1, meta = samp_meta)
b14_all_homRef <- get_samegeno_inds(geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = branch_14_lab, genotype = 0, meta = samp_meta)
b14_all_homAlt <- get_samegeno_inds(geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = branch_14_lab, genotype = 2, meta = samp_meta)

# Start with clone 14 specific SVs
b14_specific_hets <- get_shared_unique_het_inds(
  geno_mat = combo_filt_genos_noNAs, 
  branch_name_vec = c(branch_13_lab, branch_14_lab), 
  test_names = branch_14_lab, meta = samp_meta)
# [1] 32
# 32 het SVs found ONLY in clone 14

# Out of curiosity, how many are het in 14 and only 1 sample in 13?
b13_1_het <- get_samegeno_inds(geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = branch_13_lab, genotype = 1, meta = samp_meta, n_miss = 3)
length(intersect(b14_all_hets, b13_1_het))
# [1] 274
# So it looks like there are some b13 branches that are more similar to b14?

b14_135_spec_hets <- get_shared_unique_het_inds(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = c(branch_14_lab, branch_13_lab[4]), meta = samp_meta)
length(b14_135_spec_hets)
# [1] 245
# most SVs that are in clone 14 and one clone 13 branch are in 13.5

test_4het_allsamps <- get_samegeno_inds(geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab), genotype = 1, meta = samp_meta, n_miss = 4)
# 649
# 8 choose 4 is 70, 659/70 = 9.4; so should be only 9 or 10 for a set of 
# 4 samples if totally random; instead have 32, so that's promising...

# check all the possible combinations of 4 samples to see how many shared,
#   unique het SVs are in each combination
n_het_in_4_samps <- test_shared_unique_het_combos(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab), n_test = 4, 
  meta = samp_meta)

n_het_in_4_samps[[1]]
#  [1]   2   0   3   0  38   1   2   0   3   6   1   0   0   0   0   1   1   1   9
# [20]   0   1   0   0   3   0   2   0   0   1   0   0  12   0   0   0   0   0   0
# [39]   0   0   0   0   0   0   0   1   1   0   1   2   0  25   1   0   0   2   0
# [58]   0   1   2   3  12   9   1   0 449   7   8   3  32

# SURPRISE!!!! It turns out, combination 66, which is 13.5, 14.1, 14.2, and 
#   14.3 has many more unique SVs than either clone-specific combo. Second
#   highest... combination 5 which is 13.1, 13.2, 13.3, 14.5. Third highest
#   is clone 14 branches (combination 70) and clone 13 branches are below many
#   other combinations. WHY?????

table(n_het_in_4_samps[[1]])
#  0   1   2   3   6   7   8   9  12  25  32  38 449 
# 34  14   6   5   1   1   1   2   2   1   1   1   1 

b13_specific_hets <- get_shared_unique_het_inds(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_13_lab, meta = samp_meta)

length(b13_specific_hets)
# [1] 2
# 2 heterozygous SVs found ONLY in clone 13

c(branch_13_lab, branch_14_lab)[
  n_het_in_4_samps[[2]][, which(n_het_in_4_samps[[1]] == 449)]]
# [1] 13.5 14.2 14.3 14.4

c(branch_13_lab, branch_14_lab)[
  n_het_in_4_samps[[2]][, which(n_het_in_4_samps[[1]] == 38)]]
# [1] 13.1 13.2 13.3 14.5

c(branch_13_lab, branch_14_lab)[
  n_het_in_4_samps[[2]][, which(n_het_in_4_samps[[1]] == 25)]]
# [1] 13.2 14.2 14.3 14.4

####
# Look at combinations of 2 branches - the top two brances of each clone
#   should have the hightest number of shared, unique SVs

b13_top2_specific_hets <- get_shared_unique_het_inds(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_13_lab[1:2], meta = samp_meta)

length(b13_top2_specific_hets)
# [1] 29

b14_top2_specific_hets <- get_shared_unique_het_inds(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_14_lab[1:2], meta = samp_meta)
length(b14_top2_specific_hets)
# [1] 289

n_het_in_2_samps <- test_shared_unique_het_combos(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab), n_test = 2,
  meta = samp_meta)

table(n_het_in_2_samps[[1]])
#   0   1   2   3   4  13  15  17  18  19  22  29 289 
#   7   6   3   1   1   3   1   1   1   1   1   1   1
# 14.2 and 14.3 have the most by far; 13.1 and 13.2 have second most

c(branch_13_lab, branch_14_lab)[
  n_het_in_2_samps[[2]][, which(n_het_in_2_samps[[1]] == 22)]]
# [1] 13.5 14.4

c(branch_13_lab, branch_14_lab)[
  n_het_in_2_samps[[2]][, which(n_het_in_2_samps[[1]] == 19)]]
# [1] 13.3 14.5

c(branch_13_lab, branch_14_lab)[
  n_het_in_2_samps[[2]][, which(n_het_in_2_samps[[1]] == 18)]]
# [1] 13.5 14.3

c(branch_13_lab, branch_14_lab)[
  n_het_in_2_samps[[2]][, which(n_het_in_2_samps[[1]] == 17)]]
# [1] 14.2 14.4

c(branch_13_lab, branch_14_lab)[
  n_het_in_2_samps[[2]][, which(n_het_in_2_samps[[1]] == 15)]]
# [1] 14.3 14.5

c(branch_13_lab, branch_14_lab)[
  n_het_in_2_samps[[2]][, which(n_het_in_2_samps[[1]] == 13)]]
# [1] 13.2 14.3; 13.5 14.5; 14.3 14.4

###
# Look at top 3 branches

b13_top3_specific_hets <- get_shared_unique_het_inds(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab), 
  test_names = branch_13_lab[1:3], meta = samp_meta)
length(b13_top3_specific_hets)
# [1] 19

b14_top3_specific_hets <- get_shared_unique_het_inds(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab),
  test_names = branch_14_lab[1:3], meta = samp_meta)
length(b14_top3_specific_hets)
# [1] 506

n_het_in_3_samps <- test_shared_unique_het_combos(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab), n_test = 3,
  meta = samp_meta)

table(n_het_in_3_samps[[1]])
#   0   1   2   3   4   5   6   7   9  10  19  22  23 506 
#  24  12   3   2   1   1   5   1   1   2   1   1   1   1

c(branch_13_lab, branch_14_lab)[
  n_het_in_3_samps[[2]][, which(n_het_in_3_samps[[1]] == 23)]]
# [1] 13.5 14.2 14.4

c(branch_13_lab, branch_14_lab)[
  n_het_in_3_samps[[2]][, which(n_het_in_3_samps[[1]] == 22)]]
# [1] 13.5 14.2 14.3

# Top numbers other than 506 are combinations of clone 14 branches with 
#  branch 13.5

####
# look at SVs het in 5 samples, homozygous in 3 samples
# note: the function isn't exactly suited for testing shared, unique homozygous
#  genotypes, should make new function if want to pursue this farther

n_het_in_5_samps <- test_shared_unique_het_combos(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab), n_test = 5,
  meta = samp_meta)

table(n_het_in_5_samps[[1]])
#   0   1   2   3   5   8   9  14  16  24  27  43 245 
#  24  10   7   3   1   2   1   2   1   2   1   1   1

c(branch_13_lab, branch_14_lab)[
  n_het_in_5_samps[[2]][, which(n_het_in_5_samps[[1]] == 245)]]
# [1] 13.5 14.2 14.3 14.4 14.5
# 13.1, 13.2, and 13.3 have 245 unique, homozygous SVs
# would make sense if first split is between 13.5 and 13.3 OR if loss of
#  heterozygosity is more common than I assume

c(branch_13_lab, branch_14_lab)[
  n_het_in_5_samps[[2]][, which(n_het_in_5_samps[[1]] == 43)]]
# [1] 13.2 13.5 14.2 14.3 14.4

c(branch_13_lab, branch_14_lab)[
  n_het_in_5_samps[[2]][, which(n_het_in_5_samps[[1]] == 27)]]
# [1] 13.1 13.2 14.2 14.3 14.4

c(branch_13_lab, branch_14_lab)[
  n_het_in_5_samps[[2]][, which(n_het_in_5_samps[[1]] == 24)]]
# 13.1 13.2 13.3 13.5 14.5;   13.3 14.2 14.3 14.4 14.5

####
# Look at SVs in 6 samples
n_het_in_6_samps <- test_shared_unique_het_combos(
  geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = c(branch_13_lab, branch_14_lab), n_test = 6,
  meta = samp_meta)

table(n_het_in_6_samps[[1]])
#   0   1   2   3   4   5   6  10  11  16  18  20  29  35  48 203 
#   2   3   1   6   2   3   1   2   1   1   1   1   1   1   1   1

c(branch_13_lab, branch_14_lab)[
  n_het_in_6_samps[[2]][, which(n_het_in_6_samps[[1]] == 203)]]
# [1] 13.3 13.5 14.2 14.3 14.4 14.5
# 13.1 and 13.2 have 203 shared, unique homozygous SVs
# is this related to how 13.2 has an excess of unique homozygous SVs?

c(branch_13_lab, branch_14_lab)[
  n_het_in_6_samps[[2]][, which(n_het_in_6_samps[[1]] == 48)]]
# [1] 13.1 13.2 13.5 14.2 14.3 14.4
# shared homozygous SVs in 13.3 and 14.5 - doesn't really make sense




# What happens if we allow for miss-called 0's (actually 1's)
b13_all_hets_permissive <- get_samegeno_inds(geno_mat = combo_filt_genos_noNAs, 
  branch_name_vec = branch_13_lab, genotype = 1, meta = samp_meta, n_miss = 1)
b13_spec_hets_perm <- intersect(union(b13_all_hets, b13_all_hets_permissive), 
  b14_all_homRef)
length(b13_spec_hets_perm)
# [1] 24
# 22 more clone 13 specific het SVs if allow one sample to NOT be het. These
#   are either missing genotypes (called 0 but actually 1) or instances
#   of de novo mutation shared by 3 branches - should be 13.1, 13.2, 13.3
b13_libs <- get_branch_lib_names(branch_13_lab, meta = samp_meta)
b13_col_inds <- get_lib_col_inds(combo_filt_genos_noNAs, 
  lib_name_vec = b13_libs)
combo_filt_genos_noNAs[b13_spec_hets_perm, b13_col_inds]
# 2 are 1 in all 4
# 3 are 1 in all but 13.3 - maybe it's just missing from that one
# the rest (19) are 1 in all but 13.5, which is sort of what we'd expect

b14_all_hets_permissive <- get_samegeno_inds(geno_mat = combo_filt_genos_noNAs,
  branch_name_vec = branch_14_lab, genotype = 1, meta = samp_meta, n_miss = 1)
b14_spec_hets_perm <- intersect(union(b14_all_hets, b14_all_hets_permissive),
  b13_all_homRef)
length(b14_spec_hets_perm)
# 551
# A LOT!
b14_libs <- get_branch_lib_names(branch_14_lab, meta = samp_meta)
b14_col_inds <- get_lib_col_inds(combo_filt_genos_noNAs,
  lib_name_vec = b14_libs)
combo_filt_genos_noNAs[b14_spec_hets_perm[33:150], b14_col_inds]
# most of the "permissive" ones are actualy what we'd expect for real, denovo
#  SVs where they are missing in 14.5 but present in the other 3


combo_filt_genos_noNAs

branch_13_het <- 


# 1) Identify het SVs
# 2) Calculate % alt alleles in het calls
# 3) For each SV, aggregate the % alt alleles for each het call
# 4) Set desired min prob of homozygous genotype, ex: 0.001
# 5) For each SV, determine mean or min/max %alt alleles seen in het genotypes
# 5a)  depends on if is 0 (0/0) or 2 (1/1) genotype
# 6) Determine prob. of homozygous genotype given the prob of %alt alleles
#      in the het samples

het_inds <- which(geno_info_list[[1]][[1]] == 1)
test_ratio <- apply(geno_info_list[[1]][[3]], 1, function(x) x[2]/sum(x))
het_ratio <- rep(NA, times = length(geno_info_list[[1]][[1]]))
het_ratio[het_inds] <- test_ratio[het_inds]

# some rules for heterozygous genotypes: ratio bust be between 0.05

# let's see what the "bad" het inds are in each sample - these are hets called
#   using very skewed read counts


test_geno_list <- strsplit(combo_df_2[ , first_samp_ind], split = ':')
test_coverage <- as.numeric(unlist(lapply(test_geno_list, function(x) x[3])))
test_geno_calls <- unlist(lapply(test_geno_list, function(x) sum(as.numeric(unlist(strsplit(x[1], split = '/'))))))
test_allele_reads <- matrix(data = unlist(lapply(test_geno_list, function(x) as.numeric(unlist(strsplit(x[2], split = ','))))), byrow = T, ncol = 2, dimnames = list(rownames = NULL, colnames = c('n_REF', 'n_ALT')))

