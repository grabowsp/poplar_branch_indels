# Script to look at PacBio SV calling VCF files

# LOAD LIBRARIES #
## functions written to analyse PacBio SV VCFs
function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pb_SV_analysis_functions.r'
source(function_file)

# LOAD DATA #
meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')

data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'
combo_file <- 'ref.ALLData.vcf'
combo_file_tot <- paste(data_dir, combo_file, sep = '')
# combo_vcf_0 <- read.table(combo_file_tot, sep = '\t', stringsAsFactors = F)

## file containing amount of missing data in each vcf
newPB_miss_geno_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/tmp/newPB_missing_genos.txt'
newPB_miss <- read.table(newPB_miss_geno_file, header = T, sep = '\t', 
  stringsAsFactors = F)

# SET CONSTANTS #

# SET VARIABLES #

# SET OUTPUTS #
## table of type of SVs in each of the vcfs, including combo.vcf
info_out_sub <- 'tmp/newPB_SV_gen_info.txt'
info_out_file <- paste(data_dir, info_out_sub, sep = '')

## table with information about number of SVs/indels found in indiviual files
##  but NOT the combo.vcf file
indiv_sv_out_sub <- 'tmp/newPB_indiv_SV.txt'
indiv_sv_out_file <- paste(data_dir, indiv_sv_out_sub, sep = '')

## table of SV's NOT in combo.vcf showing how many are shared by different
##   numbers of individual SVs
indiv_sv_tab_out_sub <- 'tmp/newPB_indiv_SV_table.txt'
indiv_sv_tab_out_file <- paste(data_dir, indiv_sv_tab_out_sub, sep = '')

## table of sequencing output vs SVs in the vcfs
sv_output_out_sub <- 'tmp/newPB_SV_seqOut.txt'
sv_output_out_file <- paste(data_dir, sv_output_out_sub, sep = '')

##################
# For new PacBio Caller
# Get Number of SVs from each file
all_vcf_files <- system(paste('ls ', data_dir, '*vcf', sep = ''), intern = T)

processed_vcfs_all <- lapply(all_vcf_files, load_ind_pbnew_vcf)
num_svs_all <- unlist(lapply(processed_vcfs_all, nrow))

sv_type_mat <- matrix(unlist(lapply(processed_vcfs_all, 
  function(x) table(x$type))), byrow = T, ncol = 4)
colnames(sv_type_mat) <- names(table(processed_vcfs_all[[1]]$type))

vcf_info_df <- data.frame(file = gsub(data_dir, '', all_vcf_files), 
  num_SVs = num_svs_all, sv_type_mat, stringsAsFactors = F)

write.table(vcf_info_df, file = info_out_file, quote = F, sep = '\t', 
  row.names = F, col.names= T)

########
# Sequencing Output vs # SVs and missingness

# sv_info <- read.table(info_out_file, header = T, stringsAsFactors = F, 
#  sep = '\t') 
sv_info <- vcf_info_df

samp_meta$newPB_num_SVs <- NA

for(i in seq(nrow(samp_meta))){
  sv_ind <- grep(samp_meta$lib_name[i], sv_info$file)
  samp_meta$newPB_num_SVs[i] <- sv_info$num_SVs[sv_ind]
}

sv_and_output <- samp_meta[ , c('lib_name', 'tot_reads', 'singlepass_reads', 
  'X20kb_reads', 'newPB_num_SVs')]

sv_and_output$n_miss_in_combo <- NA

for(j in seq(nrow(sv_and_output))){
  miss_ind <- grep(sv_and_output$lib_name[j], newPB_miss$lib)
  sv_and_output$n_miss_in_combo[j] <- newPB_miss[miss_ind, 2]
}

write.table(sv_and_output, file = sv_output_out_file, quote = F, sep = '\t',
  row.names = F, col.names = T)

###########
# Calculate numer and percentage of Indels detected in individual files but
#   NOT in the combo.vcf file
## Process the Combo file
combo_vcf_indel <- make_combo_indel_df(combo_file_tot)

## Process Individual Files
indiv_vcf_files <- system(paste('ls ', data_dir, 'ref.P*vcf', sep = ''), 
  intern = T)

indiv_nonoverlap_sites <- lapply(indiv_vcf_files, load_to_nonoverlap_df, 
  combo_df = combo_vcf)

raw_vcfs <- lapply(indiv_vcf_files, load_ind_pbnew_vcf)
tot_num_indels <- unlist(lapply(raw_vcfs, function(x) 
  length(grep('INS|DEL', x$type))))

num_indiv_indels <- unlist(lapply(indiv_nonoverlap_sites, function(x)
  length(grep('INS|DEL', x$type))))

lib_names <- gsub('.vcf', '', gsub(paste(data_dir, 'ref.', sep = ''), '', 
  indiv_vcf_files))

samp_names <- c()
for(ln in lib_names){
  meta_ind <- which(samp_meta$lib_name == ln)
  tmp_name <- samp_meta$branch_name[meta_ind]
  samp_names <- c(samp_names, tmp_name)
}

indiv_sv_df <- data.frame(lib = lib_names, samp = samp_names, 
  tot_indels = tot_num_indels, indels_not_in_combo = num_indiv_indels, 
  stringsAsFactors = F)

indiv_sv_df$per_not_in_combo <- (indiv_sv_df$indels_not_in_comb / 
  indiv_sv_df$tot_indels)

indiv_NO_50bp_sites <- lapply(indiv_nonoverlap_sites,
  function(x) x[which(x$sv_length >= 50), ])

num_indiv_50bp_indels <-  unlist(lapply(indiv_NO_50bp_sites, function(x)
  length(grep('INS|DEL', x$type))))

raw_indel_len_vcfs <- lapply(raw_vcfs, gen_indel_raw_df)

raw_50bp_vcfs <- lapply(raw_indel_len_vcfs, 
  function(x) x[which(x$sv_length >= 50), ])

tot_num_50bp_indels <- unlist(lapply(raw_50bp_vcfs, function(x)
  length(grep('INS|DEL', x$type))))

indiv_sv_df$tot_indels_50bp <- tot_num_50bp_indels
indiv_sv_df$indels_50bp_not_in_combo <- num_indiv_50bp_indels
indiv_sv_df$per_50bp_not_in_combo <- (indiv_sv_df$indels_50bp_not_in_combo /
  indiv_sv_df$tot_indels_50bp)

# SVs in combo.vcf but not individual VCFs
indiv_miss_sites <- lapply(raw_vcfs, gen_ind_missing_df, 
  combo_df = combo_vcf_indel)

n_miss_in_indiv_from_combo <- unlist(lapply(indiv_miss_sites, nrow))

n_miss_50bp_in_indiv_from_combo <- unlist(lapply(indiv_miss_sites, 
  function(x) sum(x$sv_length >= 50)))

indiv_sv_df$indels_not_in_indiv <- n_miss_in_indiv_from_combo
indiv_sv_df$per_not_in_indiv <- (indiv_sv_df$indels_not_in_indiv / 
  indiv_sv_df$tot_indels)

indiv_sv_df$indels_50bp_not_in_indiv <- n_miss_50bp_in_indiv_from_combo
indiv_sv_df$per_50bp_not_in_indiv <- (indiv_sv_df$indels_50bp_not_in_indiv / 
  indiv_sv_df$tot_indels_50bp)

indiv_sv_df <- indiv_sv_df[order(indiv_sv_df$samp), ]

write.table(indiv_sv_df, file = indiv_sv_out_file, quote = F, sep = '\t',
  row.names = F, col.names = T)


##############
# Tally up SVs NOT in combo.vcf but found in individual files
indiv_NO_site_names_all <- unlist(lapply(indiv_nonoverlap_sites, 
  function(x) x$full_name))

indiv_NO_site_tab <- table(table(indiv_NO_site_names_all))
indiv_NO_site_tab
#     1     2     3     4     5     6     7     8     9    10 
# 12060  2125  1192   736   548   453   340   234   163    58

#indiv_NO_50bp_sites <- lapply(indiv_nonoverlap_sites, 
#  function(x) x[which(x$sv_length >= 50), ])

indiv_NO_50bp_names_all <- unlist(lapply(indiv_NO_50bp_sites, 
  function(x) x$full_name))

tab_NO50 <- table(table(indiv_NO_50bp_names_all))
tab_NO50
#    1    2    3    4    5    6    7    8    9   10 
# 1635  200  100   69   41   35   30   18   11    6

indiv_miss_site_names_all <- unlist(lapply(indiv_miss_sites, 
  function(x) x$full_name))

indiv_miss_site_tab <- table(table(indiv_miss_site_names_all))
indiv_miss_site_tab
#   1    2    3    4    5    6    7    8    9   10   11   12   14   16   18
#8316 3691 2209 1493 1192 1079 1054 1092 1534  579    2    4    5    2    1  
# 20  21 
#  1   1
# there are more than 10 because some positions contain 2 different sized
#  indels so they will be double-counted

indiv_miss_50bp_names_all <- unlist(lapply(indiv_miss_sites, 
  function(x) x$full_name[x$sv_length >= 50]))

indiv_miss_50bp_site_tab <- table(table(indiv_miss_50bp_names_all))
indiv_miss_50bp_site_tab
#    1    2    3    4    5    6    7    8    9   10   11   12   14   16   20
# 4006 1531  828  568  484  416  514  537  768  167    1    1    3    1    1
# 21
#  1

indiv_NO_tab_df <- data.frame(
  n_samps_with_SV_missing_in_combo = as.numeric(
    names(indiv_NO_site_tab))[c(1:10)], 
  all_indels = as.vector(indiv_NO_site_tab)[c(1:10)], 
  indels_50bp = as.vector(tab_NO50)[c(1:10)], 
  in_combo_not_indiv = as.vector(indiv_miss_site_tab)[c(1:10)],
  in_combo_not_indiv_50bp = as.vector(indiv_miss_50bp_site_tab)[c(1:10)],
  stringsAsFactors = F)
# note: this doesn't fully account for positions that are double-counted,
#   by those are very few postions so should affect overall conclustions

write.table(indiv_NO_tab_df, file = indiv_sv_tab_out_file, quote = F, 
  sep = '\t', row.names = F, col.names = T)

sum(indiv_NO_tab_df$all_indels[2:10])
# 5849
sum(indiv_NO_tab_df$indels_50bp[2:10])
# 510

n_indels_in_combo <- nrow(combo_vcf_indel)
# [1] 56794

n_indels_50bp_in_combo <- sum(combo_vcf_indel$sv_length >= 50)
# [1] 30859

sum(indiv_NO_tab_df$all_indels[2:10]) / n_indels_in_combo
# [1] 0.1029862

sum(indiv_NO_tab_df$indels_50bp[2:10]) / n_indels_50bp_in_combo
# [1] 0.01652678

sum(indiv_NO_tab_df$in_combo_not_indiv[2:10])
# [1] 13923
sum(indiv_NO_tab_df$in_combo_not_indiv[2:10]) / n_indels_in_combo
# [1] 0.2451491

sum(indiv_NO_tab_df$in_combo_not_indiv_50bp[2:10])
# [1] 5813
sum(indiv_NO_tab_df$in_combo_not_indiv_50bp[2:10]) / n_indels_50bp_in_combo
# [1] 0.1883729

quit(save = 'no')
########
# Sandbox

weird_sites <- names(tab_NO50)[which(tab_NO50 == 10)]

weird_site_indiv_info <- list()
for(ws in weird_sites){
  weird_site_indiv_info[[ws]] <- matrix(data = unlist(
    lapply(indiv_nonoverlap_sites, function(x) x[grep(ws, x$full_name), ])), 
    byrow = T, ncol = ncol(indiv_nonoverlap_sites[[1]]))
}

