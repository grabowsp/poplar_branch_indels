# Script for processing indel length info for
#  the 3 main SV calling pipelines we used for the poplar data

# LOAD PACKAGES #


# LOAD DATA #
samp_meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v3.0.txt'
samp_meta <- read.table(samp_meta_file, header = T, stringsAsFactors = F, 
  sep = '\t')

oldPB_data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/old_PB_SVcaller/tmp/'

newPB_data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/tmp/'

sniffles_data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/sniffles/tmp/'

# SET CONSTANTS #
bad_libs <- c('13.4', '14.1')
bad_lib_inds <- c()
for(bl in bad_libs){
  tmp_ind <- which(samp_meta$branch_name == bl)
  bad_lib_inds <- c(bad_lib_inds, tmp_ind)
}
samp_libs <- samp_meta$lib_name[-bad_lib_inds]

# SET VARIABLES #
oldPB_pre <- 'pop_branches_'
oldPB_SVLEN_combo_suf <- '_combo_length.INFO'
oldPB_SVLEN_solo_suf <- '_solo_length.INFO'

miss_oldPB_solo_lib <- 'PAZH'

newPB_pre <- 'pop_branches_newPB_'
newPB_SVLEN_combo_suf <- '_combo_length.INFO'
newPB_SVLEN_solo_suf <- '_solo_length.INFO'

snif_pre <- 'pop_branches_sniffles_'
snif_SVLEN_combo_suf <- '_combo_length.INFO'
snif_SVLEN_solo_suf <- '_solo_length.INFO'

# SET OUTPUT #
## directory for where to save combined-analysis files
combo_analysis_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/combined_analysis/'
processed_indel_file <- 'combined_processed_indel_length_info.rds'

############
# Process SV Length Data
process_PB_SVlength_data <- function(lib_vec, data_dir, data_pre, data_suf){
  tmp_length_list <- list()
  for(opbl in lib_vec){
    tmp_file <- paste(data_dir, data_pre, opbl, data_suf, sep = '')
    tmp_data <- read.table(tmp_file, header = T, stringsAsFactors = F)
    tmp_df <- data.frame(SVLEN = tmp_data$SVLEN, samp = opbl, 
      stringsAsFactors = F)
    tmp_df$SVLEN <- as.numeric(tmp_df$SVLEN)
    tmp_na_inds <- which(is.na(tmp_df$SVLEN))
    if(length(tmp_na_inds) > 0){
      tmp_df <- tmp_df[-which(is.na(tmp_df$SVLEN)),]
    }
    tmp_length_list[[opbl]] <- tmp_df
  }
  return(tmp_length_list)
}

process_SNIF_SVlength_data <- function(lib_vec, data_dir, data_pre, 
  data_suf, combo_infile = F){
  len_col_name <- 'SVLEN'
  if(combo_infile){ 
    len_col_name <- 'AVGLEN'
  }
  tmp_length_list <- list()
  for(opbl in lib_vec){
    tmp_file <- paste(data_dir, data_pre, opbl, data_suf, sep = '')
    tmp_data <- read.table(tmp_file, header = T, stringsAsFactors = F)
    del_inds <- grep('DEL', tmp_data$ALT)
    ins_inds <- grep('INS', tmp_data$ALT)
    tmp_data[del_inds, len_col_name] <- tmp_data[del_inds, len_col_name] * -1
    tmp_data_indel <- tmp_data[sort(union(del_inds, ins_inds)), ]
    tmp_df <- data.frame(SVLEN = tmp_data_indel[,len_col_name], samp = opbl,
      stringsAsFactors = F)
    tmp_length_list[[opbl]] <- tmp_df
  }
  return(tmp_length_list)
}

gen_SVlength_df <- function(SVlength_list, SV_caller){
  tmp_len_vec <- unlist(lapply(SVlength_list, function(x) x$SVLEN))
  tmp_samp_vec <- unlist(lapply(SVlength_list, function(x) x$samp))
  tmp_out_df <- data.frame(SVLEN = tmp_len_vec, samp = tmp_samp_vec,
  SV_caller = SV_caller, stringsAsFactors = F)
  return(tmp_out_df)
}

## Old PacBio Caller
old_combo_len_list <-  process_PB_SVlength_data(lib_vec = samp_libs, 
  data_dir = oldPB_data_dir, data_pre = oldPB_pre, 
  data_suf = oldPB_SVLEN_combo_suf)
 
oldPB_combo_df <- gen_SVlength_df(old_combo_len_list, 
  SV_caller = 'old_PB_combo')

old_solo_len_list <-  process_PB_SVlength_data(
  lib_vec = setdiff(samp_libs, miss_oldPB_solo_lib), 
  data_dir = oldPB_data_dir, data_pre = oldPB_pre, 
  data_suf = oldPB_SVLEN_solo_suf)
 
oldPB_solo_df <- gen_SVlength_df(old_solo_len_list, 
  SV_caller = 'old_PB_solo')

## New PacBio Caller
new_combo_len_list <- process_PB_SVlength_data(lib_vec = samp_libs, 
  data_dir = newPB_data_dir, data_pre = newPB_pre, 
  data_suf = newPB_SVLEN_combo_suf)

newPB_combo_df <- gen_SVlength_df(new_combo_len_list, 
  SV_caller = 'new_PB_combo')

new_solo_len_list <- process_PB_SVlength_data(lib_vec = samp_libs, 
  data_dir = newPB_data_dir, data_pre = newPB_pre, 
  data_suf = newPB_SVLEN_solo_suf)

newPB_solo_df <- gen_SVlength_df(new_solo_len_list, 
  SV_caller = 'new_PB_solo')

## SNIFFLES
snif_combo_len_list <- process_SNIF_SVlength_data(lib_vec = samp_libs,
  data_dir = sniffles_data_dir, data_pre = snif_pre,
  data_suf = snif_SVLEN_combo_suf, combo_infile = T)

snif_combo_df <- gen_SVlength_df(snif_combo_len_list,
  SV_caller = 'snif_combo')

snif_solo_len_list <- process_SNIF_SVlength_data(lib_vec = samp_libs,
  data_dir = sniffles_data_dir, data_pre = snif_pre,
  data_suf = snif_SVLEN_solo_suf, combo_infile = F)

snif_solo_df <- gen_SVlength_df(snif_solo_len_list,
  SV_caller = 'snif_solo')

## Combine results into single dataframe

combined_indel_length_df <- rbind(oldPB_combo_df, oldPB_solo_df, newPB_combo_df, newPB_solo_df, snif_combo_df, snif_solo_df)

tot_proc_indel_file <- paste(combo_analysis_dir, processed_indel_file, sep = '')
saveRDS(combined_indel_length_df, file = tot_proc_indel_file)

quit(save = 'no')


# SANDBOX #
#load a old-pb file
test_file <- paste(sniffles_data_dir, snif_pre, samp_libs[1], snif_SVLEN_solo_suf, sep = '')
test_data <- read.table(test_file, header = T, stringsAsFactors = F)
del_inds <- grep('DEL', test_data$ALT)
ins_inds <- grep('INS', test_data$ALT)
test_data$AVGLEN[del_inds] <- test_data$AVGLEN[del_inds] * -1
test_data_indel <- test_data[union(del_inds, ins_inds), ]


