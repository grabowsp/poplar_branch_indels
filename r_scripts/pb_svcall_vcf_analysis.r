# Script to look at PacBio SV calling VCF files

# LOAD LIBRARIES #

# LOAD DATA #
meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')


######
# SANDBOX

# New PacBio caller
# Process the Combo file
data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'
combo_file <- 'ref.ALLData.vcf'
combo_file_tot <- paste(data_dir, combo_file, sep = '')

combo_vcf_0 <- read.table(combo_file_tot, sep = '\t', stringsAsFactors = F)
combo_vcf_0$full_name <- paste(combo_vcf_0[,1], combo_vcf_0[,2], sep = '_')

combo_dup_inds <- which(duplicated(combo_vcf_0$full_name))
combo_vcf <- combo_vcf_0[-combo_dup_inds,]

combo_info_list <- strsplit(combo_vcf[ ,8], split = ';')
combo_sv_vec <- gsub('SVTYPE=', '', unlist(lapply(combo_info_list, function(x) x[1])))
combo_vcf$type <- combo_sv_vec

combo_vcf_indel <- combo_vcf[grep('DEL|INS', combo_vcf$type), ]

combo_vcf_indel$sv_length <- abs(sapply(combo_vcf_indel[,8],
  function(x) as.numeric(gsub('SVLEN=', '',
  unlist(strsplit(x, split = ';'))[3]))))

# functions to process the individual files

load_ind_pbnew_vcf <- function(in_file){
  tmp_vcf_0 <- read.table(in_file, sep = '\t', stringsAsFactors = F)
  tmp_vcf_0$full_name <- paste(tmp_vcf_0[,1], tmp_vcf_0[,2], sep = '_')
  tmp_dup_inds <- which(duplicated(tmp_vcf_0$full_name))
  tmp_vcf <- tmp_vcf_0[-tmp_dup_inds, ]
  tmp_info_list <- strsplit(tmp_vcf[ ,8], split = ';')
  tmp_sv_vec <- gsub('SVTYPE=', '', 
    unlist(lapply(tmp_info_list, function(x) x[1])))
  tmp_vcf$type <- tmp_sv_vec
  return(tmp_vcf)
}

# Test function
paxl_file <- 'ref.PAXL.vcf'
paxl_file_tot <- paste(data_dir, paxl_file, sep = '')
test_paxl <- load_ind_pbnew_vcf(paxl_file_tot)

gen_ind_uni_df <- function(indiv_df, combo_df){
  # generate data.frame of unique/not-shared positions in the individual vs 
  #  combofile
  indiv_share_inds <- which(indiv_df$full_name %in% combo_df$full_name)
  # combo_share_inds <- which(combo_df$full_name %in% indiv_df$full_name)
  indiv_bnd_inv_inds <- grep('BND|INV', indiv_df$type)
  indiv_uni <- indiv_df[-union(indiv_share_inds, indiv_bnd_inv_inds),]
  indiv_uni$dist_closest_combo <- NA
  # 
  for(chrom in unique(indiv_uni[,1])){
    tmp_combo_inds <- grep(chrom, combo_df[,1])
    tmp_samp_inds <- grep(chrom, indiv_uni[,1])
    for(i in tmp_samp_inds){
      indiv_uni$dist_closest_combo[i] <- min(
        abs(combo_df[tmp_combo_inds, 2] - indiv_uni[i,2]))
    }
  }
  indiv_uni$sv_length <- abs(sapply(indiv_uni[,8], 
    function(x) as.numeric(gsub('SVLEN=', '',
    unlist(strsplit(x, split = ';'))[3]))))
  return(indiv_uni)
}

# Test function
test_paxl_2 <- gen_ind_uni_df(indiv_df = test_paxl, combo_df = combo_vcf)

gen_nonoverlap_df <- function(uni_ind_df, combo_df, use_svsize = T, 
  dist_cut = 100, size_cut = 0){
  # function to generate data.frame of SVs in individual files that are not 
  #   present in the combo file
  if(use_svsize){
    tmp_nonover_inds <- which((uni_ind_df$dist_closest_combo - 
      (2*uni_ind_df$sv_length))> 0)
  } else {
    tmp_nonover_inds <- which(uni_ind_df$dist_closest_combo > dist_cut)
  }
  tmp_big_NO_inds <- intersect(tmp_nonover_inds, which(
    uni_ind_df$sv_length > size_cut))
  return(uni_ind_df[tmp_big_NO_inds, ])
}

# Test function
test_NO_paxl <- gen_nonoverlap_df(uni_ind_df = test_paxl_2, use_svsize = T)
# NOTE: this is with no size cutoff - will eventually want to use 50bp as cutoff

load_to_nonoverlap_df <- function(in_file, combo_df, use_svsize = T, 
  dist_cut = 100, size_cut = 0){
  # load individual file and process to get non-overlapping positions
  tmp_df_1 <- load_ind_pbnew_vcf(in_file = in_file)
  tmp_df_2 <- gen_ind_uni_df(indiv_df = tmp_df_1, combo_df = combo_df)
  tmp_df_3 <- gen_nonoverlap_df(uni_ind_df = tmp_df_2, combo_df = combo_df, 
                use_svsize = use_svsize, dist_cut = dist_cut, 
                size_cut = size_cut)
  return(tmp_df_3)
}

# Test function
paxl_full_process <- load_to_nonoverlap_df(in_file = paxl_file_tot, 
  combo_df = combo_vcf)

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

num_indiv_50bp_indels <-  unlist(lapply(indiv_NO_50bp_sites, function(x)
  length(grep('INS|DEL', x$type))))

gen_indel_raw_df <- function(indiv_df){
  # generate data.frame from raw VCF of insertions and deletions  and 
  #   include lengths
  indiv_bnd_inv_inds <- grep('BND|INV', indiv_df$type)
  indiv_indels <- indiv_df[-indiv_bnd_inv_inds,]
  indiv_indels$sv_length <- abs(sapply(indiv_indels[,8],
    function(x) as.numeric(gsub('SVLEN=', '',
    unlist(strsplit(x, split = ';'))[3]))))
  return(indiv_indels)
}

raw_indel_len_vcfs <- lapply(raw_vcfs, gen_indel_raw_df)

raw_50bp_vcfs <- lapply(raw_indel_len_vcfs, 
  function(x) x[which(x$sv_length >= 50), ])

tot_num_50bp_indels <- unlist(lapply(raw_50bp_vcfs, function(x)
  length(grep('INS|DEL', x$type))))

indiv_sv_df$tot_indels_50bp <- tot_num_50bp_indels
indiv_sv_df$indels_50bp_not_in_combo <- num_indiv_50bp_indels
indiv_sv_df$per_50bp_not_in_combo <- (indiv_sv_df$indels_50bp_not_in_combo /
  indiv_sv_df$tot_indels_50bp)

indiv_sv_df <- indiv_sv_df[order(indiv_sv_df$samp), ]

indiv_sv_out_file <- paste(data_dir, 'tmp/newPB_indiv_SV.txt', sep = '')
write.table(indiv_sv_df, file = indiv_sv_out_file, quote = F, sep = '\t',
  row.names = F, col.names = T)

#######
# Shared SVs NOT in combo VCF file
indiv_NO_site_names_all <- unlist(lapply(indiv_nonoverlap_sites, 
  function(x) x$full_name))

indiv_NO_site_tab <- table(table(indiv_NO_site_names_all))
indiv_NO_site_tab
#     1     2     3     4     5     6     7     8     9    10 
# 12062  2125  1193   734   548   453   339   235   162    58

indiv_NO_50bp_sites <- lapply(indiv_nonoverlap_sites, 
  function(x) x[which(x$sv_length >= 50), ])

indiv_NO_50bp_names_all <- unlist(lapply(indiv_NO_50bp_sites, 
  function(x) x$full_name))


tab_NO50 <- table(table(indiv_NO_50bp_names_all))
tab_NO50
#    1    2    3    4    5    6    7    8    9   10 
# 1635  200  100   69   40   35   30   18   11    6 

indiv_NO_tab_df <- data.frame(
  n_samps_with_SV_missing_in_combo = as.numeric(names(indiv_NO_site_tab)), 
  all_indels = as.vector(indiv_NO_site_tab), 
  indels_50bp = as.vector(tab_NO50), stringsAsFactors = F)

indiv_sv_tab_out_file <- paste(data_dir, 'tmp/newPB_indiv_SV_table.txt', 
  sep = '')
write.table(indiv_NO_tab_df, file = indiv_sv_tab_out_file, quote = F, 
  sep = '\t', row.names = F, col.names = T)

sum(indiv_NO_tab_df$all_indels[2:10])
# 5847
sum(indiv_NO_tab_df$indels_50bp[2:10])
# 509

n_indels_in_combo <- nrow(combo_vcf_indel)
# [1] 56464

n_indels_50bp_in_combo <- sum(combo_vcf_indel$sv_length >= 50)
# [1] 30665

sum(indiv_NO_tab_df$all_indels[2:10]) / n_indels_in_combo
# [1] 0.1035527

sum(indiv_NO_tab_df$indels_50bp[2:10]) / n_indels_50bp_in_combo
# [1] 0.01659873




weird_sites <- names(tab_NO50)[which(tab_NO50 == 10)]

weird_site_indiv_info <- list()
for(ws in weird_sites){
  weird_site_indiv_info[[ws]] <- matrix(data = unlist(
    lapply(indiv_nonoverlap_sites, function(x) x[grep(ws, x$full_name), ])), 
    byrow = T, ncol = ncol(indiv_nonoverlap_sites[[1]]))
}

#####
# Get Number of SVs from each file
all_vcf_files <- system(paste('ls ', data_dir, '*vcf', sep = ''), intern = T)

processed_vcfs_all <- lapply(all_vcf_files, load_ind_pbnew_vcf)
num_svs_all <- unlist(lapply(processed_vcfs_all, nrow))

sv_type_mat <- matrix(unlist(lapply(processed_vcfs_all, 
  function(x) table(x$type))), byrow = T, ncol = 4)
colnames(sv_type_mat) <- names(table(processed_vcfs_all[[1]]$type))

vcf_info_df <- data.frame(file = gsub(data_dir, '', all_vcf_files), 
  num_SVs = num_svs_all, sv_type_mat, stringsAsFactors = F)

info_out_file <- paste(data_dir, 'tmp/newPB_SV_gen_info.txt', sep = '')
write.table(vcf_info_df, file = info_out_file, quote = F, sep = '\t', 
  row.names = F, col.names= T)


####
# Sequencing Output vs # SVs and missingness

data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'

samp_meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
samp_meta <- read.table(samp_meta_file, header = T, stringsAsFactors = F)

sv_info_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/tmp/newPB_SV_gen_info.txt'

sv_info <- read.table(sv_info_file, header = T, stringsAsFactors = F, 
  sep = '\t') 

newPB_miss_geno_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/tmp/newPB_missing_genos.txt'
newPB_miss <- read.table(newPB_miss_geno_file, header = T, sep = '\t', 
  stringsAsFactors = F)

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

sv_output_out_file <- paste(data_dir, 'tmp/newPB_SV_seqOut.txt', sep = '')
write.table(sv_and_output, file = sv_output_out_file, quote = F, sep = '\t',
  row.names = F, col.names = T)

quit(save = 'no')
