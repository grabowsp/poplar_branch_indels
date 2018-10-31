# Funcitons used for analysis of SV's called using the PacBio pipeline

load_ind_pbnew_vcf <- function(in_file){
  # Load individual (not combo) VCF and extract the SV type from the info
  #  field; designed for VCFs produced by PacBio v2 (Sept. 2018) pipeline
  # INPUTS
  # in_file = full file path for the individual VCF to import
  # OUTPUTS
  # data.frame with the fields from the VCF and SV type column
  ###########
  tmp_vcf_0 <- read.table(in_file, sep = '\t', stringsAsFactors = F)
  tmp_vcf_0$full_name <- paste(tmp_vcf_0[,1], tmp_vcf_0[,2], sep = '_')
#  tmp_dup_inds <- which(duplicated(tmp_vcf_0$full_name))
#  tmp_vcf <- tmp_vcf_0[-tmp_dup_inds, ]
  tmp_vcf <- tmp_vcf_0
  tmp_info_list <- strsplit(tmp_vcf[ ,8], split = ';')
  tmp_sv_vec <- gsub('SVTYPE=', '', 
    unlist(lapply(tmp_info_list, function(x) x[1])))
  tmp_vcf$type <- tmp_sv_vec
  return(tmp_vcf)
}

# Test function
## data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'
## paxl_file <- 'ref.PAXL.vcf'
## paxl_file_tot <- paste(data_dir, paxl_file, sep = '')
## test_paxl <- load_ind_pbnew_vcf(paxl_file_tot)

gen_indel_raw_df <- function(indiv_df){
  # generate data.frame from raw VCF of insertions and deletions  and 
  #   include lengths
  # INPUTS #
  # indiv_df = data.frame from VCF file, can come from load_ind_pbnew_vcf() 
  #   function
  # OUTPUTS #
  # data.frame including all InDels and the their lengths
  ###################
  indiv_bnd_inv_inds <- grep('BND|INV', indiv_df$type)
  indiv_indels <- indiv_df[-indiv_bnd_inv_inds,]
  indiv_indels$sv_length <- abs(sapply(indiv_indels[,8],
    function(x) as.numeric(gsub('SVLEN=', '',
    unlist(strsplit(x, split = ';'))[3]))))
  return(indiv_indels)
}

# test_paxl_raw <- gen_indel_raw_df(test_paxl)

# NEXT STEPS: Annotate this function, including example, update other functions
#  to include info about this function; then continue with function at end of
#  page
make_combo_indel_df <- function(combo_file){
  # Function to load and process combo VCF so that have SV length for all
  #  InDels
  # INPUTS #
  # combo_file = full path for combined VCF
  # OUTPUTS #
  # data.frame for InDels that includes the length of the SVs
  #############3
  combo_vcf_0 <- read.table(combo_file, sep = '\t', stringsAsFactors = F)
  combo_vcf_0$full_name <- paste(combo_vcf_0[,1], combo_vcf_0[,2], sep = '_')
#  combo_dup_inds <- which(duplicated(combo_vcf_0$full_name))
#  combo_vcf <- combo_vcf_0[-combo_dup_inds,]
  combo_vcf <- combo_vcf_0
  combo_info_list <- strsplit(combo_vcf[ ,8], split = ';')
  combo_sv_vec <- gsub('SVTYPE=', '', unlist(lapply(combo_info_list,
    function(x) x[1])))
  combo_vcf$type <- combo_sv_vec
  combo_vcf_indel <- combo_vcf[grep('DEL|INS', combo_vcf$type), ]
  combo_vcf_indel$sv_length <- abs(sapply(combo_vcf_indel[,8],
    function(x) as.numeric(gsub('SVLEN=', '',
    unlist(strsplit(x, split = ';'))[3]))))
  return(combo_vcf_indel)
}

# data_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/'
# combo_file <- 'ref.ALLData.vcf'
# combo_file_tot <- paste(data_dir, combo_file, sep = '')
# test_combo <- make_combo_indel_df(combo_file_tot) 

gen_ind_uni_df <- function(indiv_df, combo_df){
  # generate data.frame of unique/not-shared positions in the individual vs 
  #  combofile; Combo file needs to be loaded and processed separately
  #  Output only contains insertions and deletions because these are the
  #    most prevalent and are the easiest to extract SV size from
  # INPUTS #
  # indiv_df = data.frame based on VCF; can be generated using 
  #              load_ind_pbnew_vcf()
  # combo_df = data.frame containing info from combined VCF, can be generated
  #              using the make_combo_indel_df() function
  # OUTPUTS #
  # data.frame containing ONLY INSertions and DELetions that are NOT at
  #   positions shared by the combo_df file; also includes the size of the
  #   SVs
  ###############
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
## test_paxl_2 <- gen_ind_uni_df(indiv_df = test_paxl, combo_df = combo_vcf)

gen_nonoverlap_df <- function(uni_ind_df, use_svsize = T,
  dist_cut = 100, size_cut = 0){
  # function to generate data.frame of SVs in individual files that are not 
  #   present or within a certain distance from SVs in the combo file
  # INPUTS #
  # uni_ind_df = data.frame containing SVs at positions NOT shared by the 
  #                combo file used as comparison; can be generated using
  #                gen_ind_uni_df() function
  # use_svsize = use the SV size to determine distance cutoff for calling a SV
  #                as truely missing from the combo file; if <T>, then
  #                next-closest SV in combo file must be a a distance greater 
  #                than 2X the size of the SV; if <F>, then use <dist_cut>
  # dist_cut = distance cutoff used for calling a SV as missing in combo file;
  #              if <use_svsize = F>, then next-closest SV in combo file must
  #              be greater than <dist_cut>
  # size_cut = size cutoff used for including SVs in analysis; useful if only
  #              want to look at SVs above a certain size
  # OUTPUTS #
  # data.frame including only InDels in uni_ind_df greater than the chosed 
  #   distance cutoff
  ################  
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
## test_NO_paxl <- gen_nonoverlap_df(uni_ind_df = test_paxl_2, use_svsize = T)

load_to_nonoverlap_df <- function(in_file, combo_df, use_svsize = T,
  dist_cut = 100, size_cut = 0){
  # Function that loads and processes an individual VCF to get sites in
  #   individual VCF NOT present in the combined vcf
  # INPUTS
  # in_file = full file path for the individual VCF to import
  # combo_df = data.frame from combined VCF used for comparison; can be 
  #              generated using the make_combo_indel_df() function
  # use_svsize = use the SV size to determine distance cutoff for calling a SV
  #                as truely missing from the combo file; if <T>, then
  #                next-closest SV in combo file must be a a distance greater 
  #                than 2X the size of the SV; if <F>, then use <dist_cut>
  # dist_cut = distance cutoff used for calling a SV as missing in combo file;
  #              if <use_svsize = F>, then next-closest SV in combo file must
  #              be greater than <dist_cut>
  # size_cut = size cutoff used for including SVs in analysis; useful if only
  #              want to look at SVs above a certain size
  # OUTPUTS #
  # data.frame including only InDels in uni_ind_df greater than the chosed 
  #   distance cutoff
  ##########################
  tmp_df_1 <- load_ind_pbnew_vcf(in_file = in_file)
  tmp_df_2 <- gen_ind_uni_df(indiv_df = tmp_df_1, combo_df = combo_df)
  tmp_df_3 <- gen_nonoverlap_df(uni_ind_df = tmp_df_2,
                use_svsize = use_svsize, dist_cut = dist_cut,
                size_cut = size_cut)
  return(tmp_df_3)
}

# Test function
## paxl_full_process <- load_to_nonoverlap_df(in_file = paxl_file_tot, 
##  combo_df = combo_vcf)

gen_ind_missing_df <- function(indiv_df, combo_df, use_svsize = T,
  dist_cut = 100, size_cut = 0){
  # Function to generate data.frame of SVs in the combo.vcf 
  #   but missing from the individual file
  # INPUTS #
  # indiv_df = data.frame from VCF file, can come from load_ind_pbnew_vcf() 
  #   function
  # combo_df = data.frame including info about indels from the combined VCF;
  #              can be generated using make_combo_indel_df() function
  # use_svsize = use the SV size to determine distance cutoff for calling a SV
  #                as truely missing from the individual file; if <T>, then
  #                next-closest SV in individual file must be distance greater 
  #                than 2X the size of the SV; if <F>, then use <dist_cut>
  # dist_cut = distance cutoff used for calling a SV as missing in indiv file;
  #              if <use_svsize = F>, then next-closest SV in indiv file must
  #              be greater than <dist_cut>
  # size_cut = size cutoff used for including SVs in analysis; useful if only
  #              want to look at SVs above a certain size
  # OUTPUTS #
  # data.frame of InDels that are in combo VCF but not in the individual
  #   VCF
  ####################
  combo_share_inds <- which(combo_df$full_name %in% indiv_df$full_name)
  combo_miss_df <- combo_df[-combo_share_inds, ]
  combo_miss_df$dist_closest_indiv <- NA
  # 
  for(chrom in unique(combo_miss_df[,1])){
    tmp_combo_inds <- grep(chrom, combo_miss_df[,1])
    tmp_samp_inds <- grep(chrom, indiv_df[,1])
    for(i in tmp_combo_inds){
      combo_miss_df$dist_closest_indiv[i] <- min(
        abs(indiv_df[tmp_samp_inds, 2] - combo_miss_df[i,2]))
    }
  }
  if(use_svsize){
    tmp_nonover_inds <- which((combo_miss_df$dist_closest_indiv -
      (2*combo_miss_df$sv_length))> 0)
  } else {
    tmp_nonover_inds <- which(combo_miss_df$dist_closest_indiv > dist_cut)
  }
  tmp_big_NO_inds <- intersect(tmp_nonover_inds, which(
    combo_miss_df$sv_length > size_cut))
  return(combo_miss_df[tmp_big_NO_inds, ])
}

# paxl_indiv_miss <- gen_ind_missing_df(indiv_df = test_paxl, 
#  combo_df = test_combo)


