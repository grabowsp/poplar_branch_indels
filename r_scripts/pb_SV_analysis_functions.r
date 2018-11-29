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

make_combo_indel_df <- function(combo_file, n_head_lines = 100){
  # Function to load and process combo VCF so that have SV length for all
  #  InDels
  # INPUTS #
  # combo_file = full path for combined VCF
  # n_head_lines = number of lines to scan in to get the full header from the
  #                  vcf
  # OUTPUTS #
  # data.frame for InDels that includes the length of the SVs
  #############3
  combo_vcf_0 <- read.table(combo_file, sep = '\t', stringsAsFactors = F)
  # add header to data.frame
  combo_vcf_header <- scan(combo_file_tot, nlines = n_head_lines, 
    what = 'character', sep = '\n', quiet = T)
  col_info_ind <- grep('#CHROM', combo_vcf_header, fixed = T)
  col_info <- gsub('#', '', unlist(strsplit(combo_vcf_header[col_info_ind], 
    split = '\t')), fixed = T)
  colnames(combo_vcf_0) <- col_info
  #
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
  combo_vcf_indel$full_name_long <- paste(combo_vcf_indel$full_name, 
    combo_vcf_indel$type, combo_vcf_indel$sv_length, sep = '_')
  combo_vcf_indel <- combo_vcf_indel[-which(
    duplicated(combo_vcf_indel$full_name_long)), ]
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
## test_paxl_2 <- gen_ind_uni_df(indiv_df = test_paxl, combo_df = test_combo)

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
##  combo_df = test_combo)

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

remove_missing_combo_SVs <- function(combo_df){
  # Remove any SVs from process combined VCF that have missing genotypes in any
  #   samples
  # INPUTS #
  # combo_df = data.frame including info about indels from the combined VCF;
  #              can be generated using make_combo_indel_df() function
  # OUTPUTS #
  # data.frame with same info as combo_df but with any SVs with missing
  #   genotype data removed
  ##########
  first_samp_ind <- which(colnames(combo_df) == 'FORMAT') + 1
  last_samp_ind <- which(colnames(combo_df) == 'full_name') - 1
  #
  miss_geno_inds <- c()
  for(sn in c(first_samp_ind:last_samp_ind)){
    tmp_inds <- grep('.', combo_df[,sn], fixed = T)
    miss_geno_inds <- c(miss_geno_inds, tmp_inds)
  }
  #
  miss_inds_toremove <- unique(miss_geno_inds)
  #
  combo_df_2 <- combo_df[-sort(miss_inds_toremove), ]
  return(combo_df_2)
}

# test_combo_filt <- remove_missing_combo_SVs(test_combo)

make_ind_geno_info_list <- function(geno_info_vec, sv_names = NULL){
  # Function to generate list of info about the genotypes of each sample in
  #   combo_df which is derived from a combined VCF
  # INPUTS #
  # geno_info_vec = vector of genotype info with three types of information
  #                   separated by colons: 1) genotype in X/Z form where
  #                   0 = reference allele and 1 = alternate allele; 
  #                   2) number of reads for each allele
  #                   3) total number of reads at that SV
  # sv_names = the names of the SVs; can usually come from the 
  #              'full_name_long' column of the combo_df that contains 
  #              the geno_info_vec
  # OUTPUTS #
  # list containing 3 elements:
  #  [[1]] numerical genotype with 0 = RR, 1 = RA, 2 = AA; 
  #  [[2]] vector of total coverage/number of reads at that SV; and 
  #  [[3]] matrix containing number of reads supporting each allele, with
  #    rows = SV and columns = allele, first column = R, second column = A
 ###########
  tmp_geno_list <- strsplit(geno_info_vec, split = ':')
  tmp_coverage <- as.numeric(unlist(lapply(tmp_geno_list, function(x) x[3])))
  names(tmp_coverage) <- sv_names
  tmp_num_genos <- unlist(lapply(tmp_geno_list, function(x)
    sum(as.numeric(unlist(strsplit(x[1], split = '/'))))))
  names(tmp_num_genos) <- sv_names
  tmp_allele_reads <- matrix(
    data = unlist(lapply(tmp_geno_list, function(x)
      as.numeric(unlist(strsplit(x[2], split = ','))))),
    byrow = T, ncol = 2,
    dimnames = list(rownames = sv_names, colnames = c('n_REF', 'n_ALT')))
  geno_info_list <- list()
  geno_info_list[[1]] <- tmp_num_genos
  geno_info_list[[2]] <- tmp_coverage
  geno_info_list[[3]] <- tmp_allele_reads
  return(geno_info_list)
}

# test_ind_geno_info <- make_ind_geno_info_list(test_combo_filt[,12], 
#   sv_names = test_combo_filt$full_name_long)

make_allsamp_geno_info_list <- function(combo_df){
  # Make list containing genotype info for all samples in a combo_df based
  #  on combined VCF, uses the make_ind_geno_info_list() function to generate
  #  the genotype info for each sample
  # INPUTS #
  # combo_df = data.frame including info about indels from the combined VCF;
  #              should be first generated using make_combo_indel_df() function 
  #              or then filtered using the remove_missing_combo_SVs() function
  # OUTPUTS #
  # list with each element containing sub-list with genotype info for a sample.
  # Genotype info for each sample contains 3 types of info: 
  #  [[1]] numerical genotype with 0 = RR, 1 = RA, 2 = AA; 
  #  [[2]] vector of total coverage/number of reads at that SV; and 
  #  [[3]] matrix containing number of reads supporting each allele, with
  #    rows = SV and columns = allele, first column = R, second column = A
  #################
  first_samp_ind <- which(colnames(combo_df) == 'FORMAT') + 1
  last_samp_ind <- which(colnames(combo_df) == 'full_name') - 1
  tot_info_list <- list()
  for(s_ind in c(first_samp_ind:last_samp_ind)){
    samp_name <- colnames(combo_df)[s_ind]
    tmp_samp_geno_list <- make_ind_geno_info_list(combo_df[ , s_ind], 
      sv_names = combo_df$full_name_long)
    tot_info_list[[samp_name]] <- tmp_samp_geno_list
  }
  return(tot_info_list)
}

# test_full_geno_list <- make_allsamp_geno_info_list(test_combo_filt)

call_ind_genotype <- function(ind_geno_list, min_sv_coverage = 10, 
  het_ratio_cut = 0.15, min_minor_allele_count = 2, max_hom_ratio = 0.05){
  # Function for calling genotypes at each SV using a list containg genotype
  #  info for that sample
  # INPUTS #
  # ind_geno_list = list containing genotype info for a sample. Has 3 elements:
  #                   [[1]] numerical genotype with 0 = RR, 1 = RA, 2 = AA; 
  #                   [[2]] vector of total coverage/number of reads at SV; 
  #                   [[3]] matrix containing number of reads supporting each 
  #                     allele, with rows = SV and columns = allele, first 
  #                     column = R, second column = A
  # min_sv_coverage = the required minimum coverage/total read depth for an SV;
  #                     if coverage is below <min_sv_coverage>, then call
  #                     genotype as NA
  # het_ratio_cut = the ratio of A/R or R/A above which the genotype is called
  #                   as heterozygous (or "1")
  #                   ex: if <het_ratio_cut> = 0.15, any SV's with a R/A ratio
  #                   between 0.15 and 0.85 are called as heterozygous
  # min_minor_allele_count = the required minimum number of reads for the
  #                            allele with fewer reads in order for a genotype
  #                            to be called heterozygous. If read ratio is
  #                            above <het_ratio_cut> but below 
  #                            <min_minor_allele_count>, then genotype is NA
  # max_hom_ratio = the maximum read ratio for a homozygous genotype; if
  #                   the read ratio is between <max_hom_ratio> and 
  #                   <het_ratio_cut>, then genotype is NA 
  # OUTPUTS #
  # vector of numerical genotypes using cutoffs and parameters in the function;
  #   0 = homozygous Ref; 1 = heterozygous; 2 = homozygous Alt
  ##############
  # Call initial genotypes based of R/A read ratio
  tmp_ref_ratio <- apply(ind_geno_list[[3]], 1, function(x) x[1]/sum(x))
  tmp_genos <- rep(NA, times = length(tmp_ref_ratio))
  names(tmp_genos) <- names(ind_geno_list[[1]])
  tmp_call_hom_alt <- which(tmp_ref_ratio < het_ratio_cut)
  tmp_call_hom_ref <- which(tmp_ref_ratio > (1 - het_ratio_cut))
  tmp_call_het <- setdiff(seq(length(tmp_genos)),
    union(tmp_call_hom_ref, tmp_call_hom_alt))
  tmp_genos[tmp_call_hom_ref] <- 0
  tmp_genos[tmp_call_hom_alt] <- 2
  tmp_genos[tmp_call_het] <- 1
  # Filter genotypes based on cutoffs
  tmp_too_low_inds <- which(ind_geno_list[[2]] < 10)
  tmp_min_ratio <- apply(ind_geno_list[[3]], 1, function(x) min(x)/sum(x))
  tmp_hom_inds <- setdiff(seq(length(tmp_genos)), tmp_call_het)
  tmp_too_skewed_homs <- tmp_hom_inds[
    which(tmp_min_ratio[tmp_hom_inds] > 0.05)]
  tmp_low_min_allele <- apply(ind_geno_list[[3]], 1, min)
  tmp_too_low_het <- intersect(which(tmp_low_min_allele <= 2), tmp_call_het)
  tmp_na_inds <- union(tmp_too_low_inds, 
    union(tmp_too_skewed_homs, tmp_too_low_het))
  tmp_genos[tmp_na_inds] <- NA
  return(tmp_genos)
}

# test_ind_num_geno <- call_ind_genotype(test_ind_geno_info, 
#   min_sv_coverage = 10, het_ratio_cut = 0.15, min_minor_allele_count = 2, 
#   max_hom_ratio = 0.05)

call_allsamp_genotypes <- function(geno_info_list,  min_sv_coverage = 10,
  het_ratio_cut = 0.15, min_minor_allele_count = 2, max_hom_ratio = 0.05){
  # Function to call numeric genotype for all samples in geno_info_list using
  #   the given parameters and cutoffs
  # INPUTS #
  # geno_info_list = list of genotype info for each sample generated by the
  #                    make_allsamp_geno_info_list() function
  # min_sv_coverage = the required minimum coverage/total read depth for an SV;
  #                     if coverage is below <min_sv_coverage>, then call
  #                     genotype as NA
  # het_ratio_cut = the ratio of A/R or R/A above which the genotype is called
  #                   as heterozygous (or "1")
  #                   ex: if <het_ratio_cut> = 0.15, any SV's with a R/A ratio
  #                   between 0.15 and 0.85 are called as heterozygous
  # min_minor_allele_count = the required minimum number of reads for the
  #                            allele with fewer reads in order for a genotype
  #                            to be called heterozygous. If read ratio is
  #                            above <het_ratio_cut> but below 
  #                            <min_minor_allele_count>, then genotype is NA
  # max_hom_ratio = the maximum read ratio for a homozygous genotype; if
  #                   the read ratio is between <max_hom_ratio> and 
  #                   <het_ratio_cut>, then genotype is NA 
  # OUTPUTS #
  # matrix containing filtered, numerical genotypes for each sample;
  #   columns = samples, rows = SVs 
  #####
  allsamp_genos_list <- lapply(geno_info_list, call_ind_genotype, 
    min_sv_coverage = min_sv_coverage, het_ratio_cut = het_ratio_cut, 
    min_minor_allele_count = min_minor_allele_count, 
    max_hom_ratio = max_hom_ratio)
  samp_names <- names(geno_info_list)
  allsamp_genos_mat <- matrix(data = unlist(allsamp_genos_list), 
    ncol = length(geno_info_list), byrow = F)
  colnames(allsamp_genos_mat) <- samp_names
  rownames(allsamp_genos_mat) <- names(geno_info_list[[1]][[1]])
  return(allsamp_genos_mat)
}

# test_allsamp_num_genos <- call_allsamp_genotypes(test_full_geno_list)

filt_geno_mat <- function(geno_mat, max_nas = NULL, min_length = NULL, 
  max_length = NULL, sv_type = NULL){
  # Function to filter out SV's with more than <max_nas> NA's or SV's that are
  #   longer or shorter than a size cutoff
  # INPUTS #
  # geno_mat = matrix of numeric genotypes; cols = samples, rows = SVs; best
  #              if generated using call_allsamp_genotypes() function
  # max_nas = the maximum number of NAs allowed at an SV across all samples;
  #             SV's with more NA's are removed
  # min_length = the minimum length of an SV to be included; smaller SVs are
  #                removed
  # max_length = the maximum length of an SV to be included; larger SVs are
  #                removed
  # OUTPUTS #
  # matrix of numerical genotypes
  #######3
  tmp_bad_inds <- c()
  if(length(max_nas) > 0){
    tmp_num_nas <- apply(geno_mat, 1, function(x) sum(is.na(x)))
    tmp_bad_inds <- c(tmp_bad_inds, which(tmp_num_nas > max_nas))
  }
  sv_length_vec <- as.numeric(unlist(lapply(
    strsplit(rownames(geno_mat), split = '_'), function(x) x[[length(x)]])))
  sv_type_vec <- unlist(lapply(
    strsplit(rownames(geno_mat), split = '_'), function(x) x[[length(x)-1]]))
  if(length(min_length) > 0){
    tmp_small_svs <- which(sv_length_vec < min_length)
    tmp_bad_inds <- c(tmp_bad_inds, tmp_small_svs)
  }
  if(length(max_length) > 0){
    tmp_big_svs <- which(sv_length_vec > max_length)
    tmp_bad_inds <- c(tmp_bad_inds, tmp_big_svs)
  }
  if(length(sv_type) > 0){
    tmp_svtype <- which(sv_type_vec != sv_type)
    tmp_bad_inds <- c(tmp_bad_inds, tmp_svtype)
  }
  tmp_bad_inds <- sort(tmp_bad_inds)
  geno_mat_filt <- geno_mat[-tmp_bad_inds, ]
  return(geno_mat_filt)
}

# test_filt_allsamp_genos <- filt_geno_mat(test_allsamp_num_genos, max_nas = 0)
# test_filt_genos_50 <- filt_geno_mat(test_allsamp_num_genos, max_nas = 0, 
#   min_length = 50)

tally_singleton_svs <- function(geno_mat){
  # Function for tallying singleton SV genotypes by sample
  #  note: written for SVs with only 2 genotype class, NOT 3 genotype classes
  # INPUTS #
  # geno_mat = matrix of numeric genotypes; cols = samples, rows = SVs; best
  #              if generated using call_allsamp_genotypes() function and
  #              filtered with filt_geno_mat() function
  # OUTPUT #
  # list with three elements, each containing tables of the number of 
  #   singleton genotypes for each genotype class (0, 1, 2) in each library
  ######
  n_geno_types <- apply(geno_mat, 1, function(x) length(table(x)))
  two_geno_inds <- which(n_geno_types == 2)
  #
  test_n_0 <- apply(geno_mat[two_geno_inds, ], 1, function(x) sum(x == 0))
  test_n_1 <- apply(geno_mat[two_geno_inds, ], 1, function(x) sum(x == 1))
  test_n_2 <- apply(geno_mat[two_geno_inds, ], 1, function(x) sum(x == 2))
  #
  homRef_singleton_SVs <- two_geno_inds[which(test_n_0 == 1)]
  homRef_sing_lib <- apply(geno_mat[homRef_singleton_SVs,], 1,
    function(x) colnames(geno_mat)[which(x == 0)])
  homAlt_singleton_SVs <- two_geno_inds[which(test_n_2 == 1)]
  homAlt_sing_lib <- apply(geno_mat[homAlt_singleton_SVs,], 1,
    function(x) colnames(geno_mat)[which(x == 2)])
  het_singleton_SVs <- two_geno_inds[which(test_n_1 == 1)]
  het_sing_lib <- apply(geno_mat[het_singleton_SVs,], 1,
    function(x) colnames(geno_mat)[which(x == 1)])
  #
  sing_tab_list <- list()
  sing_tab_list[['homRef_singletons']] <- table(homRef_sing_lib)
  sing_tab_list[['het_singletons']] <- table(het_sing_lib)
  sing_tab_list[['homAlt_singletons']] <- table(homAlt_sing_lib)
  return(sing_tab_list)
}

get_branch_lib_names <- function(branch_name_vec, meta){
  # Get the library names of a vector of branch names
  # INPUTS #
  # branch_name_vec = vector of branch names; should be numeric (ex: 13.1)
  #                     instead of character because that's the form in the
  #                     metadata
  # meta = sample metadata dataframe
  # OUTPUTS #
  # vector of the library names in the same order as <branch_name_vec>
  ########
  meta_ind_vec <- c()
  for(bn in branch_name_vec){
    tmp_ind <- which(meta$branch_name == bn)
    meta_ind_vec <- c(meta_ind_vec, tmp_ind)
  }
  lib_vec <- meta$lib_name[meta_ind_vec]
  return(lib_vec)
}

get_lib_col_inds <- function(geno_mat, lib_name_vec){
  # Get the column indices from geno_mat of the libraries in lib_name_vec
  # INPUTS #
  # geno_mat = matrix of numeric genotypes; cols = samples, rows = SVs; best
  #              if generated using call_allsamp_genotypes() function and
  #              filtered with filt_geno_mat() function
  # lib_name_vec = vector of the library names that want to find indices of
  #                  in <geno_mat>; can be generated with 
  #                  get_branch_lib_names() function
  # OUTPUT #
  # vector of the column indices in <geno_mat> that contain the data for the
  #   libraries in <lib_name_vec>
  ########3
  col_ind_vec <- c()
  mat_col_names <- colnames(geno_mat)
  for(lnv in lib_name_vec){
    tmp_ind <- which(mat_col_names == lnv)
    col_ind_vec <- c(col_ind_vec, tmp_ind)
  }
  return(col_ind_vec)
}

get_samegeno_inds <- function(geno_mat, branch_name_vec, genotype, meta, 
  n_miss = 0){
  # get the indices of SVs that are all the same, given genotype in a set of
  #   samples
  # INPUTS #
  # geno_mat = matrix of numeric genotypes; cols = samples, rows = SVs; best
  #              if generated using call_allsamp_genotypes() function and
  #              filtered with filt_geno_mat() function
  # branch_name_vec = vector of branch names; should be numeric (ex: 13.1)
  #                     instead of character because that's the form in the
  #                     metadata
  # genotype = the numerical genotype that are looking for in the samples
  # meta = sample metadata dataframe
  # n_miss = the number of samples that can have a different genotype. Usually
  #            want <n_miss> = 0, but if worried about miss-called genotypes,
  #            then can try increasing <n_miss>
  # OUTPUT #
  # vector of indices of the SVs in geno_mat that are all the target genotype
  #   in the target samples
  #############
  tmp_lib_names <- get_branch_lib_names(branch_name_vec, meta)
  tmp_data_col_inds <- get_lib_col_inds(geno_mat, tmp_lib_names) 
  #
  test_geno_count <- apply(geno_mat[ , tmp_data_col_inds], 1, 
    function(x) sum(x == genotype))
  target_num <- length(branch_name_vec) - n_miss
  invar_geno_inds <- which(test_geno_count == target_num)
  return(invar_geno_inds)
}

get_shared_unique_het_inds <- function(geno_mat, branch_name_vec, test_names, 
  meta){
  # Get the indices of SVs that are shared between but unique to a set of
  #  samples - they are heterozygous in <test_names> and homozygous in all
  #  the rest
  # INPUTS #
  # geno_mat = matrix of numeric genotypes; cols = samples, rows = SVs; best
  #              if generated using call_allsamp_genotypes() function and
  #              filtered with filt_geno_mat() function
  # branch_name_vec = vector of branch names; should be numeric (ex: 13.1)
  #                     instead of character because that's the form in the
  #                     metadata
  # test_names = vector of names in <branch_name_vec> that are looking for
  #                shared, unique SVs in
  # genotype = the numerical genotype that are looking for in the samples
  # meta = sample metadata dataframe
  # OUTPUT #
  # vector of indices of the SVs in geno_mat that are shared, unique
  #   heterozygous SVs in the <test_names> samples
  ########
  alt_names <- setdiff(branch_name_vec, test_names)
  tot_het_inds <- get_samegeno_inds(geno_mat = geno_mat, 
    branch_name_vec = test_names, genotype = 1, meta = meta, 
    n_miss = 0)
  other_homRef <- get_samegeno_inds(geno_mat = geno_mat, 
    branch_name_vec = alt_names, genotype = 0, meta = meta, n_miss = 0)
  other_homAlt <- get_samegeno_inds(geno_mat = geno_mat, 
    branch_name_vec = alt_names, genotype = 2, meta = meta, n_miss = 0)
  #
  spec_hets <- intersect(tot_het_inds, union(other_homRef, other_homAlt))
  return(spec_hets)
}

get_shared_unique_hom_inds <- function(geno_mat, branch_name_vec, test_names,
  meta){
  # Get the indices of SVs that are homozygous ONLY in <test_names> and
  #  heterozygous in the rest of the samples
  # INPUTS #
  # geno_mat = matrix of numeric genotypes; cols = samples, rows = SVs; best
  #              if generated using call_allsamp_genotypes() function and
  #              filtered with filt_geno_mat() function
  # branch_name_vec = vector of branch names; should be numeric (ex: 13.1)
  #                     instead of character because that's the form in the
  #                     metadata
  # test_names = vector of names in <branch_name_vec> that are looking for
  #                shared, unique homozygous SV genotypes in
  # genotype = the numerical genotype that are looking for in the samples
  # meta = sample metadata dataframe
  # OUTPUT #
  # vector of indices of the SVs in geno_mat that are shared, unique
  #   heterozygous SVs in the <test_names> samples
 #######
  alt_names <- setdiff(branch_name_vec, test_names)
  test_homRef <- get_samegeno_inds(geno_mat = geno_mat,
    branch_name_vec = test_names, genotype = 0, meta = meta, 
    n_miss = 0)
  test_homAlt <- get_samegeno_inds(geno_mat = geno_mat,
    branch_name_vec = test_names, genotype = 2, meta = meta,
    n_miss = 0)
  other_het <- get_samegeno_inds(geno_mat = geno_mat,
    branch_name_vec = alt_names, genotype = 1, meta = meta, n_miss = 0)
  #
  spec_homs <- intersect(union(test_homRef, test_homAlt), other_het)
  return(spec_homs)
}

test_shared_unique_het_combos <- function(geno_mat, branch_name_vec, n_test, 
  meta){
  # Find the number of shared, unique heterozygous SVs in all the possible
  #  combination of <n_test> samples
  ##############
  test_combos <- combn(seq(length(branch_name_vec)), m = n_test)
  n_diff_vec <- c()
  for(cc in seq(ncol(test_combos))){
    test_branches <- branch_name_vec[test_combos[ ,cc]]
    alt_branches <- setdiff(branch_name_vec, test_branches)
    test_hets <- get_samegeno_inds(geno_mat = geno_mat,
      branch_name_vec = test_branches, genotype = 1, meta = meta,
      n_miss = 0)
    alt_homRef <- get_samegeno_inds(geno_mat = geno_mat,
      branch_name_vec = alt_branches, genotype = 0, meta = meta,
      n_miss = 0)
    alt_homAlt <- get_samegeno_inds(geno_mat = geno_mat,
    branch_name_vec = alt_branches, genotype = 2, meta = meta,
      n_miss = 0)
  tmp_n_fixed_diff <- length(intersect(test_hets, 
    union(alt_homRef, alt_homAlt)))
  n_diff_vec <- c(n_diff_vec, tmp_n_fixed_diff)
  }
  combo_res_list <- list()
  combo_res_list[[1]] <- n_diff_vec
  combo_res_list[[2]] <- test_combos
  return(combo_res_list)
}

gen_mer_seqs <- function(mer_length, n_bp){
  # Generate sequences of stated length of repeated k_mers. For instance, if
  #   n_bp = 2, then generate binucleotide repeats of length <mer_length>
  #   Works for n_bp = 1 or n_bp = 2.
  # INPUTS #
  # mer_length = the total length of the binucleotide repeat sequences that
  #                are interested in
  # n_bp = the number of different nucleotides to have in the repeated
  #          sequence
  # OUTPUTS #
  # vector of repeated Kmer sequences
  ####
  nuc_combos <- combn(c('G', 'C', 'T', 'A'), m = n_bp, simplify = T)
  nuc_seqs <- apply(nuc_combos, 2, paste, sep = '', collapse = '')
  n_repeats <- mer_length / 2
  nuc_mer_seqs <- sapply(nuc_seqs, function(x)
    paste(rep(x, times = n_repeats), sep = '', collapse = ''))
  return(nuc_mer_seqs)
}

mer_counts <- function(sv_geno_df, nuc_seqs){
  # Count the number of repeat-containing (mono- or bi-nucleotide) kmers in 
  #   each SV sequence
  #   To be used for identifying suspect/bad/error-prone SVs
  # INPUTS # 
  # sv_geno_df = data.frame from VCF that includes the columns <ALT> which
  #                contains the sequence for tha ALTernate allele; <REF> which
  #                contains the sequence for the REFerence allele; and <type>
  #                which indicates whether the SV is INSertion or DELetion.
  #                The SV sequence for 'INS' is in the 'ALT' column, sequence
  #                for 'DEL' is in the 'REF' column
  # nuc_seqs = vector of repeat-containing sequences that want to count the
  #                instances of in the SV sequences
  # OUTPUT #
  # List; each element contains info for a sequence in <nuc_seqs>;
  #  sub-elements contain the count of that sequence in each 'INS' or 'DEL'
  #########
  library(stringr)
  nuc_counts <- list()
  for(bn in nuc_seqs){
    tmp_ins_count <- sapply(sv_geno_df$ALT[which(sv_geno_df$type == 'INS')],
      function(x) str_count(x, bn))
    tmp_del_count <- sapply(sv_geno_df$REF[which(sv_geno_df$type == 'DEL')],
      function(x) str_count(x, bn))
    nuc_counts[[bn]][['INS']] <- tmp_ins_count
    nuc_counts[[bn]][['DEL']] <- tmp_del_count
  }
  return(nuc_counts)
}

per_mer_length <- function(sv_geno_df, nuc_seqs){
  # Calculate the percentage of SV sequence that is comprised of mono- or 
  #   bi-nucleotide repeat sequences
  # INPUTS # 
  # sv_geno_df = data.frame from VCF that includes the columns <ALT> which
  #                contains the sequence for tha ALTernate allele; <REF> which
  #                contains the sequence for the REFerence allele; <type>
  #                which indicates whether the SV is INSertion or DELetion; and
  #                <sv_lenght> which contains the length of the SV.
  #                The SV sequence for 'INS' is in the 'ALT' column, sequence
  #                for 'DEL' is in the 'REF' column
  # nuc_seqs = vector of mono- or bi-nucleotide sequences that want to count 
  #                the instances of in the SV sequences; can be generated with
  #                gen_mer_seqs() function
  # OUTPUT #
  # List; each element contains info for a sequence in <nuc_seqs>;
  #   sub-elements contain the percentage of INS or DEL sequence that is
  #   comprised of the repeat sequence
  #######
  library(stringr)
  nuc_counts <- list()
  tmp_tot_ins_length <- sv_geno_df$sv_length[which(sv_geno_df$type == 'INS')]
  tmp_tot_del_length <- sv_geno_df$sv_length[which(sv_geno_df$type == 'DEL')]
  for(bn in nuc_seqs){
    seq_length <- nchar(bn) 
    tmp_ins_count <- sapply(sv_geno_df$ALT[which(sv_geno_df$type == 'INS')],
      function(x) str_count(x, bn))
    tmp_ins_bn_length <- tmp_ins_count * seq_length
    tmp_ins_per_bn <- tmp_ins_bn_length / tmp_tot_ins_length
    # 
    tmp_del_count <- sapply(sv_geno_df$REF[which(sv_geno_df$type == 'DEL')],
      function(x) str_count(x, bn))
    tmp_del_bn_length <- tmp_del_count * seq_length
    tmp_del_per_bn <- tmp_del_bn_length / tmp_tot_del_length
    #
    nuc_counts[[bn]][['INS']] <- tmp_ins_per_bn
    nuc_counts[[bn]][['DEL']] <- tmp_del_per_bn
  }
  return(nuc_counts)
}

per_mer_inds <- function(sv_geno_df, per_mer_list, per_cutoff){
  # Get the indices of SVs that consist of a percentage of mono- or 
  #   bi-nucleotide repeat
  #   sequences that are above <per_cutoff>
  # INPUTS #
  # sv_geno_df = data.frame from VCF that includes the column <sv_lenght> 
  #                which contains the length of the SV. Preferable same 
  #                genotype dataframe used to generate <per_mer_list>
  # per_mer_list = list containing, for multile repeat sequences,
  #                  the percentage of each SV sequence consisting of each
  #                  repeated sequence. Should be generated using the 
  #                  per_mer_length() function
  # per_cutoff = the percentage cutoff. If a percentage of an SV's sequence 
  #                that is a repeat sequence, then it's index
  #                will be retained
  # OUTPUTS #
  # List; each element contains info for a repeat sequence;
  #   sub-elements contain the indices of INS or DEL sequences that contain
  #   a specific repeat sequence above the per_cutoff threshold
  #######
  per_inds <- list()
  tot_ins_inds <- which(sv_geno_df$type == 'INS')
  tot_del_inds <- which(sv_geno_df$type == 'DEL')
  for(bnp in names(per_mer_list)){
    tmp_ins_inds <- tot_ins_inds[
                      which(per_mer_list[[bnp]][['INS']] > per_cutoff)]
    tmp_del_inds <- tot_del_inds[
                      which(per_mer_list[[bnp]][['DEL']] > per_cutoff)]
    per_inds[[bnp]][['INS']] <- tmp_ins_inds
    per_inds[[bnp]][['DEL']] <- tmp_del_inds
  }
  return(per_inds)
}

id_prob_SV_seqs <- function(sv_geno_df, mer_length = 8, per_mn_cutoff = 0.7, 
  per_bn_pure_cutoff = 0.5, per_bn_multi_cutoff = 0.6){
  # Identify SV's that have problematic sequences in terms of their content
  #  of binucleotide repeat sequences. These sequences are associated with
  #  sequenceing errors and are highly error prone.
  # sv_geno_df = data.frame from VCF that includes the columns <ALT> which
  #                contains the sequence for tha ALTernate allele; <REF> which
  #                contains the sequence for the REFerence allele; <type>
  #                which indicates whether the SV is INSertion or DELetion; and
  #                <sv_lenght> which contains the length of the SV.
  #                The SV sequence for 'INS' is in the 'ALT' column, sequence
  #                for 'DEL' is in the 'REF' column
  # mer_length = the total length of the binucleotide repeat sequences that
  #                are interested in
  # per_mn_cutoff = the percentage cutoff for a mononucleotide repeat. 
  #                     If the percentage of an SV's sequence that is the 
  #                     repeat sequence, then it's index will be retained
  # per_bn_pure_cutoff = the percentage cutoff for a single binucleotide 
  #                     repeat. If the percentage of an SV's sequence that is 
  #                     a single type of  binucleotide repeat sequence, 
  #                     then it's index will be retained
  # per_bn_multi_cutoff = the percentage cutoff for two type of binucleotide
  #                      repeat. If two binucleotide repeate sequences make up
  #                      a percentage of an SV sequence above this cutoff, then
  #                      the index is retained for that SV
  # OUTPUTS #
  # Vector of SV in sv_geno_df that are problematic
  #####
  mn_mer_seqs <- gen_mer_seqs(mer_length = 8, n_bp = 1)
  bn_mer_seqs <- gen_mer_seqs(mer_length = 8, n_bp = 2)
  #
  mn_per_mer <- per_mer_length(sv_geno_df = sv_geno_df, nuc_seqs = mn_mer_seqs)
  mn_per_inds <- unlist(per_mer_inds(sv_geno_df = sv_geno_df, 
    per_mer_list = mn_per_mer, per_cutoff = per_mn_cutoff))
  #
  bn_per_mer <- per_mer_length(sv_geno_df = sv_geno_df, nuc_seqs = bn_mer_seqs)
  bn_per_pure_inds <- unlist(per_mer_inds(sv_geno_df = sv_geno_df, 
    per_mer_list = bn_per_mer, per_cutoff = per_bn_pure_cutoff))
  #
  mult_bn_cut <- per_bn_multi_cutoff / 2
  tmp_multi_inds <- unlist(per_mer_inds(sv_geno_df = sv_geno_df,
    per_mer_list = bn_per_mer, per_cutoff = mult_bn_cut))
  tmp_mult_table <- table(unlist(tmp_multi_inds))
  bn_per_multi_inds <- as.numeric(names(tmp_mult_table))[
                         which(tmp_mult_table > 1)]
  #
  tot_rem_inds <- sort(union(mn_per_inds, 
    union(bn_per_pure_inds, bn_per_multi_inds)))
  return(tot_rem_inds)
}

dup_inds_to_remove <- function(sv_geno_df){
  # For duplicated SV positions, get the indices of the smaller SVs that will
  #   be removed
  # INPUTS #
  # sv_geno_df = data.frame from VCF that includes the columns <full_name> 
  #                which is character string of Chromosome and positon pasted
  #                together; and <sv_lenght> which contains the length of 
  #                the SV.
  # OUTPUTS #
  # vector containing the indices of sv_geno_df of duplicated SVs to be 
  #  removed; these are SVs that are in duplicated positions and that are 
  #  smaller than the largets SV at that position
  ##########
  name_table <- table(sv_geno_df$full_name)
  dup_names <- names(name_table)[which(name_table > 1)]
  rm_inds <- c()
  for(dn in dup_names){
    tmp_inds <- which(sv_geno_df$full_name == dn)
    sv_lengths <- sv_geno_df$sv_length[tmp_inds]
    keep_ind <- tmp_inds[which.max(sv_lengths)]
    rm_inds <- c(rm_inds, setdiff(tmp_inds, keep_ind))
  }
  return(rm_inds)
}

overlap_inds_to_remove <- function(sv_geno_df, dist_cut){
  # Find SVs that are close and/or overlap and get the indices of the smaller
  #   SVs that should be removed
  #######3
  chrom_vec <- unique(sv_geno_df$CHROM)
  rm_ind_vec <- c()
  for(cv in chrom_vec){
    tmp_chrom_inds <- which(sv_geno_df$CHROM == cv)
    tmp_sub_df <- sv_geno_df[tmp_chrom_inds,]
    for(i in seq(nrow(tmp_sub_df))){
      tmp_diffs <- abs(tmp_sub_df$POS[i] - tmp_sub_df$POS)
      close_inds <- which(tmp_diffs <= dist_cut)
      tmp_sv_len <- tmp_sub_df$sv_length[close_inds]
      keep_ind <- close_inds[which.max(tmp_sv_len)]
      tmp_rm_inds <- setdiff(close_inds, keep_ind)
      rm_ind_vec <- c(rm_ind_vec, tmp_chrom_inds[tmp_rm_inds])
    }
  }
  return(unique(rm_ind_vec))
}



