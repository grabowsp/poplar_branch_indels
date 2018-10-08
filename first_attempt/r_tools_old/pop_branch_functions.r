# R functions for analyzing poplar branch structural variation

proc_share_df <- function(share_df, indel_type){
  ##################################
  # Process <NAME.diff.sites_in_files> outputted by vcftools and split up 
  #   by insertions and deletions 
  #  This file contains info about shared and different variants between two
  #   vcf files
  # INPUTS
  # share_df = data.frame containing data imported from 
  #              <NAME.diff.sites_in_files>
  # indel_type = either "DEL" or "INS" depending on whether imported file
  #               contains data from insertions or deletions
  # OUTPUT
  # data.frame with columns showing indel size and type
  #################################### 
  sites = share_df
  sites$POS1 <- as.integer(sites$POS1)
  sites$POS2 <- as.integer(sites$POS2)
  indel_size_1 <- unlist( lapply(sites$REF1, function(x) 
                    length( unlist( strsplit(x, split = '') ) ) ) )
  indel_size_2 <- unlist( lapply(sites$REF2, function(x) 
                    length( unlist( strsplit(x, split = '') ) ) ) )
  indel_size_3 <- unlist( lapply(sites$ALT1, function(x) 
                    length( unlist( strsplit(x, split = '') ) ) ) )
  indel_size_4 <- unlist( lapply(sites$ALT2, function(x) 
                    length( unlist( strsplit(x, split = '') ) ) ) )
  indel_size_max <- apply(cbind(indel_size_1, indel_size_2, indel_size_3, 
                      indel_size_4), 1, max, na.rm = F)
  sites$INDEL_SIZE <- indel_size_max
  sites$INDEL_TYPE <- indel_type
  return(sites)
}

dist_to_other_vars <- function(share_df){
  ########################
  # Calculate distance to closest indel in either sample
  # INPUTS
  # share_df = data.frame that contains both INS and DEL data from a comparison
  #              of two samples
  # OUTPUT
  # data.frame with distances to variants in the different samples as well
  #  as the shortest distance
  ##########################
  share_df$close_altSamp_sameVar <- NA
  share_df$close_sameSamp_sameVar <- NA
  share_df$close_altSamp_altVar <- NA
  share_df$close_sameSamp_altVar <- NA
  #
  vartype_inds <- list()
  vartype_inds[['INS']] <- which(share_df$INDEL_TYPE == 'INS')
  vartype_inds[['DEL']] <- which(share_df$INDEL_TYPE == 'DEL')
  #
  for(chr in unique(share_df$CHROM)){
    chr_inds <- which(share_df$CHROM == chr)
    # print(chr) #COMMENT
    for(j in c('1', '2')){
      samp_inds <- intersect(chr_inds, which(share_df$IN_FILE == j))
      samp_col <- paste('POS', j, sep = '')
      other_col <- paste('POS', setdiff(c('1', '2'), j), sep = '')
      # print(j) # COMMENT
      for(k in seq(2)){
        sameVar_inds <- intersect(chr_inds, vartype_inds[[k]])
        altVar_inds <- intersect(chr_inds, vartype_inds[[setdiff(seq(2),k)]])
        #
        altSamp_sameVar_data <- share_df[sameVar_inds, other_col]
        altSamp_altVar_data <- share_df[altVar_inds, other_col]
        sameSamp_sameVar_data <- share_df[sameVar_inds, samp_col]
        sameSamp_altVar_data <- share_df[altVar_inds, samp_col]
        #
        s_t_inds <- intersect(vartype_inds[[k]], samp_inds)
        # print(k) #COMMENT
        for(i in s_t_inds){
          samp_pos <- share_df[i, samp_col]
          #
          share_df$close_altSamp_sameVar[i] <- min(abs(
            altSamp_sameVar_data - samp_pos), na.rm = T)
          share_df$close_altSamp_altVar[i] <- min(abs(
            altSamp_altVar_data - samp_pos), na.rm = T)
          share_df$close_sameSamp_sameVar[i] <- min(abs
            (setdiff(sameSamp_sameVar_data - samp_pos, 0)), na.rm = T)
          share_df$close_sameSamp_altVar[i] <- min(abs
            (sameSamp_altVar_data - samp_pos), na.rm = T)
        }
      }
    }
  }
  return(share_df)
}

sampnames_in_compfile <- function(file, file_dir, file_prefix, file_suffix,
    split_chr = '_v_'){
  # Extract the sample names from the file name of indel comparison files
  #  that end with .diff_indels. May work with other files, but written
  #  specifically for the .diff_indels files.
  ##
  # INPUTS
  # file = file that is being analyzed - should following format:
  #          FILE_DIR/FILE_PREFIX<samp1_v_samp2>FILE_SUFFIX
  # file_dir = the directory where file is found, ending with '/'
  # file_prefix = part of the filename that preceeds the sample names
  #   ex: for <popbranch_13_1_v_13_2.diff_indels>, file prefix = 'popbranch_'
  # file_suffix = part of the filename that follows the second sample
  #   ex: for <popbranch_13_1_v_13_2.diff_indels>, file_suffix = '.diff_indels'
  # split_chr = character string that separates the samples in the filename
  #   ex: for <popbranch_13_1_v_13_2.diff_indels>, split_chr = '_v_'
  # OUTPUT
  # vector with the sample names in sampe order as in the filename
  #######################
  samps_in_file <- gsub(file_suffix, '',
    gsub(paste(file_dir, file_prefix, sep = ''), '', file), fixed = T)
  names_in_file <- unlist(strsplit(samps_in_file, split = split_chr))
  return(names_in_file)
}

proc_diffindel_data <- function(file, test_samp, names_in_file){
  # Read <FILE.diff_indels> file and process results specific to "test_samp"
  #  To be used for generating files for each sample about shared and different
  #  indels.
  # INPUTS
  # file = the .diff_indels file
  # test_samp = the sample name that are processing the file for
  # names_in_file = vector of names in filename of <file>, as generated by 
  #   sampnames_in_compfile() function. <test_samp> should be one of the
  #   samples in <names_in_file>
  # OUTPUT
  # data.frame only containing indels found in test_samp, with a new "snp_name"
  #   column, with data.frame sorted by order in "snp_name" column
  ######################
  indel_info <- read.table(file, header = T, stringsAsFactors = F, sep = '\t')
  # figure out if test samp is "1" or "2" in data.frame
  test_samp_numb <- as.character(which(names_in_file == test_samp))
  test_s_inds <- union(which(indel_info$IN_FILE == test_samp_numb), 
     union(which(indel_info$IN_FILE == 'B'), which(indel_info$IN_FILE == 'O')))
  info_sub <- indel_info[test_s_inds, ]
  # add names for SNPs that can be sorted
  info_sub$snp_name <- paste(info_sub$CHROM, info_sub[, 
    paste('POS', test_samp_numb, sep = '')], sep = '_')
  info_sort <- info_sub[order(info_sub$snp_name), ]
  return(info_sort)
}

gen_samp_tot_diff_tab <- function(test_samp, file_dir, file_prefix, 
  file_suffix, split_chr = '_v_'){
  # Generate table that contains shortest distance to closest indel for each 
  ########
  file_names <- system(paste('ls ', file_dir, '*', file_suffix, sep = ''), 
    intern = T)
  # pull out the names of all the samples in the directory
  samp_comps <- gsub(file_suffix, '',
    gsub(paste(file_dir, file_prefix, sep = ''), '', file_names), fixed = T)
  all_samps <- sort(unique(unlist(strsplit(samp_comps, split = split_chr))))
  #
  samp_files <- file_names[grep(test_samp, file_names)]
  #use data from first file to generate the overall table which will be 
  #  filled-in with the rest of the files
  file_ind <- 1
  tmp_samps <- sampnames_in_compfile(samp_files[file_ind], file_dir = file_dir,
    file_prefix = file_prefix, file_suffix = file_suffix, split_chr = split_chr)
  test_info <- proc_diffindel_data(file = samp_files[file_ind], 
    test_samp = test_samp, names_in_file = tmp_samps)
  test_samp_numb <- as.character(which(tmp_samps == test_samp))
  #
  samp_tot_info <- data.frame(test_info[, c('snp_name', 'CHROM', 
    paste('POS', test_samp_numb, sep = ''), 'INDEL_TYPE', 'INDEL_SIZE')], 
    stringsAsFactors = F)
  # add column for each of the comparison samples
  samp_tot_info[, setdiff(all_samps, test_samp)] <- NA
  # add info from first file
  comp_col <- which(colnames(samp_tot_info) == setdiff(tmp_samps, test_samp))
  samp_tot_info[ , comp_col] <- test_info$min_var_dist
  # If same indel found at same position, 'B', then distance is 0
  samp_tot_info[ which(test_info$IN_FILE == 'B'), comp_col] <- 0
  # If same or similar indel is overlapping, 'O', then distance is 1
  samp_tot_info[which(test_info$IN_FILE == 'O'), comp_col] <- 1
# add data from remaining files
  for(i in c(2:length(samp_files))){
    file_ind <- i
    tmp_samps_1 <- sampnames_in_compfile(samp_files[file_ind], 
      file_dir = file_dir, file_prefix = file_prefix, 
      file_suffix = file_suffix, split_chr = split_chr)
    test_info_2 <- proc_diffindel_data(file = samp_files[file_ind], 
      test_samp = test_samp, names_in_file = tmp_samps_1)
    comp_col <- which(colnames(samp_tot_info) == 
      setdiff(tmp_samps_1, test_samp))

    if(sum(samp_tot_info$snp_name == test_info_2$snp_name) != 
     nrow(test_info_2)){
       print(paste('SNP NAMES DONT MATCH BETWEEN TOTAL TABLE and FILE', i))
    }
    samp_tot_info[ , comp_col] <- test_info_2$min_var_dist
    samp_tot_info[ which(test_info_2$IN_FILE == 'B'), comp_col] <- 0
    samp_tot_info[which(test_info_2$IN_FILE == 'O'), comp_col] <- 1
  }
  return(samp_tot_info)
}
