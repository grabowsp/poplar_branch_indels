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


