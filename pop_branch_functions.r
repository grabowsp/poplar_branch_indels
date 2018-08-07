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


