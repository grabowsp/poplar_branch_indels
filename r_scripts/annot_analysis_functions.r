# Functions used for analysis of annotation of SVs

gene_in_SV <- function(annot_df, chr_name, sv_start, sv_length){
  # Function to find full genes (or annotation entries) withing the boundary
  #   of an SV
  # INPUTS
  # annot_df = annotation dataframe: col_1 = chromosome name, col_2 = source 
  #              info, col_3 = annotation description, col_4 = start pos,
  #              col_5 = end pos, then additional columns 
  # chr_name = the name of the chromosome of the SV
  # sv_start = starting position of the SV
  # sv_length = length of SV
  # OUTPUT
  # index/indices of genes in annot_df that are completely within the SV
  ###########
  chr_inds <- which(annot_df[, 1] == chr_name)
  st_inds <- intersect(which(annot_df[, 4] >= sv_start), chr_inds)
  sv_end <- sv_start + sv_length
  en_inds <- intersect(which(annot_df[, 5] <= sv_end), chr_inds)
  overlap_inds <- intersect(st_inds, en_inds)
  return(overlap_inds)
}

gene_5overlap_SV <- function(annot_df, chr_name, sv_start){
  # Function to find genes (annotation entries) that overlap with the 
  #  5'/upstream end of a SV but are not completely within the SV. 
  # INPUTS
  # annot_df = annotation dataframe: col_1 = chromosome name, col_2 = source 
  #              info, col_3 = annotation description, col_4 = start pos,
  #              col_5 = end pos, then additional columns 
  # chr_name = the name of the chromosome of the SV
  # sv_start = starting position of the SV
  # OUTPUT
  # index of gene(s) in annot_df that overlap with 5' end of SV
  ########
  chr_inds <- which(annot_df[, 1] == chr_name)
  st_5_inds <- intersect(which(annot_df[,4] < sv_start), chr_inds)
  en_5_inds <- intersect(which(annot_df[,5] >= sv_start), chr_inds)
  overlap_inds <- intersect(st_5_inds, en_5_inds)
  return(overlap_inds)
}

gene_3overlap_SV <- function(annot_df, chr_name, sv_start, sv_length){
  # Function to find genes (annotation entries) that overlap with the 
  #  3'/upstream end of a SV but are not completely within the SV.
  # INPUTS
  # annot_df = annotation dataframe: col_1 = chromosome name, col_2 = source 
  #              info, col_3 = annotation description, col_4 = start pos,
  #              col_5 = end pos, then additional columns 
  # chr_name = the name of the chromosome of the SV
  # sv_start = starting position of the SV
  # sv_length = length of SV
  # OUTPUT
  # index of gene(s) in annot_df that overlap with 3' end of SV
  ######
  chr_inds <- which(annot_df[, 1] == chr_name)
  sv_end <- sv_start + sv_length
  st_3_inds <- intersect(which(annot_df[,4] <= sv_end), chr_inds)
  en_3_inds <- intersect(which(annot_df[,5] > sv_end), chr_inds)
  overlap_inds <- intersect(st_3_inds, en_3_inds)
  return(overlap_inds)
}

get_unique_code_seq <- function(annot_df, target_inds){
  # Function to get the unique genic positions for a series of genes or 
  #   other annotation entries in annot_df
  # INPUTS
  # annot_df = annotation dataframe: col_1 = chromosome name, col_2 = source 
  #              info, col_3 = annotation description, col_4 = start pos,
  #              col_5 = end pos, then additional columns
  # target_inds = the row indices of annot_df to be included for getting
  #                  the positions; these should all be on the same chromosome
  # OUTPUT 
  # vector of all the positions in part of the annotated entries of the
  #  selected indices
  #######
  tmp_code_vec <- c()
  for(ti in target_inds){
    tmp_p <- c(annot_df[ti,4]:annot_df[ti,5])
    tmp_code_vec <- union(tmp_code_vec, tmp_p)
  }
  return(tmp_code_vec)
}

get_perc_code_seq <- function(sv_start, sv_length, unique_code_seq){
  # Function for calculating the percentage of the SV that contains
  #   conding sequence (or whatever seq is inputted)
  # INPUTS
  # sv_start = starting position of the SV
  # sv_length = length of SV
  # unique_code_seq = vector of positions that are coding sequence (or
  #                     whatever type of sequence is of interest); is best if
  #                     generated by get_unique_code_seq function
  # OUTPUT
  # percentage of sequence in SV that is genic (or whatever the annotation is)
  ############
  sv_seq <- c(sv_start:(sv_start+sv_length-1))
  sv_with_code_seq <- intersect(sv_seq, unique_code_seq)
  per_code_seq <- length(sv_with_code_seq)/sv_length
  return(per_code_seq)
}

calc_sv_perc_code <- function(sv_name, annot_df){
  # Wrapper for calculating the percentage of an SV that contains coding
  #   (or whatever seq interested in) seq
  # INPUTS
  # sv_name = name of the SV; generally in the form CHR_POS_TYPE_SIZE
    # annot_df = annotation dataframe: col_1 = chromosome name, col_2 = source 
  #              info, col_3 = annotation description, col_4 = start pos,
  #              col_5 = end pos, then additional columns
  # OUTPUT
  # percentage of sequence in SV that is genic (or whatever the annotation is) 
  ############
  name_split <- unlist(strsplit(sv_name, split = '_'))
  if(name_split[1] == 'scaffold'){
    tmp_name <- paste(name_split[1], name_split[1], sep = '_')
    name_split <- c(tmp_name, name_split[3:5])
  }
  tmp_chr <- name_split[1]
  sv_start_pos <- as.numeric(name_split[2])
  sv_length_num <- as.numeric(name_split[4])
  annot_inds <- c(
    gene_in_SV(annot_df = annot_df, chr_name = tmp_chr, 
      sv_start = sv_start_pos, sv_length = sv_length_num),
    gene_5overlap_SV(annot_df = annot_df, chr_name = tmp_chr, 
      sv_start = sv_start_pos),
    gene_3overlap_SV(annot_df = annot_df, chr_name = tmp_chr,
      sv_start = sv_start_pos, sv_length = sv_length_num)
  )
  tmp_unique_seq <- get_unique_code_seq(annot_df = annot_df, 
    target_inds = sort(annot_inds))
  perc_code_seq <- get_perc_code_seq(sv_start = sv_start_pos, 
    sv_length = sv_length_num, unique_code_seq = tmp_unique_seq)
  return(perc_code_seq)
}

gen_window_inds <- function(max_pos, window_size, end_add = 2000, min_pos = 1){
  # Generate indices for window sizes
  #  note: currently, does not generate final window at end of chromosome
  # INPUTS
  # max_pos = the largest known position on the chromosome
  # window_size = the size of non-overlapping windows
  # end_add = the amount to add to max_pos to estimate the amount of seq
  #             at the end of the chromosome after known max_pos
  # min_pos = the minimum position on the chromosome
  # OUTPUT
  # data.frame: 1st column is beginning position of window, 2nd column is
  #  ending position of window
  ##########  
  full_max <- (max_pos + end_add)
  tmp_start <- seq(from = min_pos, to = full_max, by = window_size)
  tmp_end <- tmp_start + window_size - 1
  window_df <- data.frame(start = tmp_start, end = tmp_end, 
    stringsAsFactors = F)
  return(window_df)
}

get_in_annot_inds <- function(annot_df, chr_name, start_pos, end_pos){
  # Function to get all gene indices within an interval, including 5' and
  #   3' overlaps
  # INPUTS
  # annot_df = annotation dataframe: col_1 = chromosome name, col_2 = source 
  #              info, col_3 = annotation description, col_4 = start pos,
  #              col_5 = end pos, then additional columns
  # chr_name = name of chromosome with the interval
  # start_pos = starting position of the interval
  # end_pos = ending position of the interval
  # OUTPUT
  # indices of annotation items within the interval
  ###########
  tmp_length <- end_pos - start_pos
  annot_inds <- c(
    gene_in_SV(annot_df = annot_df, chr_name = chr_name,
      sv_start = start_pos, sv_length = tmp_length),
    gene_5overlap_SV(annot_df = annot_df, chr_name = chr_name,
      sv_start = start_pos),
    gene_3overlap_SV(annot_df = annot_df, chr_name = chr_name,
      sv_start = start_pos, sv_length = tmp_length)
  )
  inds_sorted <- sort(annot_inds)
  return(inds_sorted)
}

calc_window_perc_code <- function(annot_df, chr_name, start_pos, end_pos){
  # Calculate percentage of sequence in a window/interval that is coding (or
  #  or whatever annotation interested in)
  # INPUTS
  # annot_df = annotation dataframe: col_1 = chromosome name, col_2 = source 
  #              info, col_3 = annotation description, col_4 = start pos,
  #              col_5 = end pos, then additional columns
  # chr_name = name of chromosome with the interval
  # start_pos = starting position of the interval
  # end_pos = ending position of the interval
  # OUTPUT
  # percentage of seq in window/interval that is part of that annotation
  ################3
  tmp_annot_inds <- get_in_annot_inds(annot_df = annot_df, chr_name = chr_name,
    start_pos = start_pos, end_pos = end_pos)
  tmp_unique_seq <- get_unique_code_seq(annot_df = annot_df,
    target_inds = tmp_annot_inds)
  tmp_length <- end_pos - start_pos
  perc_code_seq <- get_perc_code_seq(sv_start = start_pos,
    sv_length = tmp_length, unique_code_seq = tmp_unique_seq)
  return(perc_code_seq)
}

calc_genomewide_perc <- function(annot_df, window_size, loud = F){
  ##########
  chr_names <- unique(annot_df[,1])
  chr_names <- chr_names[-grep('scaffold', chr_names)]
  out_wind_list <- list()
  for(pchr in chr_names){
    if(loud){print(pchr)}
    chr_inds <- which(annot_df[,1] == pchr)
    tmp_max_pos <- max(annot_df[chr_inds, 5])
    tmp_window_df <- gen_window_inds(max_pos = tmp_max_pos, 
      window_size = window_size)
    tmp_wind_perc <- apply(tmp_window_df, 1, function(x)
      calc_window_perc_code(annot_df = annot_df, chr_name = pchr,
        start_pos = x[1], end_pos = x[2]))
    out_wind_list[[pchr]] <- tmp_wind_perc
  }
  return(out_wind_list)
}

