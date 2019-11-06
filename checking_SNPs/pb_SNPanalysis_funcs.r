# Functions for analyzing SNP data using PacBio results

# Sections:
# 1) Functions for analyzing genotype probability files outputted by Sujan
# 2) Functions for analyzing the mpileup results from samtools
# 3) Functions for subsampling and generating genotypes from Illumina VCFs 

# Section 1) Functions for analyzing the genotype probability files 
#   outputted by Sujan 

gen_geno_probs <- function(geno_df){
  ## Generate the relative probabilities of each genotype in each sample
  ##   using the different genotypes' conditional probabilities
  # INPUTS
  # geno_df = data.frame that includes the conditional probabilities for 
  #             each genotype for each sample at each SNP
  # OUTPUT
  # List with three elements, each element is a matrix of relative 
  #  probabilities of a genotype class for each sample, col = samp, row = SNP;
  #  [[1]] = homozygous REF (RR); [[2]] = HET (H); [[3]] = homozygous ALT (VV)
  ########
  geno_mat_list <- list()
  branch_cols <- c(5:(ncol(geno_df)-1))
  RR_mat <- H_mat <- VV_mat <- matrix(NA, ncol = length(branch_cols),
    nrow = nrow(geno_df))
  for(bc in seq(length(branch_cols))){
    tmp_probs <- strsplit(geno_df[ , branch_cols[bc]], split = ':')
    tmp_probs <- lapply(tmp_probs, function(x) as.numeric(x))
    tmp_geno_probs <- lapply(tmp_probs, function(x) x/sum(x))
    RR_mat[,bc] <- unlist(lapply(tmp_geno_probs, function(x) x[1]))
    H_mat[,bc] <- unlist(lapply(tmp_geno_probs, function(x) x[2]))
    VV_mat[,bc] <- unlist(lapply(tmp_geno_probs, function(x) x[3]))
  }
  geno_mat_list[['RR_mat']] <- RR_mat
  geno_mat_list[['H_mat']] <- H_mat
  geno_mat_list[['VV_mat']] <- VV_mat
  return(geno_mat_list)
}
#########
get_good_snps <- function(geno_df, good_cut){
  # Get the indicices of SNPs that have at least one "good" HOM and HET
  #  genotype
  # geno_df = data.frame that includes the conditional probabilities for 
  #             each genotype for each sample at each SNP
  # good_cut = the relative probability cutoff necessary for calling a "good"
  #              genotype (ex: 0.999)
  # OUTPUT
  # vector of indices of rows from geno_df that pass the test
  ################
  tmp_geno_probs <- gen_geno_probs(geno_df = geno_df)
  #
  good_RR_snps <- which(apply(tmp_geno_probs[[1]], 1, max) > good_cut)
  good_VV_snps <- which(apply(tmp_geno_probs[[3]], 1, max) > good_cut)
  good_HOM_snps <- sort(union(good_RR_snps, good_VV_snps))
  good_HET_snps <- which(apply(tmp_geno_probs[[2]], 1, max) > good_cut)
  #
  good_both_snps <- intersect(good_HOM_snps, good_HET_snps)
  return(good_both_snps)
}
###########
gen_geno_dosages <- function(geno_df){
  # Generate the genotype dosage based on the relative probability
  #  RR = 0, RA = 1, AA = 2
  # INPUTS
  # geno_df = data.frame that includes the conditional probabilities for 
  #             each genotype for each sample at each SNP
  # OUTPUT = matrix of the estimated genotype for each sample at each 
  # SNP
  ########
  branch_cols <- c(5:(ncol(geno_df)-1))
  split_list <- apply(geno_df[, branch_cols], 1, function(x) 
    strsplit(x, split = ':'))
  tmp_probs <- lapply(split_list, function(x)
    lapply(x, function(y) as.numeric(y)/ sum(as.numeric(y))))
  tmp_dose <- lapply(tmp_probs, function(x) 
    round(unlist(lapply(x, function(y) sum(y * c(0,1,2)))) , digits = 3))
  dose_mat <- matrix(data = unlist(tmp_dose), byrow = T, 
    nrow = length(tmp_dose))
  return(dose_mat)
}
############
assign_geno_from_dosage <- function(dosage_mat, dose_dist_cut){
  # Assign at R,H,V, or N genotype based on the dosage; dosages more than
  #  dose_dist_cut from 0,1,or2 will be assigned N
  # INPUTS
  # dosage_mat = matrix of the genotype dosages; generated by 
  #                gen_geno_dosages function
  # dose_dist_cut = the distance allowed from pure dosage to be considered
  #                   a certain genotype
  # OUTPUT
  #
  #########
  tmp_geno_mat <- matrix('N', nrow = nrow(dosage_mat), ncol = ncol(dosage_mat))
  r_inds <- which(dosage_mat <= (0+dose_dist_cut))
  h_inds <- which((dosage_mat >= (1-dose_dist_cut)) & 
    (dosage_mat <= (1+dose_dist_cut)))
  v_inds <- which(dosage_mat >= (2-dose_dist_cut))
  tmp_geno_mat[r_inds] <- 'R'
  tmp_geno_mat[h_inds] <- 'H'
  tmp_geno_mat[v_inds] <- 'V'
  return(tmp_geno_mat)
}
#############
#############

# Seciton 2) Functions for analyzing the mpileup results from samtools

make_count_mat <- function(counts_df, nread){
  # Function to count up the number of alleles in a SNP above a count
  #  threshold
  # INPUTS
  # counts_df = data.frame of the counts of the 4 possible alleles
  # nread = the minimum readcount threashold for an allele to be counted
  # OUTPUT
  # matrix of the number of alleles above the readcount threshold for each
  #  sample at each SNP
  ############
  out_count_mat <-  matrix(data = as.numeric(NA), nrow = nrow(counts_df),
    ncol = ncol(counts_df))
  for(i in seq(ncol(counts_df))){
    tmp_count_list <- strsplit(counts_df[,i], split = ',')
    tmp_count_vec <- unlist(lapply(tmp_count_list, function(x)
      sum(as.numeric(x) > nread)))
    out_count_mat[,i] <- tmp_count_vec
  }
  return(out_count_mat)
}
###############
make_per_mat <- function(counts_df, per_cut){
  # Function to count up the number of alleles in a SNP above a percentage
  #  of all reads at the SNP
  # INPUTS
  # counts_df = data.frame of the counts of the 4 possible alleles
  # per_cut = the minimum readcount percentage for an allele to be counted
  # OUTPUT
  # matrix of the number of alleles above the readcount percentage 
  #  threshold for each sample at each SNP
  ############### 
  out_count_mat <-  matrix(data = as.numeric(NA), nrow = nrow(counts_df),
    ncol = ncol(counts_df))
  for(i in seq(ncol(counts_df))){
    tmp_count_list <- strsplit(counts_df[,i], split = ',')
    tmp_count_vec <- unlist(lapply(tmp_count_list, function(x)
      sum((as.numeric(x)/sum(as.numeric(x))) > per_cut)))
    out_count_mat[,i] <- tmp_count_vec
  }
  return(out_count_mat)
}
##############
make_tot_depth_mat <- function(counts_df){
  # Calculate total depth for SNP
  # INPUTS
  # counts_df = data.frame of the counts of the 4 possible alleles
  # OUTPUTS
  # matrix of the total read-depth for each sample at each SNP
  ##########
  out_depth_mat <- matrix(data = as.numeric(NA), nrow = nrow(counts_df),
    ncol = ncol(counts_df))
  for(i in seq(ncol(counts_df))){
    tmp_count_list <- strsplit(counts_df[,i], split = ',')
    tmp_count_vec <- unlist(lapply(tmp_count_list, function(x)
      sum(as.numeric(x))))
    out_depth_mat[,i] <- tmp_count_vec
  }
  return(out_depth_mat)
}
##############
##############

# Section 3) Functions for subsampling and generating genotypes from 
#    Illumina VCFs and PacBio Count data

get_vcf_allele_counts <- function(vcf_geno_vec){
  # Function for extracting the allele counts for a sample at each SNP in 
  #   a VCF
  # INPUT
  # vcf_geno_vec = vector of vcf genotypes with following format:
  #   Genotype:Seq_Depth:Ref_depth,Alt_depth. Ex: 0/1:40:20,10
  # OUTPUT
  # List of allele counts
  ##########
  tmp_counts_1 <- sapply(vcf_geno_vec, function(x)
    unlist(strsplit(x, split = ':'))[3])
  tmp_counts_list_0 <- strsplit(tmp_counts_1, split = ',')
  tmp_counts_list <- lapply(tmp_counts_list_0, as.numeric)
  return(tmp_counts_list)
}
# test_counts_list <- get_vcf_allele_counts(t14_Illum_vcf[,10])
###################
generate_subsamp_counts <- function(allele_count_list, n_reads){
  # Function to generate a subsampled allele counts based on the
  #   overall allele counts
  # INPUTS
  # allele_count_list = list of allele counts in c(R,A) format
  # n_reads = the number of reads to subsample
  # OUTPUT
  # List of subsampled allele counts
  ###############
  full_readcount_list <- lapply(allele_count_list, function(x)
    c(rep('R', times = as.numeric(x)[1]), rep('A', times = as.numeric(x[2])))
  )
  subsamp_readcount_list <- lapply(full_readcount_list, function(x)
    if(length(x) < n_reads){return(NA)} else{sample(x, size = n_reads)}
  )
  subsamp_count_list <- lapply(subsamp_readcount_list, function(x)
    c(sum(x == 'R'), sum(x == 'A')))
  return(subsamp_count_list)
}
# test_subsamp_counts <- generate_subsamp_counts(test_counts_list, n_reads = 30)
###################
generate_genotype_dosages <- function(allele_count_list, seq_err = 0.01){
  # Function to generate dosage genotypes using conditional probabilities of
  #  the 3 types of genotypes and the allele counts
  #  Works for full allele counts or subsampled counts
  # INPUTS
  # allele_count_list = list of allele counts in c(R,A) format
  # seq_err = the estimated sequencing error rate used for calculating
  #            conditional probabilities of the genotypes
  # OUTPUT
  # vector of dosage genotypes: 0 = RR; 1 = RA; 2 = AA
  #######
  hom_ref_prob_list <- lapply(allele_count_list, function(x)
    dbinom(x[2], sum(x), prob = seq_err))
  het_prob_list <- lapply(allele_count_list, function(x)
    dbinom(x[2], sum(x), prob = 0.5))
  hom_alt_prob_list <- lapply(allele_count_list, function(x)
    dbinom(x[1], sum(x), prob = seq_err))
  #
  prob_df <- data.frame(p_RR = unlist(hom_ref_prob_list),
    p_RA = unlist(het_prob_list), p_AA = unlist(hom_alt_prob_list),
    stringsAsFactors = F)
  #
  dose_vec <- apply(prob_df, 1, function(x) (1*x[2] + 2*x[3])/sum(x))
  return(dose_vec)
}
# test_dose_vec <- generate_genotype_dosages(test_subsamp_counts)
##################
generate_subsamp_dosages <- function(vcf_geno_vec, n_reads, seq_err = 0.01){
  # Function to go from vcf genotypes to subsampled dosages
  # INPUTS
  # vcf_geno_vec = vector of vcf genotypes with following format:
  #   Genotype:Seq_Depth:Ref_depth,Alt_depth. Ex: 0/1:40:20,10
  # n_reads = the number of reads to subsample
  # seq_err = the estimated sequencing error rate used for calculating
  #            conditional probabilities of the genotypes
  # OUTPUT
  # vector of dosage genotypes: 0 = RR; 1 = RA; 2 = AA
  ################
  test_counts_list <- get_vcf_allele_counts(vcf_geno_vec)
  test_subsamp_counts <- generate_subsamp_counts(test_counts_list, 
    n_reads = n_reads)
  test_dose_vec <- generate_genotype_dosages(test_subsamp_counts, 
    seq_err = seq_err)
  return(test_dose_vec)
}
#test_subsamp_dose <- generate_subsamp_dosages(t14_Illum_vcf[,10], n_reads = 30)
###################
subsamp_PB_counts <- function(PB_count_vec, n_reads){
  # Function to subsample PacBio reads
  # INPUTS
  # PB_count_vec = vector of the allele counts for one sample from a "count"
  #                 files that I generated from the mpileup results of the
  #                 PacBio output; in format 'R,A1,A2,A3'
  # n_reads = the number of reads to subsample
  # OUTPUT
  # List of the subsample allele counts in same format at input
  #######################
  count_list <- sapply(PB_count_vec, function(x) 
    as.numeric(unlist(strsplit(x, split = ',')))
  )
  full_allele_vec_list <- list()
  for(i in seq(length(count_list))){
    tmp_vec <- paste(rep('R', times = count_list[[i]][1]))
    for(j in c(2:length(count_list[[i]]))){
      tmp_allele <- paste('A', (j-1), sep = '')
      tmp_vec <- c(tmp_vec, rep(tmp_allele, times = count_list[[i]][j])) 
    }
    full_allele_vec_list[[i]] <- tmp_vec
  }
  names(full_allele_vec_list) <- names(count_list)
  subsamp_allele_list <- lapply(full_allele_vec_list, function(x)
    if(length(x) < n_reads){return(NA)} else{sample(x, size = n_reads)}
  )
  subsamp_count_list <- lapply(subsamp_allele_list, function(x)
    c(sum(x == 'R'), sum(x == 'A1'), sum(x == 'A2'), sum(x == 'A3')))
  return(subsamp_count_list)
}
# test_pb_subsamp <- subsamp_PB_counts(t14_test_counts_0[,1], 25)
####################
####################



