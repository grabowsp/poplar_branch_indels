# Script used for analysing the pileups made from the PacBio data and 
#   Sujan's SNPs sets

# Load allele count file
t14_d42_s2_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth42_allele_counts.txt'
t14_d42_s2_c_0 <- read.table(t14_d42_s2_in, header = F, sep = '\t', 
  stringsAsFactors = F)

t14_d42_s2_vcf_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/tree14_depth42_v2.vcf'

t14_d42_s2_vcf <- read.table(t14_d42_s2_vcf_in, stringsAsFactors= F)

# Load the sample order, lib-to-branch names, and set branch-name order
samp_ord_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pb_pileup_results/sample_order.txt'

samp_ord <- unlist(
  read.table(samp_ord_in, header = F, stringsAsFactors = F))

lib_map_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pop_branch_pb_lib_names.txt'
lib_map <- read.table(lib_map_in, header = T, sep = '\t', stringsAsFactors = F)

lib_map$b_name <- gsub('.', '_', as.character(lib_map$branch_name), fixed = T)

l_ord <- c()
for(i in seq(length(samp_ord))){
  tmp_ind <- which(lib_map$lib_name == samp_ord[i])
  l_ord <- c(l_ord, tmp_ind)
}

b_ord <- lib_map$b_name[l_ord]

####
make_count_mat <- function(counts_df, nread){
  # Function to count up the number of alleles in a SNP above a count
  #  threshold
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

make_per_mat <- function(counts_df, per_cut){
  # Function to count up the number of alleles in a SNP above a percentage
  #  of all reads at the SNP
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

make_tot_depth_mat <- function(counts_df){
  # Calculate total depth for SNP
  #########3
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

######

t14_d42_s2_mat1 <- make_per_mat(counts_df = t14_d42_s2_c_0, per_cut = 0.1)
colnames(t14_d42_s2_mat1) <- b_ord

t14_d42_s2_mat2 <- make_per_mat(counts_df = t14_d42_s2_c_0, per_cut = 0.2)
colnames(t14_d42_s2_mat1) <- b_ord

t14_d42_s2_depth <- make_tot_depth_mat(counts_df = t14_d42_s2_c_0)

t14_d42_s2_meandepth <- apply(t14_d42_s2_depth, 1, mean)

t14_d42_s2_high_inds <- which(t14_d42_s2_meandepth > 100)

# t14_d42_s2_mat2[high_inds,]
# t14_d42_s2_mat2[-high_inds,]

t14_d42_s2_n_genos <- apply(t14_d42_s2_mat2, 1, function(x) 
  length(unique(x)))

sum(t14_d42_s2_n_genos == 1)
# 125

sum(t14_d42_s2_mat2[t14_d42_s2_n_genos == 1,1] == 1)
# 107 are homozygous
# 18 are heterozygous

sum(t14_d42_s2_n_genos[t14_d42_s2_high_inds] == 1)
# 32 (of 38) SNPs invariant

sum(t14_d42_s2_n_genos[-t14_d42_s2_high_inds] == 1)
# 93 (of 146) SNPs invariant

t14_d42_s2_down_dist <- abs((t14_d42_s2_vcf[c(1:(nrow(t14_d42_s2_vcf)-1)),2] - 
  t14_d42_s2_vcf[c(2:nrow(t14_d42_s2_vcf)),2]))  

t14_d42_s2_dist_mat <- matrix(NA, nrow = nrow(t14_d42_s2_vcf), ncol = 2)
t14_d42_s2_dist_mat[c(1:(nrow(t14_d42_s2_dist_mat)-1)), 1] <- t14_d42_s2_down_dist

t14_d42_s2_up_dist <- abs((t14_d42_s2_vcf[c(2:nrow(t14_d42_s2_vcf)),2] - 
  t14_d42_s2_vcf[c(1:(nrow(t14_d42_s2_vcf)-1)),2]))

t14_d42_s2_dist_mat[c(2:nrow(t14_d42_s2_dist_mat)), 2] <- t14_d42_s2_up_dist

t14_d42_s2_dist_1000_inds <- sort(
  unique(union(which(t14_d42_s2_dist_mat[,1] < 1000), 
  which(t14_d42_s2_dist_mat[,2] < 1000))))

length(t14_d42_s2_dist_1000_inds)
# 80

t14_d42_s2_dist_100_inds <- sort(
  unique(union(which(t14_d42_s2_dist_mat[,1] < 100),
  which(t14_d42_s2_dist_mat[,2] < 100)))

length(t14_d42_s2_dist_100_inds)
# 73

# next steps:
# number of SNPs that are invariant
# which SNPs have "good" patterns in tree 14
# which SNPs have "bad" patterns in tree 14
# which SNPs have "good" patterns in tree13
# which SNPs have "bad" patterns in tree13



