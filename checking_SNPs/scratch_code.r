func_file <- paste('/home/grabowsky/tools/workflows/poplar_branch_indels/',
  'checking_SNPs/', 'pb_SNPanalysis_funcs.r', sep = '')
source(func_file)

t14_Illum_vcf_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree14.good_positions.v1.vcf', sep = '')
t14_Illum_vcf <- read.table(t14_Illum_vcf_file, stringsAsFactors = F)

unlist(strsplit(t14_Illum_vcf[1,10], split = ':'))[3]

s1_counts <- sapply(t14_Illum_vcf[,10], function(x)
  unlist(strsplit(x, split = ':'))[3])

s1_counts_list <- strsplit(s1_counts, split = ',')

tmp_readcount <- c(rep('R', times = as.numeric(s1_counts_list[[1]][1])),
  rep('A', times = as.numeric(s1_counts_list[[1]][2])))

full_readcount_list <- lapply(s1_counts_list, function(x)
  c(rep('R', times = as.numeric(x)[1]), rep('A', times = as.numeric(x[2])))
)

subsamp_readcount_list <- lapply(full_readcount_list, function(x)
  if(length(x) < 30){return(NA)}else{sample(x, size = 30)})

subsamp_count_list <- lapply(subsamp_readcount_list, function(x) 
  c(sum(x == 'R'), sum(x == 'A')))

# Use ALT to calc prob of R:R
hom_ref_prob_list <- lapply(subsamp_count_list, function(x) 
  dbinom(x[2], sum(x), prob = 0.01))
het_prob_list <- lapply(subsamp_count_list, function(x)
  dbinom(x[2], sum(x), prob = 0.5))
hom_alt_prob_list <- lapply(subsamp_count_list, function(x) 
  dbinom(x[1], sum(x), prob = 0.01))

prob_df <- data.frame(p_RR = unlist(hom_ref_prob_list), 
  p_RA = unlist(het_prob_list), p_AA = unlist(hom_alt_prob_list), 
  stringsAsFactors = F)

dose_vec <- apply(prob_df, 1, function(x) (1*x[2] + 2*x[3])/sum(x))


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



