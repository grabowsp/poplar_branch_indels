# Script for using subsampling to try to get a less biased somatic SNP result

## LOAD FUNCTIONS AND PACKAGES
func_file <- paste('/home/grabowsky/tools/workflows/poplar_branch_indels/',
  'checking_SNPs/', 'pb_SNPanalysis_funcs.r', sep = '')
source(func_file)

###############

## LOAD DATA

# PacBio allele counts
t13_PB_counts_in <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_v2SNPs_PB_full_allele_counts.txt',
  sep = '')
t13_PB_counts_0 <- read.table(t13_PB_counts_in, header = F, sep = '\t',
  stringsAsFactors = F)

# PacBio VCF made frmo mpileup results
t13_PB_vcf_in <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_v2SNPs_PB_full.vcf',
  sep = '')
t13_PB_vcf_0 <- read.table(t13_PB_vcf_in, stringsAsFactors = F)

# load Illunina SNP info based on file from Sujan
t13_Illum_snp_info_file <- paste('/home/grabowsky_scratch/',
  'poplar_branch_files/snps_v2/sujan_092519/', 'tree13_pval_file_SNPinfo.txt',
  sep = '')
t13_Illum_snp_info <- read.table(t13_Illum_snp_info_file, header = F,
  stringsAsFactors = F)

# load Illumina-based VCF
t13_Illum_vcf_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13.good_positions.v1.vcf', sep = '')
t13_Illum_vcf <- read.table(t13_Illum_vcf_file, stringsAsFactors = F)

# load the sample order
samp_ord_in <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree14_PB_sample_order.txt',
  sep = '')
samp_ord <- unlist(
  read.table(samp_ord_in, header = F, stringsAsFactors = F))

# load info about mapping library to sample names
lib_map_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/check_snps_with_pb/pop_branch_pb_lib_names.txt'
lib_map <- read.table(lib_map_in, header = T, sep = '\t', stringsAsFactors = F)

#################################
## ANALYSIS SCRIPT

# Remove INDELS and Scaffold SNPs from mpileup VCF and add snp_name
indel_inds <- grep('INDEL', t13_PB_vcf_0[,8])
scaff_inds <- grep('scaff', t13_PB_vcf_0[,1])
remove_inds <- union(indel_inds, scaff_inds)
t13_vcf <- t13_PB_vcf_0[-remove_inds,]
t13_vcf$snp_name <- paste(t13_vcf[ ,1], t13_vcf[ ,2], sep = '_')

# Remove INDELS and Scaffold SNPs from the count matrix
t13_PB_counts <- t13_PB_counts_0[-remove_inds, ]

# Add the Illumina Ref and Alt alleles to the PB mpileup VCF
t13_Illum_snp_info$snp_name <- paste(t13_Illum_snp_info[,1],
  t13_Illum_snp_info[,2], sep = '_')

illum_inds <- which(t13_Illum_snp_info$snp_name %in% t13_vcf$snp_name)

t13_Illum_sub_1 <- t13_Illum_snp_info[illum_inds, ]
illum_ord_1 <- order(t13_Illum_sub_1[,2])
illum_ord_2 <- order(t13_Illum_sub_1[illum_ord_1,1])

t13_Illum_sub_2 <- t13_Illum_sub_1[illum_ord_1[illum_ord_2], ]

# Check that the names match up, and if so, add Illumina Ref/Alt data to VCF
if( sum(t13_Illum_sub_2$snp_name == t13_vcf$snp_name) ==
 nrow(t13_vcf) ){
  t13_vcf$Illum_Ref <- t13_Illum_sub_2[,3]
  t13_vcf$Illum_Alt <- t13_Illum_sub_2[,4]
  } else {print('Illumina and PacBio snp_names do not match')}

# Calculate the PacBio Sequencing Depth
t13_depth_mat <- make_tot_depth_mat(counts_df = t13_PB_counts)
low_20_inds <- which(t13_depth_mat < 20)

# Assign PacBio Genotypes to each sample
t13_geno_mat_1 <- make_per_mat(counts_df = t13_PB_counts, per_cut = 0.1)

t13_geno_mat_20 <- t13_geno_mat_1
t13_geno_mat_20[low_20_inds] <- NA

# Generate PacBio Genotype Patterns
t13_char_geno_20 <- matrix(NA, nrow = nrow(t13_geno_mat_20), 
  ncol = ncol(t13_geno_mat_20))
t13_char_geno_20[which(is.na(t13_geno_mat_20))] <- 'N'
for(i in seq(4)){
  tmp_inds <- which(t13_geno_mat_20 == i)
  t13_char_geno_20[tmp_inds] <- as.character(i)
}

t14_geno_string <- apply(t13_char_geno_20[ , c(8,7,6,5)], 1, function(x)
  paste(x, collapse = ':'))
t13_geno_string <- apply(t13_char_geno_20[ , c(4,3,2,1)], 1, function(x)
  paste(x, collapse = ':'))

t13_vcf$tree14_PB_genos_20 <- t14_geno_string
t13_vcf$tree13_PB_genos_20 <- t13_geno_string

#########
# order of branches in Illumina VCF is 13.1, 13.2, 13.3, 13.5

##### CONTINUE ADAPTING FROM HERE #######

# Order of branches in Illumina VCF is: 14.5, 14.4, 14.3, 14.2
# Generate sub-sample dosage genotypes for different thresholds
sub_dose_list <- list()
for(i in c(20,25,30,35,40,45)){
  for(j in seq(3)){
    tmp_col_name <- paste('sub_illum_dose', i, j, sep = '_')
    tmp_illum_dose_mat <- apply(t13_Illum_vcf[ , c(13,12,11,10)], 2, 
      function(x) generate_subsamp_dosages(x, n_reads = i))
    tmp_dose_genos <- assign_geno_from_dosage(tmp_illum_dose_mat, 
      dose_dist_cut = 0.1)
    tmp_pattern_1 <- apply(tmp_dose_genos, 1, function(x) 
      paste(x, collapse = ''))
    tmp_pattern_2 <- gsub('R|V', '1', tmp_pattern_1)
    tmp_pattern_2 <- gsub('H', '2', tmp_pattern_2)
    tmp_pattern_2 <- unlist(lapply(
      strsplit(tmp_pattern_2, split = ''), function(x)
      paste(x, collapse = ':')
    ))
    sub_dose_list[[
      paste(tmp_col_name, '_pattern_1', sep = '')]] <- tmp_pattern_1
    sub_dose_list[[
      paste(tmp_col_name, '_pattern_2', sep = '')]] <- tmp_pattern_2
  }
  print(i)
}

# Find the Illumina SNPs that match the PB SNPs
t13_Illum_vcf$snp_name <- paste(t13_Illum_vcf[,1], t13_Illum_vcf[,2], 
  sep = '_')

illum_in_pb <- which(t13_Illum_vcf$snp_name %in% t13_vcf$snp_name)
pb_in_illum <- which(t13_vcf$snp_name %in% t13_Illum_vcf$snp_name)

# sum(t14_Illum_vcf$snp_name[illum_in_pb] == t14_vcf$snp_name[pb_in_illum])
# 24776 - the length of both vectors

# Generate data frame for overlap in PB and Illumina SNPs
t13_sub_df <- data.frame(Chr = t13_vcf[pb_in_illum, 1], 
  pos = t13_vcf[pb_in_illum,2], 
  snp_name = t13_vcf$snp_name[pb_in_illum],
  PB_ref = t13_vcf[pb_in_illum, 4], PB_alt = t13_vcf[pb_in_illum, 5],
  Illum_ref = t13_vcf$Illum_Ref[pb_in_illum],
  Illum_alt = t13_vcf$Illum_Alt[pb_in_illum],
  tree14_PB_genos_20 = t13_vcf$tree14_PB_genos_20[pb_in_illum],
  tree13_PB_genos_20 = t13_vcf$tree13_PB_genos_20[pb_in_illum],
  stringsAsFactors = F)

# add the subsampled doseage genotypes
for(i in seq(length(sub_dose_list))){
  tmp_name <- names(sub_dose_list)[i]
  t13_sub_df[, tmp_name] <- sub_dose_list[[i]][illum_in_pb]
}

# Look for matches in all the sublists

sing_het_string <- '2:1:1:1|1:2:1:1|1:1:2:1|1:1:1:2'
good_order_strings <- '1:2:2:2|1:1:2:2'

pb_full_inds <- setdiff(seq(nrow(t13_sub_df)), 
  grep('N', t13_sub_df$tree13_PB_genos_20))

# Need to make the below code into a loop for each subsamples genotype set
subsamp_match_list <- list()

for(i in c(20,25,30,35,40,45)){
  for(j in seq(3)){
    tmp_pre <- paste('sub_illum_dose_', i, '_', j, sep = '')
    tmp_col <- paste(tmp_pre, '_pattern_2', sep = '')
    tmp_col_1 <- paste(tmp_pre, '_pattern_1', sep = '')
    tmp_dose_full <- setdiff(seq(nrow(t13_sub_df)),
      grep('N', t13_sub_df[, tmp_col]))
    tmp_pb_dose_full <- intersect(pb_full_inds, tmp_dose_full)
    tmp_match_full <- intersect(tmp_pb_dose_full,
      which(t13_sub_df$tree13_PB_genos_20 == t13_sub_df[[tmp_col]]))
    tmp_13_2 <- grep('2', t13_sub_df$tree14_PB_genos_20)
    tmp_match_2 <- setdiff(tmp_match_full, tmp_13_2)
    tmp_one_het_inds <- grep(sing_het_string, t13_sub_df[, tmp_col])
    tmp_good_ord_inds <- grep(good_order_strings, t13_sub_df[, tmp_col])
    tmp_match_3 <- intersect(tmp_match_2, 
      c(tmp_one_het_inds, tmp_good_ord_inds))
#    subsamp_match_list[[tmp_pre]] <- list() 
    subsamp_match_list[[tmp_pre]][[1]] <- tmp_match_3
    subsamp_match_list[[tmp_pre]][[2]] <- table(
      t13_sub_df[tmp_match_3, tmp_col_1])
  }
}

for(i in seq(length(subsamp_match_list))){
  print(names(subsamp_match_list)[i])
  print(subsamp_match_list[[i]][[2]])
}

tree13_save_object <- list()
tree13_save_object[['tree13_subsampled_dataframe']] <- t13_sub_df
tree13_save_object[['tree13_subsampled_matches']] <- subsamp_match_list

t13_subsamp_data_out <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_subsample_info.rds',
  sep = '')
saveRDS(tree13_save_object, file = t13_subsamp_data_out)

######################
# save SNP positions
t13_subsamp_data_in <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_subsample_info.rds',
  sep = '')
tree13_info <- readRDS(t13_subsamp_data_in)

t13_vcf <- tree13_info[[1]]
t13_sub_info <- tree13_info[[2]]

sub_out_dir <- '/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/downsample/'

sub_pos_names <- paste('tree13.somatic.down_', 
  rep(c(20,25,30,35,40,45), each = 3),
  '.rep_', rep(c(1,2,3), times = 6), '.positions', sep = '')

for(i in seq(length(sub_pos_names))){
  tmp_out_name <- paste(sub_out_dir, sub_pos_names[i], sep = '')
  tmp_inds <- t13_sub_info[[i]][[1]]
  tmp_pos_mat <- t13_vcf[tmp_inds, c(1:2)]
  write.table(tmp_pos_mat, file = tmp_out_name, quote = F, sep = '\t', 
    row.names = F, col.names = F)
}

####################

# Figure out how many SNPs are included in each depth threshold

func_file <- paste('/home/grabowsky/tools/workflows/poplar_branch_indels/',
  'checking_SNPs/', 'pb_SNPanalysis_funcs.r', sep = '')
source(func_file)

t13_subsamp_data_in <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_subsample_info.rds',
  sep = '')
tree13_info <- readRDS(t13_subsamp_data_in)

t13_vcf <- tree13_info[[1]]
t13_sub_info <- tree13_info[[2]]

pb_n_inds <- grep('N', t13_vcf$tree13_PB_genos_20)
pb_full_inds <- setdiff(seq(nrow(t13_vcf)), pb_n_inds)

t14_2_inds <- grep('2', t13_vcf$tree14_PB_genos_20)


sub_colnames <- colnames(t13_vcf)[seq(from = 10, by = 2, length.out = 18)]

overlap_length_vec <- c()

for(i in seq(length(sub_colnames))){
  tmp_n_inds <- grep('N', t13_vcf[, sub_colnames[i]])
  tmp_full_1 <- setdiff(seq(nrow(t13_vcf)), tmp_n_inds)
  tmp_overlap <- setdiff(intersect(tmp_full_1, pb_full_inds), t14_2_inds)
  overlap_length_vec <- c(overlap_length_vec, length(tmp_overlap))
}

mean(overlap_length_vec[1:3])
# Depth 20: 3025

mean(overlap_length_vec[4:6])
# Depth 25: 3363

mean(overlap_length_vec[7:9])
# Depth 30: 2863.3333

mean(overlap_length_vec[10:12])
# Depth 35: 2437.667

mean(overlap_length_vec[13:15])
# Depth 40: 1871.667

mean(overlap_length_vec[16:18])
# Depth 45: 1368


