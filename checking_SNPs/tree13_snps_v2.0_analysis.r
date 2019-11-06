# Analysis of files that Sujan generated using the new criteria

# Function file
func_file <- paste('/home/grabowsky/tools/workflows/poplar_branch_indels/', 
  'checking_SNPs/', 'pb_SNPanalysis_funcs.r', sep = '')
source(func_file)

##########
# Tree14
t13_full_file <- t14_file <- paste('/home/grabowsky_scratch/', 
  'poplar_branch_files/snps_v2/sujan_092519/', 'tree13_pvalues_patterns.txt', 
  sep = '')

t13_full <- read.table(t13_full_file, header = F, stringsAsFactors = F)

nrow(t13_full)
# 3,245,234

test_df <- t13_full

# find genotypes assigned an N in sujan's pattern
nodata_inds <- grep('N', test_df[,9])
length(nodata_inds)
# 356,831
t_d1 <- test_df[-nodata_inds,]

# check for genotypes that did not come up as a N for some reason
rogue_nas <- grep('-:-:-', apply(t_d1[,5:8], 1, function(x)
  paste(x, sep = '_', collapse = '_')))
length(rogue_nas)
# 54,618

# total number of missing genotypes
sum(length(nodata_inds), length(rogue_nas))
# 411,449

sum(length(nodata_inds), length(rogue_nas))/nrow(test_df)
# [1] 0.1267856
# 12.7% of the SNPs have missing data

# Missing genotypes by sample
apply(test_df[, c(5:8)], 2, function(x) length(grep('-:-:-', x)))
#     V5     V6     V7     V8 
# 375197  87947  10814 104897 

t_d2 <- t_d1[-rogue_nas, ]
nrow(t_d2)
# 2,833,785

no_var_inds <- grep('HHHH|VVVV|RRRR', t_d2[,9])
length(no_var_inds)
# 2,749,086
length(no_var_inds)/nrow(t_d2)
# 97.0 %
sum(t_d2[,9] == 'RRRR')
# 280
sum(t_d2[,9] == 'VVVV')
# 917
sum(t_d2[,9] == 'HHHH')
# 2,747,889

t_d3 = t_d2[-no_var_inds, ]
nrow(t_d3)
# 84,699

# use the cutoff for a "good" genotype as 10,000X better than others
test_cut <- 1- (1/10000)
t_good_snps <- get_good_snps(geno_df = t_d3, good_cut = test_cut)

length(t_good_snps)
# 21,583

############
# assign genotypes for "good" SNPs
t13_dosage_mat <- gen_geno_dosages(t_d3[t_good_snps, ])
colnames(t13_dosage_mat) <- c('b13.1','b13.2','b13.3','b13.5')
t13_dosage_df <- data.frame(chr = t_d3[t_good_snps, 1],
  pos = t_d3[t_good_snps, 2], t13_dosage_mat, stringsAsFactors = F)

############

table(t_d3[t_good_snps, 9])
# Can do this later

# need to adjust these values for tree13 values
#miss_data_vec <- c(55777,36854,274073,19508)
#H3R1_vec <- c(2620, 2538, 4899, 1941)
#cor.test(miss_data_vec, H3R1_vec)
# r = 0.9912; (0.637, 0.9998); Not linear, but certainly there is a correlation

#sum(H3R1_vec)
# 11,998
#sum(H3R1_vec)/length(t14_good_snps)
# 48.3%

#H1R3_vec <- c(921, 1001, 1136, 1154)
#sum(H1R3_vec)
# 4212
#sum(H1R3_vec)/length(t14_good_snps)
# 16.9%

# Generate list of SNP positions to include
t13_out_pos_file <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree13_good_positions_v1.txt', sep = '')

write.table(t_d3[t_good_snps, c(1:2)], file = t13_out_pos_file, quote = F,
  sep = '\t', row.names = F, col.names = F)

t_good_df <- t_d3[t_good_snps, ]

# Try to generate pseudo-VCF for Jerry's approach
## Will generate the genotype info here and have to add the header lines later
psvcf <- data.frame(CHROM = t_good_df[, 1], POS = t_good_df[, 2], 
  ID = paste(t_good_df[, 1], t_good_df[, 2], sep = '_'), 
  REF = t_good_df[, 3], ALT = t_good_df[, 4],
  stringsAsFactors = F
  )

psvcf$QUAL <- '.'
psvcf$FILTER <- '.'
psvcf$INFO <- '.'
psvcf$FORMAT <- 'GT:DP:AD'
psvcf$tree13_2 <- '0/1:10:5,5'

pseudo_vcf_out <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree13_good_pos_vcf_meat.txt', sep = '')

write.table(psvcf, file = pseudo_vcf_out, quote = F, sep = '\t', 
  row.names = F, col.names = F)

########

# Write genotype dosages for good SNPs

t13_good_dosages_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_good_positions_v1.dosages', sep = '')

write.table(t13_dosage_df, file = t13_good_dosages_file, quote = F, sep = '\t',
  row.names = F, col.names = T)

