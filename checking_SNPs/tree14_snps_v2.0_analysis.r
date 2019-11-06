# Analysis of files that Sujan generated using the new criteria

# Function file
func_file <- paste('/home/grabowsky/tools/workflows/poplar_branch_indels/', 
  'checking_SNPs/', 'pb_SNPanalysis_funcs.r', sep = '')
source(func_file)

##########
# Tree14
t14_full_file <- t14_file <- '/home/grabowsky_scratch/poplar_branch_files/snps_v2/sujan_092519/tree14_pvalues_patterns.txt'

t14_full <- read.table(t14_full_file, header = F, stringsAsFactors = F)

nrow(t14_full)
# 3,258,363

# find genotypes assigned an N in sujan's pattern
nodata_inds <- grep('N', t14_full[,9])
length(nodata_inds)
# 240821
t14_d1 <- t14_full[-nodata_inds,]

# check for genotypes that did not come up as a N for some reason
rogue_nas <- grep('-:-:-', apply(t14_d1[,5:8], 1, function(x)
  paste(x, sep = '_', collapse = '_')))
length(rogue_nas)
# 44,841

# total number of missing genotypes
sum(length(nodata_inds), length(rogue_nas))
# 285662

sum(length(nodata_inds), length(rogue_nas))/nrow(t14_full)
# [1] 0.0876704
# 8.8% of the SNPs have missing data

# Missing genotypes by sample
apply(t14_full[, c(5:8)], 2, function(x) length(grep('-:-:-', x)))
#    V5     V6     V7     V8 
# 55777  36854 274073  19508

t14_d2 <- t14_d1[-rogue_nas, ]
nrow(t14_d2)
# 2,972,701

no_var_inds <- grep('HHHH|VVVV|RRRR', t14_d2[,9])
length(no_var_inds)
# 2,889,949
length(no_var_inds)/nrow(t14_d2)
# 97.2%
sum(t14_d2[,9] == 'RRRR')
# 674
sum(t14_d2[,9] == 'VVVV')
# 1331
sum(t14_d2[,9] == 'HHHH')
# 2,887,944

t14_d3 = t14_d2[-no_var_inds, ]
nrow(t14_d3)
# 82,752

# use the cutoff for a "good" genotype as 10,000X better than others
test_cut <- 1- (1/10000)
t14_good_snps <- get_good_snps(geno_df = t14_d3, good_cut = test_cut)

length(t14_good_snps)
# 24,857

############
# assign genotypes for "good" SNPs based on probabilities and dosages

t14_dosage_mat <- gen_geno_dosages(t14_d3[t14_good_snps, ])
colnames(t14_dosage_mat) <- c('b14.2','b14.3','b14.4','b14.5')
t14_dosage_df <- data.frame(chr = t14_d3[t14_good_snps, 1], 
  pos = t14_d3[t14_good_snps, 2], t14_dosage_mat, stringsAsFactors = F)

############

table(t14_d3[t14_good_snps, 9])
#HHHR HHHV HHRH HHRR HHVH HHVV HRHH HRHR HRRH HRRR HVHH HVHV HVVH HVVV RHHH RHHR 
#2620   19 2538 1478   14  167 4899 1690 1741 1154   14   24   17  169 1941  857 
#RHRH RHRR RHVV RRHH RRHR RRRH VHHH VHHV VHVH VHVV VVHH VVHV VVVH 
# 822 1136    1 1525 1001  921    2    4    4   23    9   36   31

miss_data_vec <- c(55777,36854,274073,19508)
H3R1_vec <- c(2620, 2538, 4899, 1941)
cor.test(miss_data_vec, H3R1_vec)
# r = 0.9912; (0.637, 0.9998); Not linear, but certainly there is a correlation

sum(H3R1_vec)
# 11,998
sum(H3R1_vec)/length(t14_good_snps)
# 48.3%

H1R3_vec <- c(921, 1001, 1136, 1154)
sum(H1R3_vec)
# 4212
sum(H1R3_vec)/length(t14_good_snps)
# 16.9%

# Generate list of SNP positions to include
t14_out_pos_file <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree14_good_positions_v1.txt', sep = '')

write.table(t14_d3[t14_good_snps, c(1:2)], file = t14_out_pos_file, quote = F,
  sep = '\t', row.names = F, col.names = F)

t14_good_df <- t14_d3[t14_good_snps, ]

for(cn in unique(t14_good_df[,1])){
  tmp_out_file <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree14_good_pos_', cn, '_v1.txt', sep = '')
  tmp_inds <- which(t14_good_df[,1] == cn)
  write.table(t14_good_df[tmp_inds, c(1:2)], file = tmp_out_file, quote = F,
    sep = '\t', row.names = F, col.names = F)
}

# Try to generate pseudo-VCF for Jerry's approach
## Will generate the genotype info here and have to add the header lines later
psvcf <- data.frame(CHROM = t14_good_df[, 1], POS = t14_good_df[, 2], 
  ID = paste(t14_good_df[, 1], t14_good_df[, 2], sep = '_'), 
  REF = t14_good_df[, 3], ALT = t14_good_df[, 4],
  stringsAsFactors = F
  )

psvcf$QUAL <- '.'
psvcf$FILTER <- '.'
psvcf$INFO <- '.'
psvcf$FORMAT <- 'GT:DP:AD'
psvcf$tree14_2 <- '0/1:10:5,5'

pseudo_vcf_out <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree14_good_pos_vcf_meat.txt', sep = '')

write.table(psvcf, file = pseudo_vcf_out, quote = F, sep = '\t', 
  row.names = F, col.names = F)

psvcf_chr1 <- psvcf[which(psvcf[,1] == 'Chr01'), ]

pseudo_vcf_chr1_out <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree14_good_pos_vcf_Chr01_meat.txt', sep = '')

write.table(psvcf_chr1, file = pseudo_vcf_chr1_out, quote = F, sep = '\t', 
  row.names = F, col.names = F)

###
psvcf <- read.table(pseudo_vcf_out, header = F, sep = '\t', 
  stringsAsFactors = F)
####

# Write genotype dosages for good SNPs

t14_good_dosages_file <- paste('/home/grabowsky_scratch/poplar_branch_files/', 
  'snps_v2/sujan_092519/', 'tree14_good_positions_v1.dosages', sep = '')

write.table(t14_dosage_df, file = t14_good_dosages_file, quote = F, sep = '\t',
  row.names = F, col.names = T)


