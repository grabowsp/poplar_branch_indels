# Script for checking the distributions of somatic SNPs and genotypes in
#  the branches of Tree13 and Tree14

# Function file
func_file <- paste('/home/grabowsky/tools/workflows/poplar_branch_indels/',
  'checking_SNPs/', 'pb_SNPanalysis_funcs.r', sep = '')
source(func_file)

##########
# Tree14
t14_full_file <- t14_file <- paste('/home/grabowsky_scratch/',
  'poplar_branch_files/snps_v2/sujan_092519/', 'tree14_pvalues_patterns.txt',
  sep = '')

t14_full <- read.table(t14_full_file, header = F, stringsAsFactors = F)

nrow(t14_full)
# 3258363

test_df <- t14_full

# find genotypes assigned an N in sujan's pattern
nodata_inds <- grep('N', test_df[,9])

t_d1 <- test_df[-nodata_inds,]

# check for genotypes that did not come up as a N for some reason
rogue_nas <- grep('-:-:-', apply(t_d1[,5:8], 1, function(x)
  paste(x, sep = '_', collapse = '_')))

t_d2 <- t_d1[-rogue_nas, ]
nrow(t_d2)
# 2972701 

table(t_d2$V9)

#   HHHH   HHHR    HHHV    HHRH    HHRR    HHVH    HHVV    HRHH    HRHR    HRRH 
#2887944   8430     229    7695    2557     134     250   36143    4761    4458 
#   HRRR   HVHH    HVHV    HVVH    HVVV    RHHH    RHHR    RHRH    RHRR    RHRV 
#   2142    656      80      67     298    4278    1391    1352    1353       1 
#   RHVV   RRHH    RRHR    RRRH    RRRR    RVRR    VHHH    VHHV    VHVH    VHVV 
#      1   3115    1616    1430     674       1      36      15      11      33 
#   VVHH    VVHV    VVVH    VVVV 
#     44     101      74    1331 

# HRRR: 2142 
# RHRR: 1353
# RRHR: 1616
# RRRH: 1430

######################

t14_info_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree14_v2SNPs_combined_info.rds',
  sep = '')
t14_info <- readRDS(t14_info_file)
nrow(t14_info)
# 21517

table(t14_info$tree14_Illum_genos_1)
# HRRR: 1147
# RHRR: 1133
# RRHR: 999
# RRRH: 913

table(t14_info$t14_dose_pattern_1)
# HRRR: 485
# RHRR: 1377
# RRHR: 498
# RRRH: 530

pb_full_inds <- which(t14_info$t14_PB_NAs == 0)
length(pb_full_inds)
# 16982

table(t14_info$tree14_Illum_genos_1[pb_full_inds])
# HRRR: 709
# RHRR: 646
# RRHR: 574
# RRRH: 509

table(t14_info$t14_dose_pattern_1[pb_full_inds])
# HRRR: 329
# RHRR: 983
# RRHR: 320
# RRRH: 343

table(t14_info$tree14_PB_genos[pb_full_inds])
# HRRR (2:1:1:1): 459
# RHRR (1:2:1:1): 512
# RRHR (1:1:2:1): 443
# RRRH (1:1:1:2): 553

dose_na_inds <- grep('N', t14_info$t14_dose_pattern_1)

dose_pb_full_inds <- setdiff(pb_full_inds, dose_na_inds)
length(dose_pb_full_inds)
# 13,684

table(t14_info$tree14_PB_genos[dose_pb_full_inds])
# HRRR (2:1:1:1): 359
# RHRR (1:2:1:1): 409
# RRHR (1:1:2:1): 348
# RRRH (1:1:1:2): 466

t13_2_inds <- grep('2', t14_info$tree13_PB_genos)

full_t14_not13_inds <- setdiff(dose_pb_full_inds, t13_2_inds)
length(full_t14_not13_inds)
# 6604

table(t14_info$t14_dose_pattern_1[full_t14_not13_inds])
# HRRR: 163
# RHRR: 743
# RRHR: 198
# RRRH: 227

table(t14_info$tree14_PB_genos[full_t14_not13_inds])
# HRRR (2:1:1:1): 195
# RHRR (1:2:1:1): 262
# RRHR (1:1:2:1): 192
# RRRH (1:1:1:2): 249

match_inds_1 <- match_full_inds <- intersect(full_t14_not13_inds,
  which(t14_info$tree14_PB_genos == t14_info$t14_dose_pattern_2))
length(match_inds_1)
# 124

table(t14_info$t14_dose_pattern_1[match_inds_1])
# HRRR (2:1:1:1): 13
# RHRR (1:2:1:1): 37
# RRHR (1:1:2:1): 22
# RRRH (1:1:1:2): 23

t14_info[intersect(match_inds_1, which(t14_info$t14_dose_pattern_1 == 'RRHR')),
  c('V1', 'V2', 't14_dose_pattern_1')]



