# Script for checking the distributions of somatic SNPs and genotypes in
#  the branches of Tree13 and Tree14

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

test_df <- t13_full

# find genotypes assigned an N in sujan's pattern
nodata_inds <- grep('N', test_df[,9])

t_d1 <- test_df[-nodata_inds,]

# check for genotypes that did not come up as a N for some reason
rogue_nas <- grep('-:-:-', apply(t_d1[,5:8], 1, function(x)
  paste(x, sep = '_', collapse = '_')))

t_d2 <- t_d1[-rogue_nas, ]
nrow(t_d2)
# 2833785

table(t_d2$V9)
#   HHHH   HHHR    HHHV    HHRH    HHRR    HHVH    HHVV    HRHH    HRHR    HRRH 
#2747889  31634    1034   12048    5617     213      87     435     489     329 
#   HRRR   HVHH    HVHV    HVVH    HVVV    RHHH    RHHR    RHHV    RHRH    RHRR 
#    605      3       8       5      32   15488    7444       1    4389    2412 
#   RRHH   RRHR    RRRH    RRRR    VHHH    VHHV    VHVH    VHVV    VVHH    VVHV 
#    341    725     531     280     288      96      47     333       6      32 
#   VVVH    VVVV 
#     27     917 

# HRRR: 605
# RHRR: 2412
# RRHR: 725
# RRRH: 531

######################

t13_info_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_v2SNPs_combined_info.rds',
  sep = '')
t13_info <- readRDS(t13_info_file)
nrow(t13_info)
# 21517

table(t13_info$tree13_Illum_genos_1)
#HHHR HHHV HHRH HHRR HHVH HHVV HRHH HRHR HRRH HRRR HVHV HVVH HVVV RHHH RHHR RHRH 
# 4737   25 2894 2039   24   26  342  404  282  531   4    2   17 3538 2550 1755 
# RHRR RRHH RRHR RRRH VHHH VHHV VHVH VHVV VVHH VVHV VVVH 
#  811  299  627  466   15   31    8   51    4   15   20 

# HRRR: 531
# RHRR: 811
# RRHR: 627
# RRRH: 466

table(t13_info$t13_dose_pattern_1)
# Too many to paste
# HRRR: 656
# RHRR: 299
# RRHR: 722
# RRRH: 753

pb_full_inds <- which(t13_info$t13_PB_NAs == 0)
length(pb_full_inds)
# 12504

table(t13_info$tree13_Illum_genos_1[pb_full_inds])
# HRRR: 221
# RHRR: 380
# RRHR: 247
# RRRH: 193

table(t13_info$t13_dose_pattern_1[pb_full_inds])
# HRRR: 380
# RHRR: 190
# RRHR: 394
# RRRH: 437

table(t13_info$tree13_PB_genos[pb_full_inds])
# HRRR (2:1:1:1): 418
# RHRR (1:2:1:1): 244
# RRHR (1:1:2:1): 521
# RRRH (1:1:1:2): 255

dose_na_inds <- grep('N', t13_info$t13_dose_pattern_1)

dose_pb_full_inds <- setdiff(pb_full_inds, dose_na_inds)
length(dose_pb_full_inds)
# 10,151

table(t13_info$tree13_PB_genos[dose_pb_full_inds])
# HRRR (2:1:1:1): 342
# RHRR (1:2:1:1): 203
# RRHR (1:1:2:1): 415
# RRRH (1:1:1:2): 210

t14_2_inds <- grep('2', t13_info$tree14_PB_genos)

full_t13_not14_inds <- setdiff(dose_pb_full_inds, t14_2_inds)
length(full_t13_not14_inds)
# 5023

table(t13_info$t13_dose_pattern_1[full_t13_not14_inds])
# HRRR: 252
# RHRR: 123
# RRHR: 264
# RRRH: 303

table(t13_info$tree13_PB_genos[full_t13_not14_inds])
# HRRR (2:1:1:1): 192
# RHRR (1:2:1:1): 87
# RRHR (1:1:2:1): 245
# RRRH (1:1:1:2): 94

match_inds_1 <- match_full_inds <- intersect(full_t13_not14_inds,
  which(t13_info$tree13_PB_genos == t13_info$t13_dose_pattern_2))
length(match_inds_1)
# 113

table(t13_info$t13_dose_pattern_1[match_inds_1])
# HRRR (2:1:1:1): 21
# RHRR (1:2:1:1): 20
# RRHR (1:1:2:1): 28
# RRRH (1:1:1:2): 22

t13_info[intersect(match_inds_1, which(t13_info$t13_dose_pattern_1 == 'RRHR')),
  c('V1', 'V2', 't13_dose_pattern_1')]
