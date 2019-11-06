# Script to combine the "good" positions for Tree13 and Tree14

t13_pos_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree13_good_positions_v1.txt', sep = '')
t13_pos <- read.table(t13_pos_file, header = F, sep = '\t', 
  stringsAsFactors = F)

t14_pos_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'tree14_good_positions_v1.txt', sep = '')
t14_pos <- read.table(t14_pos_file, header = F, sep = '\t', 
  stringsAsFactors = F)

###########

t13_pos$snp_name <- paste(t13_pos[,1], t13_pos[,2], sep = '_')
t14_pos$snp_name <- paste(t14_pos[,1], t14_pos[,2], sep = '_')

combo_pos <- merge(t13_pos, t14_pos,  all = T)

nrow(combo_pos)
# 45,203

nrow(t13_pos) + nrow(t14_pos)
# 46440

length(intersect(t13_pos$snp_name, t14_pos$snp_name))
# 1237

# Check that all positions are accounted for
nrow(combo_pos) == ((nrow(t13_pos) + nrow(t14_pos)) - 
  length(intersect(t13_pos$snp_name, t14_pos$snp_name)))
# [1] TRUE

# Note: the positions seem to already be sorted

combo_out_pos_file <- paste('/home/grabowsky_scratch/poplar_branch_files/',
  'snps_v2/sujan_092519/', 'trees13_14_combo_good_positions_v1.txt', sep = '')

write.table(combo_pos[, c(1:2)], file = combo_out_pos_file, quote = F, 
  sep = '\t', row.names = F, col.names = F)
