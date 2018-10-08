# Script to generate commands for vcftools to separate the INSertion and 
#   DELetions into separate files
#  and type info

# LOAD FILES #
meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v1.0.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t')


# SET VARIABLES #
indel_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info/'

# SET CONSTANTS #

# LOAD LIBARIES AND PACKAGES #


########################
branch_names <- gsub('.', '_', meta$branch_name, fixed = T)

for(i in branch_names){
  info_file <- paste(indel_dir, 'branch_', i, '_chrom.INFO', sep = '')
  info <- read.table(info_file, header = T, stringsAsFactors = F, sep = '\t')
  chr_inds <- grep('Chr', info$CHROM)
  del_inds <- which(info$SVTYPE == 'DEL')
  del_chr_inds <- intersect(del_inds, chr_inds)
  ins_inds <- which(info$SVTYPE == 'INS')
  ins_chr_inds <- intersect(ins_inds, chr_inds)
  del_pos_file <- paste(indel_dir, 'branch_', i, '_del_positions.txt', sep = '')
  write.table(info[del_chr_inds, c(1:2)], file = del_pos_file, quote = F, 
    sep = '\t', row.names = F, col.names = F)
  ins_pos_file <- paste(indel_dir, 'branch_', i, '_ins_positions.txt', sep = '')
  write.table(info[ins_chr_inds, c(1:2)], file = ins_pos_file, quote = F, 
    sep = '\t', row.names = F, col.names = F)
}

quit(save = 'no')

