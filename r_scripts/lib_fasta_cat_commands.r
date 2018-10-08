# Script to make the 'cat' commands to merge the fasta files for each
#  sample

# LOAD PACKAGES #

# LOAD DATA #
seq_meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/pop_branch_sequencer_info.tsv'
seq_meta <- read.table(seq_meta_file, header = T, stringsAsFactors = F, 
  sep = '\t')

samp_meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v3.0.txt'
samp_meta <- read.table(samp_meta_file, header = T, stringsAsFactors = F, 
  sep = '\t')

# SET VARIABLES #
merge_fasta_dir <- '/home/f1p1/tmp/poplar_branches/lib_mapping/'

# SET CONSTANTS #

# SET OUTPUT #
out_file <- paste(merge_fasta_dir, '/cat_popbranch_subread_fastas.sh', sep = '')

###########
cat_vec <- c()
for(libname in samp_meta$lib_name){
  tmp_seq_inds <- grep(libname, seq_meta$Description)
  tmp_fasta_dirs <- seq_meta$Path[tmp_seq_inds]
  tmp_fasta_files <- paste(tmp_fasta_dirs, '/*subreads.fasta.gz', sep = '')
  tmp_out_file <- paste(merge_fasta_dir, libname, '_merged_fasta.gz', sep = '')
  tmp_cat_com <- paste('cat', paste(tmp_fasta_files, collapse = ' '), '>', 
    tmp_out_file, sep = ' ')
  cat_vec <- c(cat_vec, tmp_cat_com)
}

top_line <- '#!/bin/bash'
tot_com <- c(top_line, cat_vec)

write.table(tot_com, file = out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)

quit(save = 'no')

