# Script for analyzing and processing results from pairwise comparisons
#  that are outputted from the --diff-site option in vcftools of

# INPUT FILES #
r_function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_tools/pop_branch_functions.r'

# SET VARIABLES #
file_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/shared_loci/'
in_file_prefix <- 'popbranch_'

# SET CONSTANTS #
in_file_suffix <- '.diff.sites_in_files'

# LOAD LIBRARIES AND PACKAGES #
source(r_function_file)

#################
syst_command <- paste('ls ', file_dir, in_file_prefix, '*', in_file_suffix, sep = '')

diff_files <- system(syst_command, intern = T)

samp_comps <- unique(gsub('_DEL|_INS', '', gsub(in_file_suffix, '', 
  gsub(paste(file_dir, in_file_prefix, sep = ''), '', diff_files), fixed = T)))

for(sc in samp_comps){
  test_del_file <- diff_files[grep(paste(sc, '_DEL', sep = ''), diff_files)]
  test_ins_file <- diff_files[grep(paste(sc, '_INS', sep = ''), diff_files)]
  # proc_share_df() is function to get indel size and insert indel type column
  ins_df <- proc_share_df(share_df = read.table(test_ins_file, header = T, 
              stringsAsFactors = F), indel_type = 'INS')
  del_df <- proc_share_df(share_df = read.table(test_del_file, header = T, 
              stringsAsFactors = F), indel_type = 'DEL')
  combo_df <- rbind(ins_df, del_df)
  # dist_to_other_vars() is, for variants found only in 1 samples, to find
  #  the next closest variant in same/diff sample and of same/diff indel type
  test_dist_df <- dist_to_other_vars(combo_df)
  # get the shortest distance to next closest variant in either sample
  test_dist_df$min_var_dist <- apply(
    test_dist_df[, c('close_altSamp_sameVar', 'close_sameSamp_sameVar', 
    'close_altSamp_altVar', 'close_sameSamp_altVar')], 1, min)
  out_file <- paste(file_dir, in_file_prefix, sc, '.diff_indels', sep = '')
  write.table(test_dist_df, file = out_file, quote = F, sep = '\t', 
    row.names = F, col.names = T)
}

quit(save = 'no')
