# Script for generating tables that contain the distance to the next closest
#  indel for all pairwise comparisons for a sample

# INPUT FILES #
r_funct_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_tools/pop_branch_functions.r'

# SET VARIABLES #
file_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/shared_loci/'
file_prefix <- 'popbranch_'
file_suffix <- '.diff_indels'
split_chr <- '_v_'

# SET CONSTANTS #

# LOAD LIBRARIES AND PACKAGES #
source(r_funct_file)

############
# get names of all files to be processed
file_names <- system(paste('ls ', file_dir, '*', file_suffix, sep = ''), 
  intern = T)

# get the names of all the samples
samp_comps <- gsub(file_suffix, '',
  gsub(paste(file_dir, file_prefix, sep = ''), '', file_names), fixed = T)
all_samps <- sort(unique(unlist(strsplit(samp_comps, split = '_v_'))))

# make and save a table for each of the samples
for(smpnm in all_samps){
  print(paste('starting', smpnm))
  tot_tab <- gen_samp_tot_diff_tab(test_samp = smpnm, file_dir = file_dir, 
    file_prefix = file_prefix, file_suffix = file_suffix)
  out_file <- paste(file_dir, file_prefix, smpnm, '.indel_dist_tot', sep = '')
  write.table(tot_tab, file = out_file, quote = F, sep = '\t', 
    row.names = F, col.names = T)
  print(paste('finishing', smpnm))
}

quit(save = 'no')
