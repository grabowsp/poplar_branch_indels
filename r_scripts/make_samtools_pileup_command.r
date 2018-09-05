# Script to help make the mpileup command for SAM Tools for the Poplar data

# LOAD FILES #
meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v3.0.txt'
meta_0 <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t')

# SET VARIABLES #
bad_samps <- c('13.4', '14.1')
# SET CONSTANTS #
base_file_dir <- '/home/smrtlink/userdata/jobs_root.local'
bam_file_string <- 'tasks/pbsvtools.tasks.gather_align-1/alignments.bam'

bam_sort_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams'

sam_out_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/multisamp_indel_calling'

ref_path <- '/home/smrtlink/references/Ptrichocarpa_444_v3_0/sequence/Ptrichocarpa_444_v3_0.fasta'

# SET OUTPUT #
mpileup_file_name <- 'pop_branch.mpileup'

sh_sort_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/sort_bams.sh'

sh_com_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/pop_samtools_mpileup.sh'

# LOAD PACKAGES #


####################
# NEED TO REMOVE BAD LIBRARIES
bad_inds <- c()
for(bs in bad_samps){
  tmp_ind <- which(meta_0$branch_name == bs)
  bad_inds <- c(bad_inds, tmp_ind)
}

meta <- meta_0[-bad_inds, ]

# MAKE COMMANDS FOR SORTING BAM FILES
bam_files_full <- paste(base_file_dir, sprintf('%03d', meta$root_dir_num), 
  sprintf('%06d', meta$job_num), bam_file_string, sep = '/')

out_sort_files <- paste('branch', meta$branch_name, '_sorted.bam', sep = '')

sort_coms <- paste('samtools sort', bam_files_full, '-o', out_sort_files, 
  sep = ' ')

sh_sort_vec <- c()
sh_sort_vec <- c(sh_sort_vec, '#!/bin/bash')
sh_sort_vec <- c(sh_sort_vec, paste('cd ', bam_sort_dir, sep = ''))
sh_sort_vec <- c(sh_sort_vec, sort_coms)

write.table(sh_sort_vec, file = sh_sort_file, quote = F, sep = '\t', 
  row.names = F, col.names = F)

#####################

sorted_bam_files <- paste(bam_sort_dir, out_sort_files, sep = '/')

bam_files_string <- paste(sorted_bam_files, collapse = ' ')
mpile_com <- paste('samtools mpileup -f', ref_path, '-o', mpileup_file_name, 
  bam_files_string, sep = ' ')

sh_com_vec <- c()
sh_com_vec <- c(sh_com_vec, '#!/bin/bash')
sh_com_vec <- c(sh_com_vec, paste('cd', sam_out_dir, sep = ' '))
sh_com_vec <- c(sh_com_vec, mpile_com)

write.table(sh_com_vec, file = sh_com_file, quote = F, sep = '\t', 
  row.names = F, col.names = F)

quit(save = 'no')

