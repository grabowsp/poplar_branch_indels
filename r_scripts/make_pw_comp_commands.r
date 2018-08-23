# Generate commands to do pairwise comparisons of InDel variation usin
#   vcftools

# IMPORT FILES #
del_vcf_files <- system('ls /home/t4c1/WORK/grabowsk/data/poplar_branches/struc_vcfs/*chr_del.recode.vcf', intern = T)

ins_vcf_files <- system('ls /home/t4c1/WORK/grabowsk/data/poplar_branches/struc_vcfs/*chr_ins.recode.vcf', intern = T)

# OUTPUT FILES #
del_comm_file_out <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/shared_loci/share_diff_comp_commands.sh'

ins_comm_file_out <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/shared_loci/share_diff_INS_comp_commands.sh'

# SET VARIABLES #
comm_str_1 <- 'vcftools --vcf '
comm_str_2 <- ' --diff '
comm_str_3 <- ' --diff-site --out '

in_file_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/struc_vcfs/'
# SET CONSTANTS #

# LOAD LIBRARIES #

############

gen_comp_commands <- function(file_vec, indel){
  samp_names <- gsub('_chr_ins.recode.vcf', '', 
                  gsub('_chr_del.recode.vcf', '', 
                    gsub(paste(in_file_dir, 'branch_', sep = ''), '', file_vec,
                    fixed = T), 
                  fixed = T), 
                fixed = T)
  indel_type <- indel
  comm_vec <- c()
  first_file <- 1
  in_range = T
  while(in_range){
    test_inds <- c((first_file+1):length(file_vec))
    out_pre_vec <- paste('popbranch', samp_names[first_file], 'v', 
      samp_names[test_inds], indel_type, sep = '_')
    tmp_comms <- paste(comm_str_1, file_vec[first_file], comm_str_2, 
      file_vec[test_inds], comm_str_3, out_pre_vec, sep = '')
    comm_vec <- c(comm_vec, tmp_comms)
    first_file <- first_file + 1
    if(first_file == length(file_vec)){in_range = F}
  }
  return(comm_vec)
}

del_comm_vec <- gen_comp_commands(del_vcf_files, indel = 'DEL')
write.table(del_comm_vec, file = del_comm_file_out, quote = F, sep = ' ', 
  row.names = F, col.names = F)

ins_comm_vec <- gen_comp_commands(ins_vcf_files, indel = 'INS')
write.table(ins_comm_vec, file = ins_comm_file_out, quote = F, sep = ' ', 
  row.names = F, col.names = F)

quit(save = 'no')
