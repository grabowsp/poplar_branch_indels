# R script for generating vcftools commands to extract the Indel Size and Type
#  info from the poplar branch single-sample files

# LOAD FILES #
meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v1.0.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t')

# SET VARIABLES #
chrom_name_suf <- 'chromosomes_only.recode'

vcf_com_1 <- 'vcftools --vcf '
vcf_com_2 <- ' --get-INFO SVTYPE --get-INFO SVLEN --out '

vcf_info_com_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info/get_chrom_only_indel_info.sh'

# SET CONSTANTS #

# LOAD LIBRARIES #

############3
vcf_names <- gsub('struct', chrom_name_suf, meta$local_file)

out_names <- paste('branch', gsub('.', '_', meta$branch_name, fixed = T), 
               'chrom', sep = '_')

#paste(vcf_com_1, vcf_names, vcf_com_2, out_names, sep = '')

vcf_info_commands <- c(paste('#!/bin/csh'), paste(vcf_com_1, vcf_names, 
                       vcf_com_2, out_names, sep = ''))

write.table(vcf_info_commands, file = vcf_info_com_file, quote = F, 
  sep = ' ', row.names = F, col.names = F)

quit(save = 'no')
