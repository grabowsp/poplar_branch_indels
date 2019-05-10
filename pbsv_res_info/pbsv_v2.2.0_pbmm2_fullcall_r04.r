# Instructions for making file with information that will be used for
#     analyzing results from pbsv. 
#   A list will be saved as an .rds file to be
#     loaded by r scripts and then used to extract the necessary info for each
#     analysis
#   The .rds file will be saved in the same directory and with the same
#     name as the .vcf but with 'data_info.rds' at the end.
#     ex: '/FULL/PATH/TO/DATA/FILE/FILE_NAME.vcf_data_info.rds'

# info_list = list that will contain info required by the R scripts used for
#              analyzing the SV results from vcf files. The elements of the
#              list have names required to be properly inputted for the
#              R scripts 
info_list <- list()

# [[data_dir]] = the full path to the directory containing the data file.
#           ***   This MUST end with a slash in order to work
info_list[['data_dir']] <- '/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/'

# [['vcf_short']] = the name of the vcf file to be analyzed
info_list[['vcf_short']] <- 'PtStettler14.pbmm2.ppsv_v2.2_1.full.call.r04.vcf'

# [['meta_in']] = the full path and file name of the metadata that has info
#                   about the libraries and samples in the analysis. The 
#                   boilerplate includes the most recent metadata file
info_list[['meta_in']] <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'

# [['lib_order']] = a vector of the library names IN THE ORDER THAT THEY ARE
#                     INPUTED into the pbsv when calling genotypes
info_list[['lib_order']] <- c('PAYZ', 'PAZF', 'PAXN', 'PAZG', 'PAXL', 'PBAT', 
  'PAYK', 'PBAW', 'PAZH', 'PBAU')

# [['branch_13_lab']] = a vector of the sample names of the Tree 13 branches
#                       in their order going from base of the tree to the top.
#                       Names should be in numeric form. Boilerplate contains
#                       order for all 5 Tree13 branches
info_list[['branch_13_lab']] <- c(13.4, 13.5, 13.3, 13.2, 13.1)

# [['branch_14_lab']] = a vector of the sample names of the Tree 14 branches
#                       in their order going from base of the tree to the top.
#                       Names should be in numeric form. Boilerplate contains
#                       order for all 5 Tree14 branches
info_list[['branch_14_lab']] <- c(14.5, 14.1, 14.4, 14.3, 14.2)

#####
### Info for BND analysis ###
# [['bnd_receiv_dist']] = the maximum distance between 5' and 3' locations
#                         of BND entries to be considered as possible
#                         BND-INS
info_list[['bnd_receiv_dist']] <- 200

# [['insert_max_size']]] = maximum size of potential INS for a BND-INS
#                          to be included
info_list[['insert_max_size']] <- 1e7

# [['cipos_cut']] = the cutoff for CIPOS (POS confidence interval) for a 
#                   BND-INS to be considered "good"
info_list[['cipos_cut']] <- 200

###########
# rds_out_file = the full path and file name where list will be saved. Should
#                 end with '.rds'. Will just add 'data_info.rds' to the end
#                 of the full vcf file name
rds_out_file <- paste(info_list[['data_dir']], info_list[['vcf_short']], 
  '_data_info.rds', sep = '')

saveRDS(info_list, file = rds_out_file)

quit(save = 'no')

