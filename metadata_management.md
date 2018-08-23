# Metadata Management
## Add Local VCF File Names
R code
```
meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/Poplar_Branch_Metadata.tsv'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t')

vcf_files <- system('ls /home/t4c1/WORK/grabowsk/data/poplar_branches/struc_vcfs/*.vcf', intern = T)

meta$local_file <- NA
for(i in c(1:nrow(meta))){
  file_ind <- grep(gsub('.', '_', meta$branch_name[i], fixed = T), vcf_files)
  meta$local_file[i] <- vcf_files[file_ind]
}

meta_out_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v1.0.txt'
write.table(meta, file = meta_out_file, quote = F, sep = '\t', row.names = F, col.names = T)

```
## Updated Library Info - Aug. 23 2018
Updated VCFs for 14.2(PAXN), 14.1(PAYK), and 14.3(PAXL)
R code:
```
meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v1.0.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t')

meta$job_num[which(meta$lib_name == 'PAXN')] <- 2243
meta$job_num[which(meta$lib_name == 'PAYK')] <- 2245
meta$job_num[which(meta$lib_name == 'PAXL')] <- 2269

meta_out_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v1.0.txt'
write.table(meta, file = meta_out_file, quote = F, sep = '\t', row.names = F, col.names = T)

``` 
## Add Sequencing Output Info
### Calculate Seq Depth with vcftools
#### shell script
```
bash /home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/calc_seq_depth.sh
```
### Import Depth Info
R script:
```
meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v1.0.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t')

meta$n_sites <- NA
meta$mean_depth <- NA

for(i in seq(nrow(meta))){
  depth_file <- gsub('_struct.vcf', '.idepth', meta$local_file[i])
  tmp_d_info <- read.table(depth_file, header = T, stringsAsFactors = F, sep = '\t')
  meta$n_sites[i] <- tmp_d_info$N_SITES[1]
  meta$mean_depth[i] <- tmp_d_info$MEAN_DEPTH[1]
}

meta_out_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v2.0.txt'
write.table(meta, file = meta_out_file, quote = F, sep = '\t', row.names = F, col.names = T)

```
### Import Tally of Insertions and Deletions
R script:
```
meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v2.0.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t')

meta$n_del <- NA
meta$n_ins <- NA

for(i in seq(nrow(meta))){
  inds_file <- gsub('_struct.vcf', '_del_positions.txt', meta$local_file[i])
  tmp_wc_command <- paste('wc -l ', inds_file, sep = '')
  tmp_n_del_long <- system(tmp_wc_command, intern = T)
  tmp_n_del <- as.numeric(unlist(strsplit(tmp_n_del_long, split = ' '))[1])
  meta$n_del[i] <- tmp_n_del
}

for(i in seq(nrow(meta))){
  inds_file <- gsub('_struct.vcf', '_ins_positions.txt', meta$local_file[i])
  tmp_wc_command <- paste('wc -l ', inds_file, sep = '')
  tmp_n_ins_long <- system(tmp_wc_command, intern = T)
  tmp_n_ins <- as.numeric(unlist(strsplit(tmp_n_ins_long, split = ' '))[1])
  meta$n_ins[i] <- tmp_n_ins
}

meta_out_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v3.0.txt'
write.table(meta, file = meta_out_file, quote = F, sep = '\t', row.names = F, col.names = T)

```


