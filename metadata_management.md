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
write.table(meta, file = meta_out_file, quote = F, sep = '\t', row.names = F, 
  col.names = T)

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

## Add Sequencing Stats from PacBio Sequencer Metadata
```
meta_in <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v3.0.txt'
samp_meta <- read.table(meta_in, header = T, stringsAsFactors = F, sep = '\t')

seq_meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/pop_branch_sequencer_info.tsv'

seq_meta <- read.table(seq_meta_file, header = T, stringsAsFactors = F, 
  sep = '\t')

samp_meta[, c('tot_reads', 'tot_bp', 'singlepass_reads', 'singlepass_bp', 
  'X20kb_reads', 'X20kb_bp')] <- NA

for(i in seq(nrow(samp_meta))){
  test_lib <- samp_meta$lib_name[i]
  seq_inds <- grep(test_lib, seq_meta$Description)
  samp_meta$tot_bp[i] <- sum(as.numeric(gsub(',', '', seq_meta[seq_inds, 
    'Total.Basepairs'])))
  samp_meta$tot_reads[i] <- sum(as.numeric(gsub(',', '', seq_meta[seq_inds,    
    'Total.Reads'])))
  samp_meta$singlepass_bp[i] <- sum(as.numeric(gsub(',', '', 
    seq_meta[seq_inds, 'Single.Pass.Basepairs'])))
  samp_meta$singlepass_reads[i] <- sum(as.numeric(gsub(',', '', 
    seq_meta[seq_inds, 'Single.Pass.Reads'])))
  samp_meta$X20kb_bp[i] <- sum(as.numeric(gsub(',', '', seq_meta[seq_inds,
    'X20kb.Basepairs'])))
  samp_meta$X20kb_reads[i] <- sum(as.numeric(gsub(',', '', seq_meta[seq_inds,
    'X20kb.Reads'])))
}

meta_out_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v4.0.txt'
write.table(samp_meta, file = meta_out_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)
```

