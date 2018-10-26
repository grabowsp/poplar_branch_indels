# Code for looking at the VCFs made by SNIFFLES and SURVIVOR. In Oct 2018, I noticed some weird
#  results - InDel sizes looked odd compared to the PacBio pipeline results and the 
#  number of supporting reads wasn't consistent across the pipeline. These results are some
#  preliminary commands to help explore the VCFs

file_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/sniffles/'
combo_file <- 'popbranch.ngmlr.sniffles.survivor.supervisedmerged1kbdist.vcf'
combo_file_full <- paste(file_dir, combo_file, sep = '')

combo_vcf <- read.table(combo_file_full, sep = '\t', stringsAsFactors = F)

info_list <- strsplit(combo_vcf[,8], split = ';')

# gsub('AVGLEN=', '', info_list[[1]][grep('AVGLEN', info_list[[1]])])

avg_length_vec <- as.numeric(unlist(lapply(info_list, function(x) 
  gsub('AVGLEN=', '', x[grep('AVGLEN', x)]))))

combo_info_df <- data.frame(chr = combo_vcf[,1], pos = combo_vcf[,2], type = combo_vcf[,5], 
  avg_len = avg_length_vec, stringsAsFactors = F)

# test <- sapply(combo_vcf[,10], function(x) as.numeric(unlist(strsplit(x, split = ':'))[3]))

for(sn in seq(10)){
  vcf_col <- 9 + sn
  out_col_name <- paste('samp_', sn, '_len', sep = '')
  combo_info_df[,out_col_name] <- sapply(combo_vcf[,vcf_col], function(x) 
    as.numeric(unlist(strsplit(x, split = ':'))[3]))
}

sd_vec <- apply(combo_info_df[, c(5:14)], 1, sd)

###

paxl_file <- 'PAXL.14.5v1.0Ref.ngmlr.sorted.sniffles.vcf'
paxl_file_full <- paste(file_dir, paxl_file, sep = '')
paxl_vcf <- read.table(paxl_file_full, sep = '\t', stringsAsFactors = F)

