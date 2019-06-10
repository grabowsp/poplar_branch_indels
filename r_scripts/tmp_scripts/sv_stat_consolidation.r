# code used for consolidating pbsv stat results to make results for Jeremy

ngmlr_indel_r1_file <- '/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r01.vcf_SVstats.txt'

ngmlr_indel_list <- list()

for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, '.vcf_SVstats.txt', sep = '')
  ngmlr_indel_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F) 
}

ngmlr_indel_df <- data.frame(matrix(
  data = unlist(lapply(ngmlr_indel_list, function(x) x[, 2])), 
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(ngmlr_indel_df) <- ngmlr_indel_list[[1]][,1]

pbmm2_indel_r1_file <- '/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.pbmm2.ppsv_v2.2_1.full.call.r01.vcf_SVstats.txt'

pbmm2_indel_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.pbmm2.ppsv_v2.2_1.full.call.r0', i, '.vcf_SVstats.txt', sep = '')
  pbmm2_indel_list[[i]] <- read.table(tmp_name, sep = '\t', header = T, 
    stringsAsFactors = F)
}

pbmm2_indel_df <- data.frame(matrix(
  data = unlist(lapply(pbmm2_indel_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(pbmm2_indel_df) <- pbmm2_indel_list[[1]][,1]

ngmlr_dup_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, '.vcf_DUPstats.txt', sep = '')
  ngmlr_dup_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
}

ngmlr_dup_df <- data.frame(matrix(
  data = unlist(lapply(ngmlr_dup_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(ngmlr_dup_df) <- ngmlr_dup_list[[1]][,1]

pbmm2_dup_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.pbmm2.ppsv_v2.2_1.full.call.r0', i, '.vcf_DUPstats.txt', sep = '')
  pbmm2_dup_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
}

pbmm2_dup_df <- data.frame(matrix(
  data = unlist(lapply(pbmm2_dup_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(pbmm2_dup_df) <- pbmm2_dup_list[[1]][,1]

ngmlr_inv_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, '.vcf_INVstats.txt', sep = '')
  ngmlr_inv_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
}
ngmlr_inv_df <- data.frame(matrix(
  data = unlist(lapply(ngmlr_inv_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(ngmlr_inv_df) <- ngmlr_inv_list[[1]][,1]

pbmm2_inv_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.pbmm2.ppsv_v2.2_1.full.call.r0', i, '.vcf_INVstats.txt', sep = '')
  pbmm2_inv_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
} 
pbmm2_inv_df <- data.frame(matrix(
  data = unlist(lapply(pbmm2_inv_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(pbmm2_inv_df) <- pbmm2_inv_list[[1]][,1]

ngmlr_bnd_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, '.vcf_BNDstats.txt', sep = '')
  ngmlr_bnd_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
}

ngmlr_bnd_df <- data.frame(matrix(
  data = unlist(lapply(ngmlr_bnd_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(ngmlr_bnd_df) <- ngmlr_bnd_list[[1]][,1]

pbmm2_bnd_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.pbmm2.ppsv_v2.2_1.full.call.r0', i, '.vcf_BNDstats.txt', sep = '')
  pbmm2_bnd_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
}
pbmm2_bnd_df <- data.frame(matrix(
  data = unlist(lapply(pbmm2_bnd_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(pbmm2_bnd_df) <- pbmm2_bnd_list[[1]][,1]

ngmlr_perm_indel_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, '.vcf_SVstats_permissive.txt', sep = '')
  ngmlr_perm_indel_list[[i]] <- read.table(tmp_name, sep = '\t', header = T,
    stringsAsFactors = F)
}
ngmlr_perm_indel_df <- data.frame(matrix(
  data = unlist(lapply(ngmlr_perm_indel_list, function(x) x[, 2])),
  nrow = 4, byrow = T), stringsAsFactors = F)
colnames(ngmlr_perm_indel_df) <- ngmlr_perm_indel_list[[1]][,1]

########

ngmlr_dup_geno_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, '.vcf_goodDUP_genos.rds', sep = '')
  ngmlr_dup_geno_list[[i]] <- readRDS(tmp_name)
}

ngmlr_indel_geno_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, '.vcf_goodINDEL_genos.rds', sep = '')
  ngmlr_indel_geno_list[[i]] <- readRDS(tmp_name)
}

indel_name_tab <- table(unlist(
  lapply(ngmlr_indel_geno_list, function(x) rownames(x))))

which(indel_name_tab == 4)

ngmlr_perm_indel_geno_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i, '.vcf_goodINDEL_genos_permissive.rds', 
  sep = '')
  ngmlr_perm_indel_geno_list[[i]] <- readRDS(tmp_name)
}

perm_indel_name_tab <- table(unlist(
  lapply(ngmlr_perm_indel_geno_list, function(x) rownames(x))))
which(perm_indel_name_tab == 4)
