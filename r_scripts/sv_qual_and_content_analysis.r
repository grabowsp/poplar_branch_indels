# Script for getting the names of the SV's used for analysis of large SVs

# uses the results from calculating statistics

library(ggplot2)

annot_analysis_function_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/annot_analysis_functions.r'
source(annot_analysis_function_file)

annot_file <- '/home/t4c1/WORK/sujan/Ptricocarpa/PtrichocarpaStettler14_532_v1.1/annotation/PtrichocarpaStettler14_532_v1.1.gene.gff3'

annot <- read.table(annot_file, header = F, stringsAsFactors = F)

annot_gene_inds <- which(annot[,3] == 'gene')

gene_annot <- annot[annot_gene_inds, ]

annot_cds_inds <- which(annot[,3] == 'CDS')

cds_annot <- annot[annot_cds_inds, ]

tdf_file <- paste('/home/f1p1/tmp/poplar_branches/ref_stuff/', 
  'poplar_var_14.5_V1_chromosomes/poplar_14_5_v1_tandemrepeat.bed', sep = '')

tdf_bed <- read.table(tdf_file, header = F, sep = '\t', stringsAsFactors = F) 

tdf_annot <- data.frame(
  chr = tdf_bed[,1], platform = 'tdf', annot = 'tandem_repeat', 
  start = tdf_bed[,2], end = tdf_bed[,3], dot_1 = '.', strand = tdf_bed[,6],
  dot_2 = '.', name = tdf_bed[,4], stringsAsFactors = F  
)

all_indel_names_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/',
    'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i,
    '.8branch.vcf_filtered_INDEL_names.rds', sep = '')
  all_indel_names_list[[i]] <- readRDS(tmp_name)
}

all_indel_df_list <- list()
for(i in seq(length(all_indel_names_list))){
  tmp_split <- strsplit(all_indel_names_list[[i]], split = '_')
  tmp_chr <- unlist(lapply(tmp_split, function(x) x[1]))
  tmp_pos <- as.numeric(unlist(lapply(tmp_split, function(x) x[2])))
  tmp_type <- unlist(lapply(tmp_split, function(x) x[3]))
  tmp_size <- as.numeric(unlist(lapply(tmp_split, function(x) x[4])))
  tmp_df <- data.frame(name = all_indel_names_list[[i]], chr = tmp_chr, 
    pos = tmp_pos, type = tmp_type, size = tmp_size, stringsAsFactors = F)
  all_indel_df_list[[i]] <- tmp_df
}

big_indel_name_list <- lapply(all_indel_df_list, function(x) 
  x$name[which(x$size > 50000)])

table(unlist(big_indel_name_list))
# The same 11 come up in all 4 reps, though some have slightly different
#  names

big_indel_inds <- which(all_indel_df_list[[1]]$size > 50000)

# Load duplicate info

all_dup_names_list <- list()
for(i in seq(4)){
  tmp_name <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/',
    'PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r0', i,
    '.8branch.vcf_filtered_DUP_names.rds', sep = '')
  all_dup_names_list[[i]] <- readRDS(tmp_name)
}

all_dup_df_list <- list()
for(i in seq(length(all_dup_names_list))){
  tmp_split <- strsplit(all_dup_names_list[[i]], split = '_')
  tmp_chr <- unlist(lapply(tmp_split, function(x) x[1]))
  tmp_pos <- as.numeric(unlist(lapply(tmp_split, function(x) x[2])))
  tmp_type <- unlist(lapply(tmp_split, function(x) x[3]))
  tmp_size <- as.numeric(unlist(lapply(tmp_split, function(x) x[4])))
  tmp_df <- data.frame(name = all_dup_names_list[[i]], chr = tmp_chr,
    pos = tmp_pos, type = tmp_type, size = tmp_size, stringsAsFactors = F)
  tmp_df <- tmp_df[-which(is.na(tmp_df$size)), ]
  all_dup_df_list[[i]] <- tmp_df
}

big_dup_name_list <- lapply(all_dup_df_list, function(x)
  x$name[which(x$size > 50000)])

table(unlist(big_dup_name_list))
# same 4 in all 4 reps

#################
# Caclculate percentage of genic sequence in different sized windows and
#   in all Indels

test_vcf <- all_indel_df_list[[1]]

genic_window_50k <- calc_genomewide_perc(annot_df = gene_annot, 
  window_size = 50000, loud = F)

genic_window_100k <- calc_genomewide_perc(annot_df = gene_annot, 
  window_size = 100000, loud = F)

genic_window_10k <- calc_genomewide_perc(annot_df = gene_annot, 
  window_size = 10000, loud = F)

window_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplare_genic_window_calcs.rda', sep = '/')

save(genic_window_50k, genic_window_100k, genic_window_10k, 
  file = window_calc_out)

all_indel_calcs <- sapply(test_vcf$name, calc_sv_perc_code, 
  annot_df = gene_annot)

indel_genic_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs', 
  'poplar_indel_per_genic_seq.rds', sep = '/')

saveRDS(all_indel_calcs, file = indel_genic_calc_out)

dup_vcf <- all_dup_df_list[[1]]
small_dup_inds <- which(dup_vcf$size < 20)
dup_vcf_1 <- dup_vcf[-small_dup_inds, ]

all_dup_calcs <- sapply(dup_vcf_1$name, calc_sv_perc_code,
  annot_df = gene_annot)

dup_genic_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_genic_seq.rds', sep = '/')

saveRDS(all_dup_calcs, file = dup_genic_calc_out)

######################
# Calculate percentage of CDS sequence in windows and SVs

test_vcf <- all_indel_df_list[[1]]

cds_indel_calcs <- sapply(test_vcf$name, calc_sv_perc_code,
  annot_df = cds_annot)

indel_cds_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_cds_seq.rds', sep = '/')

saveRDS(cds_indel_calcs, file = indel_cds_calc_out)

dup_vcf <- all_dup_df_list[[1]]
small_dup_inds <- which(dup_vcf$size < 20)
dup_vcf_1 <- dup_vcf[-small_dup_inds, ]

cds_dup_calcs <- sapply(dup_vcf_1$name, calc_sv_perc_code,
  annot_df = cds_annot)

dup_cds_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_cds_seq.rds', sep = '/')

saveRDS(cds_dup_calcs, file = dup_cds_calc_out)

cds_window_50k <- calc_genomewide_perc(annot_df = cds_annot,
  window_size = 50000, loud = T)

cds_window_100k <- calc_genomewide_perc(annot_df = cds_annot,
  window_size = 100000, loud = T)

cds_window_10k <- calc_genomewide_perc(annot_df = cds_annot,
  window_size = 10000, loud = T)

cds_window_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_cds_window_calcs.rda', sep = '/')

save(cds_window_50k, cds_window_100k, cds_window_10k,
  file = cds_window_calc_out)

#####
# Calculate percentage of tandem repeat sequence for windows and SVs as
#   a measure of repetitive sequence context

tr_window_50k <- calc_genomewide_perc(annot = tdf_annot,
  window_size = 50000, loud = T)

tr_window_100k <- calc_genomewide_perc(annot = tdf_annot,
  window_size = 100000, loud = T)

tr_window_10k <- calc_genomewide_perc(annot = tdf_annot,
  window_size = 10000, loud = T)

tr_window_file_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_tr_window_calcs.rda', sep = '/')

save(tr_window_50k, tr_window_100k, tr_window_10k, file = tr_window_file_out)

indel_tr_calcs <- sapply(test_vcf$name, calc_sv_perc_code,
  annot_df = tdf_annot)

indel_tr_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_tr_seq.rds', sep = '/')

saveRDS(indel_tr_calcs, file = indel_tr_calc_out)

dup_tr_calcs <- sapply(dup_vcf_1$name, calc_sv_perc_code,
  annot_df = tdf_annot)

dup_tr_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_tr_seq.rds', sep = '/')

saveRDS(dup_tr_calcs, file = dup_tr_calc_out)


##########
# Generate genic sequence boxplot for SVs

indel_genic_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_genic_seq.rds', sep = '/')

all_indel_calcs <- readRDS(indel_genic_calc_out)

dup_genic_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_genic_seq.rds', sep = '/')

all_dup_calcs <- readRDS(dup_genic_calc_out)

sv_genic_calcs <- c(all_indel_calcs, all_dup_calcs)

sv_name_info <- strsplit(names(sv_genic_calcs), split = '_')
sv_type_vec <- unlist(lapply(sv_name_info, function(x) x[3]))
sv_size_vec <- as.numeric(unlist(lapply(sv_name_info, function(x) x[4])))

#del_inds <- which(sv_type_vec == 'DEL')
#ins_inds <- which(sv_type_vec == 'INS')

inds_50k <- which(sv_size_vec > 50000)
inds_10k <- setdiff(which(sv_size_vec > 10000), inds_50k)
inds_1k <- setdiff(which(sv_size_vec > 1000), c(inds_50k, inds_10k))
inds_20 <- setdiff(which(sv_size_vec >= 20), c(inds_50k, inds_10k, inds_1k))

sv_perc_genic_df <- data.frame(name = names(sv_genic_calcs),
  perc_genic = sv_genic_calcs, type = NA, size = NA, size_class = NA,
  tot_class = NA, stringsAsFactors = F)

sv_perc_genic_df$type <- sv_type_vec
sv_perc_genic_df$size <- sv_size_vec

sv_perc_genic_df$size_class[inds_50k] <- '>50kbp'
sv_perc_genic_df$size_class[inds_10k] <- '10kbp-50kbp'
sv_perc_genic_df$size_class[inds_1k] <- '1kbp-10kbp'
sv_perc_genic_df$size_class[inds_20] <- '20bp-1kbp'

sv_perc_genic_df$tot_class <- paste(sv_perc_genic_df$type, 
  sv_perc_genic_df$size_class, sep = ' ')

# add tandem repeat and CDS info
indel_tr_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_tr_seq.rds', sep = '/')

indel_tr_calcs <- readRDS(indel_tr_calc_out)

dup_tr_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_tr_seq.rds', sep = '/')

dup_tr_calcs <- readRDS(dup_tr_calc_out)

sv_tr_calcs <- c(indel_tr_calcs, dup_tr_calcs)

sv_perc_genic_df$perc_tr <- sv_tr_calcs

indel_cds_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_cds_seq.rds', sep = '/')

indel_cds_calcs <- readRDS(indel_cds_calc_out)

dup_cds_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_cds_seq.rds', sep = '/')

dup_cds_calcs <- readRDS(dup_cds_calc_out)

sv_cds_calcs <- c(indel_cds_calcs, dup_cds_calcs)

sv_perc_genic_df$perc_cds <- sv_cds_calcs

# remove insertion indices because info I have doesn't actually look at
#  insertions
ins_inds <- which(sv_perc_genic_df$type == 'INS')
sv_perc_genic_df_2 <- sv_perc_genic_df[-ins_inds, ]

########
# generate window dataframe

window_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplare_genic_window_calcs.rda', sep = '/')

load(window_calc_out)

window_10k_vals <- unlist(genic_window_10k)

window_calc_df <- data.frame(
  name = paste('window', seq(length(window_10k_vals)), sep = '_'),
  perc_genic = window_10k_vals, type = '10k window', size = NA, 
  size_class = NA, tot_class = '10k windows', stringsAsFactors = F)

tr_window_file_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_tr_window_calcs.rda', sep = '/')
load(tr_window_file_out)

window_calc_df$perc_tr <- sample(unlist(tr_window_10k), nrow(window_calc_df))

cds_window_file_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_cds_window_calcs.rda', sep = '/')

load(cds_window_file_out)

window_calc_df$perc_cds <- unlist(cds_window_10k)

# combine SV and window data.frames

all_perc_genic_df <- rbind(sv_perc_genic_df_2, window_calc_df)

all_perc_genic_df$tot_class_2 <- NA
for(tc in unique(all_perc_genic_df$tot_class)){
  tc_inds <- which(all_perc_genic_df$tot_class == tc)
  tmp_num <- length(tc_inds)
  add_text <- paste('(n= ', tmp_num, ')', sep = '')
  new_text <- paste(tc, add_text, sep = ' ')
  all_perc_genic_df$tot_class_2[tc_inds] <- new_text
}

tot_class_levels <- unique(all_perc_genic_df$tot_class_2)[
  c(9,4,8,3,7,2,6,1,5)]

all_perc_genic_df$tot_class_2 <- factor(all_perc_genic_df$tot_class_2,
  levels = tot_class_levels)

bp_colors <- rep(NA, times = length(tot_class_levels))
bp_colors[grep('window', tot_class_levels)] <- 'gray70'
bp_colors[grep('DEL', tot_class_levels)] <- 'red2'
bp_colors[grep('DUP', tot_class_levels)] <- 'green2'
# bp_colors[grep('INS', tot_class_levels)] <- 'blue2'

all_boxplot_p <- ggplot(all_perc_genic_df,
  aes(x = tot_class_2, y = perc_genic)) +
  geom_boxplot(fill = bp_colors) +
#  scale_fill_manual(values = bp_colors) +
  stat_summary(fun.y=mean, geom='point', shape = 23, size=4,
    position = position_dodge(0.9), fill = bp_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Percent genic sequence') +
  xlab('SV type and size class') +
  ggtitle('Percent genic sequence within SVs of different sizes')

all_boxplot_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_SV_Percent_Genic_boxplot_with_windows.pdf', sep = '/')

pdf(file = all_boxplot_out)
all_boxplot_p
dev.off()

##########
# tandem repeat boxplot
tr_boxplot_p <- ggplot(all_perc_genic_df,
  aes(x = tot_class_2, y = perc_tr)) +
  geom_boxplot(fill = bp_colors) +
#  scale_fill_manual(values = bp_colors) +
  stat_summary(fun.y=mean, geom='point', shape = 23, size=4,
    position = position_dodge(0.9), fill = bp_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Percent tandem repeat sequence') +
  xlab('SV type and size class') +
  ggtitle('Percent tandem repeat sequence within SVs of different sizes')

tr_boxplot_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_SV_Percent_TandemRepeat_boxplot_with_windows.pdf', sep = '/')

pdf(file = tr_boxplot_out)
tr_boxplot_p
dev.off()

##########
# CDS Boxplot

cds_boxplot_p <- ggplot(all_perc_genic_df,
  aes(x = tot_class_2, y = perc_cds)) +
  geom_boxplot(fill = bp_colors) +
#  scale_fill_manual(values = bp_colors) +
  stat_summary(fun.y=mean, geom='point', shape = 23, size=4,
    position = position_dodge(0.9), fill = bp_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Percent coding sequence') +
  xlab('SV type and size class') +
  ggtitle('Percent coding sequence within SVs of different sizes')

cds_boxplot_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_SV_Percent_CDS_boxplot_with_windows.pdf', sep = '/')

pdf(file = cds_boxplot_out)
cds_boxplot_p
dev.off()


# windows = gray
# DEL = red
# DUP = green
# INS = blue

######
# Test for linear relationship between % genic, tandem repeat, or coding seq
#   and SV size

indel_genic_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_genic_seq.rds', sep = '/')

all_indel_calcs <- readRDS(indel_genic_calc_out)

dup_genic_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_genic_seq.rds', sep = '/')

all_dup_calcs <- readRDS(dup_genic_calc_out)

sv_genic_calcs <- c(all_indel_calcs, all_dup_calcs)

sv_name_info <- strsplit(names(sv_genic_calcs), split = '_')
sv_type_vec <- unlist(lapply(sv_name_info, function(x) x[3]))
sv_size_vec <- as.numeric(unlist(lapply(sv_name_info, function(x) x[4])))

sv_perc_genic_df <- data.frame(name = names(sv_genic_calcs),
  perc_genic = sv_genic_calcs, type = NA, size = NA,
  stringsAsFactors = F)

sv_perc_genic_df$type <- sv_type_vec
sv_perc_genic_df$size <- sv_size_vec

indel_tr_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_tr_seq.rds', sep = '/')

indel_tr_calcs <- readRDS(indel_tr_calc_out)

dup_tr_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_tr_seq.rds', sep = '/')

dup_tr_calcs <- readRDS(dup_tr_calc_out)

sv_tr_calcs <- c(indel_tr_calcs, dup_tr_calcs)

sv_perc_genic_df$perc_tr <- sv_tr_calcs

indel_cds_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_cds_seq.rds', sep = '/')

indel_cds_calcs <- readRDS(indel_cds_calc_out)

dup_cds_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_cds_seq.rds', sep = '/')

dup_cds_calcs <- readRDS(dup_cds_calc_out)

sv_cds_calcs <- c(indel_cds_calcs, dup_cds_calcs)

sv_perc_genic_df$perc_cds <- sv_cds_calcs

# remove insertion indices because info I have doesn't actually look at
#  insertions
ins_inds <- which(sv_perc_genic_df$type == 'INS')
sv_perc_genic_df_2 <- sv_perc_genic_df[-ins_inds, ]

del_df <- sv_perc_genic_df_2[which(sv_perc_genic_df_2$type == 'DEL'), ]
dup_df <- sv_perc_genic_df_2[which(sv_perc_genic_df_2$type == 'DUP'), ]

genic_v_size_del <- lm(del_df$perc_genic ~ log(del_df$size))
summary(genic_v_size_del)
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       0.251572   0.012072  20.840  < 2e-16 ***
# log(del_df$size) -0.019643   0.002482  -7.914 2.74e-15 ***
# Residual standard error: 0.3613 on 10431 degrees of freedom
# Multiple R-squared:  0.005969,	Adjusted R-squared:  0.005873 
# F-statistic: 62.63 on 1 and 10431 DF,  p-value: 2.742e-15
#### Significant

genic_v_size_dup <- lm(dup_df$perc_genic ~ log(dup_df$size))
summary(genic_v_size_dup)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       0.301510   0.050401   5.982  3.7e-09 ***
# log(dup_df$size) -0.008192   0.009216  -0.889    0.374
#### Not significant

tr_v_size_del <- lm(del_df$perc_tr ~ log(del_df$size))
summary(tr_v_size_del)
#               Estimate Std. Error t value Pr(>|t|)    
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       0.40444    0.01274   31.74   <2e-16 ***
# log(del_df$size) -0.04045    0.00262  -15.44   <2e-16 ***
# Residual standard error: 0.3814 on 10431 degrees of freedom
# Multiple R-squared:  0.02234,	Adjusted R-squared:  0.02225 
# F-statistic: 238.3 on 1 and 10431 DF,  p-value: < 2.2e-16
#### Significant

tr_v_size_dup <- lm(dup_df$perc_tr ~ log(dup_df$size))
summary(tr_v_size_dup)
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       0.079006   0.021927   3.603 0.000339 ***
# log(dup_df$size) -0.005984   0.004009  -1.493 0.136068
#### Not significant

cds_v_size_del <- lm(del_df$perc_cds ~ log(del_df$size))
summary(cds_v_size_del)
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)       0.0164225  0.0036624   4.484  7.4e-06 ***
# log(del_df$size) -0.0004056  0.0007530  -0.539     0.59 
#### Not significant

cds_v_size_dup <- lm(dup_df$perc_cds ~ log(dup_df$size))
summary(cds_v_size_dup)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -0.032773   0.016212  -2.021   0.0437 *  
# log(dup_df$size)  0.012684   0.002964   4.279 2.17e-05 ***
# Residual standard error: 0.1339 on 628 degrees of freedom
# Multiple R-squared:  0.02833,	Adjusted R-squared:  0.02678 
# F-statistic: 18.31 on 1 and 628 DF,  p-value: 2.173e-05
### Significant, but in a weird way - I didn't expect a positive relationship

###############
# Randomly select SVs for manual inspection

test_vcf <- all_indel_df_list[[1]]
dup_vcf <- all_dup_df_list[[1]]
small_dup_inds <- which(dup_vcf$size < 20)
dup_vcf_1 <- dup_vcf[-small_dup_inds, ]

sv_vcf <- rbind(test_vcf, dup_vcf_1)

sv_inds_50k <- which(sv_vcf$size > 50000)
sv_inds_10k <- setdiff(which(sv_vcf$size > 10000), sv_inds_50k)
sv_inds_1k <- setdiff(which(sv_vcf$size > 1000), c(sv_inds_50k, sv_inds_10k))
sv_inds_20 <- setdiff(which(sv_vcf$size >= 20), 
  c(sv_inds_50k, sv_inds_10k, sv_inds_1k))

sv_vcf$size_class <- NA
sv_vcf$size_class[sv_inds_50k] <- '>50kb'
sv_vcf$size_class[sv_inds_10k] <- '10kb-50kb'
sv_vcf$size_class[sv_inds_1k] <- '1kb-10kb'
sv_vcf$size_class[sv_inds_20] <- '20bp-1kb'

manual_check_inds <- c()

for(svt in c('DEL', 'DUP')){
  tmp_inds <- intersect(which(sv_vcf$type == svt), sv_inds_50k)
  manual_check_inds <- c(manual_check_inds, tmp_inds)
}

for(svt in c('DEL', 'DUP')){
  tmp_inds <- sample(intersect(which(sv_vcf$type == svt), sv_inds_10k), 10)
  manual_check_inds <- c(manual_check_inds, tmp_inds)
}

for(svt in c('DEL', 'DUP', 'INS')){
  tmp_inds_1 <- sample(intersect(which(sv_vcf$type == svt), sv_inds_1k), 10)
  tmp_inds_2 <- sample(intersect(which(sv_vcf$type == svt), sv_inds_20), 10)
  manual_check_inds <- c(manual_check_inds, tmp_inds_1, tmp_inds_2)
}

manual_check_df <- sv_vcf[manual_check_inds, ]

manual_check_df_2 <- manual_check_df[order(manual_check_df$type), ]

man_check_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_SVs_to_manualy_inspect.txt', sep = '/')

write.table(manual_check_df_2, file = man_check_out, quote = F, sep = '\t',
  row.names = F, col.names = T)


man_check_res_in <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_sv_man_inspection_results.tsv', sep = '/')

man_check_res <- read.delim(man_check_res_in, header = T, sep = '\t',
  stringsAsFactors = F)

# Tally up the total number of each quality class, then broken down by
#  type, size, and type:size

table(man_check_res$support)
# moderate   strong     weak 
#      14       69       12

table(paste(man_check_res$type, man_check_res$support, sep = ' '))
# DEL moderate   DEL strong     DEL weak DUP moderate   DUP strong     DUP weak 
#            8           30            3            6           19            9 
#  INS strong 
#          20

table(paste(man_check_res$size_class, man_check_res$support, sep = ' '))
# 10kb-50kb moderate   10kb-50kb strong     10kb-50kb weak  1kb-10kb moderate 
#                  4                  8                  8                  1 
#   1kb-10kb strong      1kb-10kb weak  20bp-1kb moderate    20bp-1kb strong 
#                27                  2                  2                 28 
#    >50kb moderate       >50kb strong         >50kb weak 
#                 7                  6                  2

table(paste(man_check_res$type, man_check_res$size_class, 
  man_check_res$support, sep = ' '))
# DEL 10kb-50kb moderate   DEL 10kb-50kb strong     DEL 10kb-50kb weak 
#                     3                      5                      2 
#   DEL 1kb-10kb strong    DEL 20bp-1kb strong     DEL >50kb moderate 
#                    10                     10                      5 
#      DEL >50kb strong         DEL >50kb weak DUP 10kb-50kb moderate 
#                     5                      1                      1 
#  DUP 10kb-50kb strong     DUP 10kb-50kb weak  DUP 1kb-10kb moderate 
#                     3                      6                      1 
#   DUP 1kb-10kb strong      DUP 1kb-10kb weak  DUP 20bp-1kb moderate 
#                     7                      2                      2 
#   DUP 20bp-1kb strong     DUP >50kb moderate       DUP >50kb strong 
#                     8                      2                      1 
#        DUP >50kb weak    INS 1kb-10kb strong    INS 20bp-1kb strong 
#                     1                     10                     10


##################
# Test differences in mean and distributions of genic, cds, and tr for SVs
#   vs windows

window_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplare_genic_window_calcs.rda', sep = '/')

load(window_calc_out)

genic_window_10k_vals <- unlist(genic_window_10k)

indel_genic_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_genic_seq.rds', sep = '/')

all_indel_calcs <- readRDS(indel_genic_calc_out)

dup_genic_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_genic_seq.rds', sep = '/')

all_dup_calcs <- readRDS(dup_genic_calc_out)

sv_genic_calcs <- c(all_indel_calcs, all_dup_calcs)

ks.test(x = sv_genic_calcs[grep('DEL', names(sv_genic_calcs))], 
  y = genic_window_10k_vals, alternative = 'less')
# D^- = 0.13608, p-value < 2.2e-16
# alternative hypothesis: the CDF of x lies below that of y

ks.test(x = sv_genic_calcs[grep('DEL', names(sv_genic_calcs))],
  y = genic_window_10k_vals, alternative = 'less')
# D^- = 0.20273, p-value < 2.2e-16
# alternative hypothesis: the CDF of x lies below that of y

#############

tr_window_file_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_tr_window_calcs.rda', sep = '/')

load(tr_window_file_out)

tr_window_10k_vals <- unlist(tr_window_10k)

indel_tr_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_tr_seq.rds', sep = '/')

indel_tr_calcs <- readRDS(indel_tr_calc_out)

dup_tr_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_tr_seq.rds', sep = '/')

dup_tr_calcs <- readRDS(dup_tr_calc_out)

sv_tr_calcs <- c(indel_tr_calcs, dup_tr_calcs)

ks.test(x = sv_tr_calcs[grep('DEL', names(sv_tr_calcs))],
  y = tr_window_10k_vals, alternative = 'greater')
# D^+ = 0.58139, p-value < 2.2e-16
# alternative hypothesis: the CDF of x lies above that of y

ks.test(x = sv_tr_calcs[grep('DUP', names(sv_tr_calcs))],
  y = tr_window_10k_vals, alternative = 'greater')
# D^+ = 0.71785, p-value < 2.2e-16
# alternative hypothesis: the CDF of x lies above that of y

##############3

cds_window_file_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_cds_window_calcs.rda', sep = '/')

load(cds_window_file_out)

cds_window_10k_vals <- unlist(cds_window_10k)

indel_cds_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_indel_per_cds_seq.rds', sep = '/')

indel_cds_calcs <- readRDS(indel_cds_calc_out)

dup_cds_calc_out <- paste('/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs',
  'poplar_duplicate_per_cds_seq.rds', sep = '/')

dup_cds_calcs <- readRDS(dup_cds_calc_out)

sv_cds_calcs <- c(indel_cds_calcs, dup_cds_calcs)

ks.test(x = sv_cds_calcs[grep('DEL', names(sv_cds_calcs))],
  y = cds_window_10k_vals, alternative = 'less')
# D^- = 0.011553, p-value = 0.1109
# alternative hypothesis: the CDF of x lies below that of y

ks.test(x = sv_cds_calcs[grep('DUP', names(sv_cds_calcs))],
  y = cds_window_10k_vals, alternative = 'less')
#D^- = 0.019, p-value = 0.6391
#alternative hypothesis: the CDF of x lies below that of y






