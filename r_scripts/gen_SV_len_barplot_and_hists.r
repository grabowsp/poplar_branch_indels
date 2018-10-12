# Script for generating barplots (number of insertions and deletions for each
#  SV calling method) and histograms (size distribution of indels)

# LOAD PACKAGES #
library(ggplot2)
library(gridExtra)

# LOAD DATA #
indel_len_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/combined_analysis/combined_processed_indel_length_info.rds'
indel_len_df <- readRDS(indel_len_file)

# SET CONSTANTS #


# SET VARIABLES #


# SET OUTPUT #
fig_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/combined_analysis/figs/'
indel_bar_file <- 'SV_indel_barplot_DiffCallers.png'

caller_len_hist_file <- 'SV_len_histogram_by_SVcaller.pdf'
samp_call_comp_len_hist_file <- 'SV_len_histogram_withinSamp_SVcaller.pdf'
############
# Barplot of number of INDELS: want insertions and deletions for each approach

indel_len_df$type <- NA
indel_len_df$type[indel_len_df$SVLEN < 0] <- 'DEL'
indel_len_df$type[indel_len_df$SVLEN > 0] <- 'INS'
indel_len_df$typeAndcall <- paste(indel_len_df$type, indel_len_df$SV_caller, 
  sep = '_')

# should remove Indels smaller than 50bp because that was the filter used for
#   SNIFFLES so is not fair comparison if include smaller indels
small_indel_inds <- which(abs(indel_len_df$SVLEN) < 50)
indel_len_df_2 <- indel_len_df[-small_indel_inds, ]

# Need to generate info to make dataframe in proper format to make barplot 
indel_tally_list <- tapply(indel_len_df_2$typeAndcall, indel_len_df_2$samp, 
  table)
tal_vec_list <- lapply(indel_tally_list, function(x) as.vector(x))
tal_names_list <- lapply(indel_tally_list, names)
tal_samp_vec <- c()
for(tls in seq(length(indel_tally_list))){
  tmp_vec <- rep(names(indel_tally_list)[tls], 
    times = length(indel_tally_list[[tls]]))
  tal_samp_vec <- c(tal_samp_vec, tmp_vec)
}

indel_tally_df <- data.frame(tally = unlist(tal_vec_list), 
  type = unlist(tal_names_list), samp = tal_samp_vec, stringsAsFactors = F)
# adjust order of bars in plot
indel_tally_df <- within(indel_tally_df, type <- factor(type, 
  levels = unique(indel_tally_df$type)[c(9,3,10,4,7,1,8,2,11,5,12,6)]))

indel_len_barplot <- ggplot(indel_tally_df, aes(x = type, y = tally, 
  fill = type)) +
  facet_grid(~samp, scales = 'free_x', space = 'free_x') +
  geom_col(position = 'dodge') +
  labs(title = "Number Indels > 50bp by SV caller", 
  x = 'Indel Type x SV caller', y = 'Number') +
  theme(axis.text.x = element_text(angle = 90), 
  plot.title = element_text( hjust = 0.5))

png(file = paste(fig_dir, indel_bar_file, sep = ''), width = 480*3, 
  height = 480)
indel_len_barplot
dev.off()
##########

# Histograms of Indel lengths
indel_len_df$type <- NA
indel_len_df$type[indel_len_df$SVLEN < 0] <- 'DEL'
indel_len_df$type[indel_len_df$SVLEN > 0] <- 'INS'

## Ensure that patterns are consistent across samples
callers <- unique(indel_len_df$SV_caller)
size_cut <- 50

call_hist_list <- list()
for(svc in callers){
  tmp_data <- indel_len_df[which(indel_len_df$SV_caller == svc), ]
  samp_names <- unique(tmp_data$samp)
  samp_list <- list()
  for(svs in samp_names){
    samp_tmp_data_full <- tmp_data[which(tmp_data$samp == svs), ]
    samp_tmp_data <- samp_tmp_data_full[
      which(abs(samp_tmp_data_full$SVLEN) >= size_cut), ]
    tmp_g <- ggplot(data = samp_tmp_data) +
      geom_density(aes(x = abs(SVLEN)), fill = 'grey70') +
      scale_x_log10() +
      facet_wrap(~type) +
      labs(title = paste(svc, svs, sep = ' '), x = 'log10(SV Length)') +
      geom_vline(xintercept = 60, linetype = 'dotted', size = 0.75) +
      geom_vline(xintercept = 200, linetype = 'dotted', size = 0.75) +
      geom_vline(xintercept = 1000, linetype = 'dotted', size = 0.75) +
      geom_vline(xintercept = 4800, linetype = 'dotted', size = 0.75)
    samp_list[[svs]] <- tmp_g
  }
  call_hist_list[[svc]] <- samp_list
}

#test_hist_pdf_1 <- 'test_indel_hist_1.pdf'
pdf(file = paste(fig_dir, caller_len_hist_file, sep = ''), width = 7, 
  height = 10)
for(i in seq(length(call_hist_list))){
  do.call(grid.arrange, c(call_hist_list[[i]], ncol = 2))
}
dev.off()

## Compare patterns within sample across SV callers
tot_samps <- unique(indel_len_df$samp)
size_cut <- 50

samp_hist_list <- list()
for(sbc_samp in tot_samps){
  tmp_data <- indel_len_df[which(indel_len_df$samp == sbc_samp), ]
  caller_names <- unique(tmp_data$SV_caller)
  call_list <- list()
  for(sbc_call in caller_names){
    call_tmp_data_full <- tmp_data[which(tmp_data$SV_caller == sbc_call), ]
    call_tmp_data <- call_tmp_data_full[
      which(abs(call_tmp_data_full$SVLEN) >= size_cut), ]
    tmp_g <- ggplot(data = call_tmp_data) +
      geom_density(aes(x = abs(SVLEN)), fill = 'grey70') +
      scale_x_log10(limits = c(50, 1e5)) +
      facet_wrap(~type) +
      labs(title = paste(sbc_call, sbc_samp, sep = ' '), 
        x = 'log10(SV Length)') +
      geom_vline(xintercept = 60, linetype = 'dotted', size = 0.75) +
      geom_vline(xintercept = 200, linetype = 'dotted', size = 0.75) +
      geom_vline(xintercept = 1000, linetype = 'dotted', size = 0.75) +
      geom_vline(xintercept = 4800, linetype = 'dotted', size = 0.75)
    call_list[[sbc_call]] <- tmp_g
  }
  samp_hist_list[[sbc_samp]] <- call_list
}

pdf(file = paste(fig_dir, samp_call_comp_len_hist_file, sep = ''), width = 7,
  height = 8)
for(i in seq(length(samp_hist_list))){
  do.call(grid.arrange, c(samp_hist_list[[i]], ncol = 2))
}
dev.off()
###########

test_data <- indel_len_df[which(indel_len_df$samp == samp_names[1]),]
test_data_sub1 <- indel_len_df[which(test_data$SV_caller == callers[1]), ]
test_g <- ggplot(data = test_data_sub1) +
  geom_density(aes(x = abs(SVLEN)), fill = 'grey70') +
  scale_x_log10() +
  facet_wrap(~type)

test_hist_file_1 <- 'test_indel_hist_1.png'
png(file = paste(fig_dir, test_hist_file_1, sep = ''))
test_g
dev.off()


quit(save = 'no')

