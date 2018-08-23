# Script to generate figs showing info about unique indels in each poplar
#  library
# Generate following figures:
#  1) Barplots showing the number of unique insertions and deletions in each
#    library
#  2) Histograms of the number of unique insertions and deletions in each
#    library using different distance cutoffs

# LOAD FILES #
## data imported below

# SET VARIABLES #
## info for importing data
file_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/shared_loci/'
file_prefix <- 'popbranch_'
file_suffix <- '.indel_dist_tot'

## characters that are between the sample names in the datafile names
split_chr <- '_v_'

## distances to be used as cutoffs
dist_vec <- c(1000, 2000, 5000, 10000)

# Figure name info
fig_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info/figs/'

combo_bar_pre <- 'pop_branch_Unique_InDel_combo_tally'
combo_bar_file <- paste(fig_dir, combo_bar_pre, '.pdf', sep = '')
combo_bar_png_file <- paste(fig_dir, combo_bar_pre, '.png', sep = '')

uniq_indel_hist_pre <- 'unique_indels_size_hist_'
uniq_del_hist_pre <- 'unique_deletions_size_hist_'

# SET CONSTANTS #

# LOAD LIBRARIES #
library(ggplot2)
library(gridExtra)

###############
# Import data
file_names <- system(paste('ls ', file_dir, '*', file_suffix, sep = ''),
  intern = T)

samp_comps <- gsub(file_suffix, '',
  gsub(paste(file_dir, file_prefix, sep = ''), '', file_names), fixed = T)
all_samps <- sort(unique(unlist(strsplit(samp_comps, split = '_v_'))))

dist_list <- list()
for(smpnm in all_samps){
  file_ind <- grep(smpnm, file_names)
  dist_list[[smpnm]] <- read.table(file_names[file_ind], header = T, 
    stringsAsFactors = F, sep = '\t')
}

# Functions for generating dataframes with showing minimum distance to
#  next-closest indel in any sample
make_min_sizeDF <- function(dist_df, dist_cut){
  dist_df$min_dist <- apply(dist_df[, c(6:ncol(dist_df))], 1, min)
  sub_tab <- dist_df[which(dist_df$min_dist > dist_cut), 
               c(1:5, which(colnames(dist_df) == 'min_dist'))]
}

make_min_dist_indel_tableDF <- function(dist_list, dist_cut){
  sub_dist_list <- lapply(dist_list, make_min_sizeDF, dist_cut = dist_cut)
  sub_dist_tables <- lapply(sub_dist_list, function(x) table(x$INDEL_TYPE))
  sub_dist_df <- data.frame(branch = rep(names(sub_dist_tables), each = 2), 
    num = unlist(sub_dist_tables), type = rep(c('DEL', 'INS'), 
    times = length(sub_dist_tables)), stringsAsFactors = F)
  return(sub_dist_df)
}

# Calculate min distance dataframes for each distance cutoff
uniq_dist_list <- list()
for(i in dist_vec){
  dist_name <- paste((i/1000), 'K', sep = '')
  uniq_dist_list[[dist_name]] <- make_min_dist_indel_tableDF(dist_list, 
    dist_cut = i)
}

# Generate ggplot2 objects for each dataframe
uniq_dist_plot_list <- list()
for(i in seq(length(uniq_dist_list))){
  uniq_dist_plot_list[[i]] <- ggplot(data = uniq_dist_list[[i]], 
    aes(x = branch, y = num, fill = type)) + 
    geom_bar(stat = 'identity', position = position_dodge()) + 
    ggtitle(paste('Number private indels, minimum',  
      names(uniq_dist_list)[i], 'bp\nfrom closest indel in any sample', 
      sep = ' ')) + 
    theme(plot.title = element_text(size = 16))
}

# Generate Barplot showing number of Unique Insertions and Deletions in each 
#  library
pdf(file = combo_bar_file, width = 9, height = 9)
do.call(grid.arrange, c(uniq_dist_plot_list, ncol = 2))
dev.off()

png(filename = combo_bar_png_file, width = 900, height = 900)
do.call(grid.arrange, c(uniq_dist_plot_list, ncol = 2))
dev.off()

# Generate histograms showing the size distribution of the unique insertions
#  and deletions in each library
for(dv in dist_vec){
  dist_cut = dv
  sub_dist_list <- lapply(dist_list, make_min_sizeDF, dist_cut = dist_cut)
  dist_name <- paste((dist_cut/1000), 'K', sep = '')
  info_list <- sub_dist_list
  #
  indel_len_plot_list <- list()
  for(i in seq(length(info_list))){
    indel_len_plot_list[[i]] <- ggplot(data = info_list[[i]]) +
      geom_histogram(aes(x = INDEL_SIZE), fill = 'grey50') +
      scale_x_log10() +
      facet_wrap(~INDEL_TYPE) +
      ggtitle(paste('branch', gsub('_', '.', names(info_list[i])), 
               ' unique indels ', dist_name, 
               'bp\nfrom closest indel in any sample', sep = '')) +
      xlab('log10(InDel Size)') 
  }
  plot_file_name <- paste(fig_dir, uniq_indel_hist_pre, dist_name, 
    '_dist.pdf', sep = '')
  pdf(file = plot_file_name, width = 7, height = 14)
  do.call(grid.arrange, c(indel_len_plot_list, ncol = 2))
  dev.off()
  png_file_name <- paste(fig_dir, uniq_indel_hist_pre, dist_name, 
    '_dist.png', sep = '')
  png(filename = png_file_name, width = 1400, height = 700)
  do.call(grid.arrange, c(indel_len_plot_list, ncol = 5))
  dev.off()
}

# Historgrams for DELETIONS ONLY showing size distribution of unique deletions
for(dv in dist_vec){
  dist_cut = dv
  sub_dist_list <- lapply(dist_list, make_min_sizeDF, dist_cut = dist_cut)
  sub_dist_list_del <- lapply(sub_dist_list, 
    function(x) x[which(x$INDEL_TYPE == 'DEL'), ])
  dist_name <- paste((dist_cut/1000), 'K', sep = '')
  info_list <- sub_dist_list_del
  #
  indel_len_plot_list <- list()
  for(i in seq(length(info_list))){
    indel_len_plot_list[[i]] <- ggplot(data = info_list[[i]]) +
      geom_histogram(aes(x = INDEL_SIZE), fill = 'grey50') +
      scale_x_log10() +
   #   facet_wrap(~INDEL_TYPE) +
      ggtitle(paste('branch', gsub('_', '.', names(info_list[i])), 
               ' unique DELETIONS ', dist_name, 
               'bp\nfrom closest indel in any sample', sep = '')) +
      xlab('log10(Del Size)') 
  }
  plot_file_name <- paste(fig_dir, uniq_del_hist_pre, dist_name, 
    '_dist.pdf', sep = '')
  pdf(file = plot_file_name, width = 7, height = 14)
  do.call(grid.arrange, c(indel_len_plot_list, ncol = 2))
  dev.off()
  png_file_name <- paste(fig_dir, uniq_del_hist_pre, dist_name, 
    '_dist.png', sep = '')
  png(filename = png_file_name, width = 1400, height = 700)
  do.call(grid.arrange, c(indel_len_plot_list, ncol = 5))
  dev.off()

}

quit(save = 'no')
