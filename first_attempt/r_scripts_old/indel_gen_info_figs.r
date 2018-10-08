# Script to generate figures with info about number and size distribution
#  of indels in each library
#  Generate following figures
#  1) Barplot showing number of Insertions and Deletions in each library
#  2) Histogram/density plots showing size distributions of Insertions
#      and Deletions in each library

# LOAD FILES #
meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v1.0.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t')

# SET VARIABLES #
bar_fig_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info/figs/pop_branch_InDel_tally.png'

combo_len_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info/figs/pop_branch_size_hists_combo.pdf'

# SET CONSTANTS #

# LOAD LIBRARIES #
library(ggplot2)
library(gridExtra)

##################
## GENERATE OBJECT WITH INDEL INFO FOR ALL SAMPLES
# Re-constitute file names to import each file
adj_branch_names <- gsub('.', '_', meta$branch_name, fixed = T)
info_files_pre <- paste('branch_', adj_branch_names, '.INFO', sep = '')
file_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/indel_info/'
info_files <- paste(file_dir, info_files_pre, sep = '')

# generate list containing info for each sample
info_list <- list()
for(i in seq(length(info_files))){
  info_list[[adj_branch_names[i]]] <- read.table(info_files[i], header = T, 
    stringsAsFactors = F, sep = '\t')
}

## BARPLOT SHOWING NUMBER OF INSERTIONS AND DELETIONS IN EACH FILE
info_type_tables <- lapply(info_list, function(x) table(x$SVTYPE))

info_type_df <- data.frame(branch = rep(names(info_type_tables), each = 2), num = unlist(info_type_tables), type = rep(c('DEL', 'INS'), times = length(info_type_tables)), stringsAsFactors = F)

type_gg <- ggplot(data = info_type_df, aes(x = branch, y = num, fill = type)) + geom_bar(stat = 'identity', position = position_dodge())

png(bar_fig_file)
type_gg
dev.off()

## HISTOGRAMS OF INDEL SIZES IN EACH SAMPLE
indel_len_plot_list <- list()
for(i in seq(length(info_list))){
  indel_len_plot_list[[i]] <- ggplot(data = info_list[[i]]) + 
    geom_density(aes(x = abs(SVLEN)), fill = 'grey50') + 
    scale_x_log10() + 
    facet_wrap(~SVTYPE) + 
    ggtitle(paste('branch', gsub('_', '.', names(info_list[i])), sep = '_')) + 
    xlab('log10(InDel Size)') + 
    geom_vline(xintercept = 60, linetype = 'dotted', size = 0.75) + 
    geom_vline(xintercept = 200, linetype = 'dotted', size = 0.75) + 
    geom_vline(xintercept = 1000, linetype = 'dotted', size = 0.75) + 
    geom_vline(xintercept = 4800, linetype = 'dotted', size = 0.75)
}

pdf(file = combo_len_file, height = 10, width = 8)
do.call(grid.arrange, c(indel_len_plot_list, ncol = 2))
dev.off()

quit(save = 'no')
