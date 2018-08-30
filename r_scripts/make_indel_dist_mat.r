# Script for generating a distance matrix based on shared indels in the
#  poplar branch data

# LOAD FILES #
meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v3.0.txt'
meta_0 <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t')

# SET VARIABLES #
rm_samps <- c('13.4', '14.1')

file_dir <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/shared_loci/'

indel_dist_pre <- 'popbranch_'
indel_dist_suf <- '.indel_dist_tot'

#dist_cut <- 100
dist_cut_vec <- c(0, 1, 100, 1000, 5000, 10000)
# OUTPUT INFO #


# SET CONSTANTS #


# LOAD LIBRARIES #
library(ggplot2)
library(gridExtra)

#################
# remove problematic samples
meta_rm_rows <- c()
for(rms in seq(length(rm_samps))){
  tmp_ind <- which(meta_0$branch_name == rm_samps[rms])
  meta_rm_rows <- c(meta_rm_rows, tmp_ind)
} 

meta <- meta_0[-meta_rm_rows, ]

# generate scaffold for distance matrix
adj_branch_names <- paste('X', gsub('.', '_', meta$branch_name, fixed = T), 
  sep = '')

tot_dist_mat <- matrix(0, nrow = nrow(meta), ncol = nrow(meta))
rownames(tot_dist_mat) <- adj_branch_names
colnames(tot_dist_mat) <- gsub('X', 'b_', adj_branch_names)

dist_mat_list <- list()
for(dc in seq(length(dist_cut_vec))){
  dist_mat_list[[dc]] <- tot_dist_mat
#  dist_mat_list[[dc]] <- matrix(0, nrow = nrow(meta), ncol = nrow(meta))
#  rownames(dist_mat_list[[dc]]) <- adj_branch_names
#  colnames(dist_mat_list[[dc]]) <- gsub('X', 'b_', adj_branch_names)
}
names(dist_mat_list) <- paste(dist_cut_vec/1000, 'K', sep = '')

ins_mat_list <- del_mat_list <- dist_mat_list

data_rm_names <- paste('X', gsub('.', '_', rm_samps, fixed = T), sep = '')

for(i in seq(nrow(meta))){
  tmp_name <- gsub('.', '_', meta$branch_name[i], fixed = T)
  tmp_name <- paste(file_dir, indel_dist_pre, tmp_name, indel_dist_suf, 
    sep = '')
  tmp_data_0 <- read.table(tmp_name, header = T, stringsAsFactors = F, 
    sep = '\t')
  # remove problematic samples
  bad_cols <- c()
  for(k in data_rm_names){
    tmp_ind <- which(colnames(tmp_data_0) == k)
    bad_cols <- c(bad_cols, tmp_ind) 
  }
  tmp_data <- tmp_data_0[ , -bad_cols]
  #
  for(j in seq(length(dist_mat_list))){
    dist_cut <- dist_cut_vec[j]
    tmp_sim_vec <- apply(tmp_data[, c(6:ncol(tmp_data))], 2,
      function(x) sum(x <= dist_cut))
    tmp_dist_vec <- 1 - (tmp_sim_vec / nrow(tmp_data))
    dist_mat_list[[j]][names(tmp_dist_vec),i] <- tmp_dist_vec
    #
    tmp_I_sim_vec <- apply(tmp_data[which(tmp_data$INDEL_TYPE == 'INS'), 
      c(6:ncol(tmp_data))], 2, function(x) sum(x <= dist_cut))
    tmp_I_dist_vec <- 1 - (tmp_I_sim_vec / sum(tmp_data$INDEL_TYPE == 'INS'))
    ins_mat_list[[j]][names(tmp_I_dist_vec),i] <- tmp_I_dist_vec
    #
    tmp_D_sim_vec <- apply(tmp_data[which(tmp_data$INDEL_TYPE == 'DEL'), 
      c(6:ncol(tmp_data))], 2, function(x) sum(x <= dist_cut))
    tmp_D_dist_vec <- 1 - (tmp_D_sim_vec / sum(tmp_data$INDEL_TYPE == 'DEL'))
    del_mat_list[[j]][names(tmp_D_dist_vec),i] <- tmp_D_dist_vec
  }
}

make_stand_dist_mat <- function(unstand_dist_mat){
  stand_dist_mat <- matrix(0, nrow = nrow(unstand_dist_mat), 
    ncol = ncol(unstand_dist_mat))
  dimnames(stand_dist_mat) <- dimnames(unstand_dist_mat)
  for(i in seq(ncol(unstand_dist_mat))){
    stand_dist_mat[,i] <- unstand_dist_mat[,i]/max(unstand_dist_mat[,i])
  }
  return(stand_dist_mat)
}

make_sym_dist_mat <- function(unsym_dist_mat){
  sym_dist_mat <- matrix(0, nrow = nrow(unsym_dist_mat),
    ncol = ncol(unsym_dist_mat))
  dimnames(sym_dist_mat) <- dimnames(unsym_dist_mat)
  for(x in seq(ncol(unsym_dist_mat))){
    for(y in seq(nrow(unsym_dist_mat))){
      tmp_val <- mean(c(unsym_dist_mat[x,y], unsym_dist_mat[y,x]))
      sym_dist_mat[y,x] <- tmp_val
    }
  }
  return(sym_dist_mat)
}

stand_dist_list <- lapply(dist_mat_list, make_stand_dist_mat)
sym_dist_list <- lapply(stand_dist_list, make_sym_dist_mat)

plot_list <- list()
for(pl in seq(length(sym_dist_list))){
  test_mds <- data.frame(cmdscale(sym_dist_list[[pl]], k = 5),
                stringsAsFactors = F)
  colnames(test_mds) <- paste('PCo', seq(ncol(test_mds)), sep = '_')
  test_mds$tree <- paste('tree_', meta$tree, sep = '')
  test_mds$branch <- as.character(meta$branch)
  x_range <- range(test_mds$PCo_1)[2] - range(test_mds$PCo_1)[1]
  x_per_adj <- 0.2
  x_lim_adj <- c(min(test_mds$PCo_1) - (x_range * x_per_adj), 
                 max(test_mds$PCo_1) + (x_range * x_per_adj))
  y_range <- range(test_mds$PCo_2)[2] - range(test_mds$PCo_2)[1]
  y_per_adj <- 0.2
  y_lim_adj <- c(min(test_mds$PCo_2) - (y_range * y_per_adj), 
                 max(test_mds$PCo_2) + (y_range * y_per_adj))
  plot_list[[pl]] <- ggplot(data = test_mds, aes(x = PCo_1, y = PCo_2)) +
                     geom_point(aes(color = tree)) +
                     geom_text(aes(label = branch), hjust = 0, vjust = 0) +
                     ggtitle(paste('Indel PCoA with indels\nmerged w/in ', 
                       names(sym_dist_list)[pl], sep = '')) +
                     xlim(x_lim_adj) +
                     ylim(y_lim_adj)
}

indel_mds_plot_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/figs/indel_PCoA_combo.pdf'

pdf(file = indel_mds_plot_file, width = 8, height = 10)
do.call(grid.arrange, c(plot_list, ncol = 2))
dev.off()

# Plot for insertions-only
ins_stand_list <- lapply(ins_mat_list, make_stand_dist_mat)
ins_sym_list <- lapply(ins_stand_list, make_sym_dist_mat)

ins_plot_list <- list()
for(pl in seq(length(ins_sym_list))){
  test_mds <- data.frame(cmdscale(ins_sym_list[[pl]], k = 5),
                stringsAsFactors = F)
  colnames(test_mds) <- paste('PCo', seq(ncol(test_mds)), sep = '_')
  test_mds$tree <- paste('tree_', meta$tree, sep = '')
  test_mds$branch <- as.character(meta$branch)
  x_range <- range(test_mds$PCo_1)[2] - range(test_mds$PCo_1)[1]
  x_per_adj <- 0.2
  x_lim_adj <- c(min(test_mds$PCo_1) - (x_range * x_per_adj), 
                 max(test_mds$PCo_1) + (x_range * x_per_adj))
  y_range <- range(test_mds$PCo_2)[2] - range(test_mds$PCo_2)[1]
  y_per_adj <- 0.2
  y_lim_adj <- c(min(test_mds$PCo_2) - (y_range * y_per_adj), 
                 max(test_mds$PCo_2) + (y_range * y_per_adj))
  ins_plot_list[[pl]] <- ggplot(data = test_mds, aes(x = PCo_1, y = PCo_2)) +
                     geom_point(aes(color = tree)) +
                     geom_text(aes(label = branch), hjust = 0, vjust = 0) +
                     ggtitle(paste('Insertions-only PCoA\nwith indels ',
                       'merged w/in ', names(ins_sym_list)[pl], sep = '')) +
                     xlim(x_lim_adj) +
                     ylim(y_lim_adj)
}

insert_mds_plot_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/figs/insertions_PCoA_combo.pdf'

pdf(file = insert_mds_plot_file, width = 8, height = 10)
do.call(grid.arrange, c(ins_plot_list, ncol = 2))
dev.off()

# Plot for deletions-only
del_stand_list <- lapply(del_mat_list, make_stand_dist_mat)
del_sym_list <- lapply(del_stand_list, make_sym_dist_mat)

del_plot_list <- list()
for(pl in seq(length(del_sym_list))){
  test_mds <- data.frame(cmdscale(del_sym_list[[pl]], k = 5),
                stringsAsFactors = F)
  colnames(test_mds) <- paste('PCo', seq(ncol(test_mds)), sep = '_')
  test_mds$tree <- paste('tree_', meta$tree, sep = '')
  test_mds$branch <- as.character(meta$branch)
  x_range <- range(test_mds$PCo_1)[2] - range(test_mds$PCo_1)[1]
  x_per_adj <- 0.2
  x_lim_adj <- c(min(test_mds$PCo_1) - (x_range * x_per_adj),
                 max(test_mds$PCo_1) + (x_range * x_per_adj))
  y_range <- range(test_mds$PCo_2)[2] - range(test_mds$PCo_2)[1]
  y_per_adj <- 0.2
  y_lim_adj <- c(min(test_mds$PCo_2) - (y_range * y_per_adj),
                 max(test_mds$PCo_2) + (y_range * y_per_adj))
  del_plot_list[[pl]] <- ggplot(data = test_mds, aes(x = PCo_1, y = PCo_2)) +
                     geom_point(aes(color = tree)) +
                     geom_text(aes(label = branch), hjust = 0, vjust = 0) +
                     ggtitle(paste('Deletions-only PCoA\nwith indels ',
                       'merged w/in ', names(del_sym_list)[pl], sep = '')) +
                     xlim(x_lim_adj) +
                     ylim(y_lim_adj)
}

delet_mds_plot_file <- '/home/grabowsky/tools/workflows/poplar_branch_indels/figs/deletions_PCoA_combo.pdf'

pdf(file = delet_mds_plot_file, width = 8, height = 10)
do.call(grid.arrange, c(del_plot_list, ncol = 2))
dev.off()


quit(save = 'no')
