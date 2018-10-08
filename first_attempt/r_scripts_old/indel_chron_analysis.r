# Look at chronology of indels - patterns of shared indels in the chronology
#  of the branches

# Goals
# - Quantify shared indels as go from oldest-to-youngest branches
# - Compare to random-order shared indels
# - identify genotype-specific indels (different from reference)
# - identify clone-specific indels

# LOAD FILES #
test_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/shared_loci/popbranch_13_5.indel_dist_tot'
test_data <- read.table(test_file, header = T, sep = '\t', stringsAsFactors = F)

# SET VARIABLES #


# SET CONSTANTS #


# SET OUTPUT INFO #


# LOAD PACKAGES #


##########
sum(apply(test_data[,c(6:14)], 1, sum) == 0)
# [1] 523
# 523 indels are shared between 13.5 and all other libraries

summary(test_data$INDEL_SIZE[which(apply(test_data[,c(6:14)], 1, sum) == 0)])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   51.0    56.5    81.0   237.9   209.5 11241.0 

table(test_data$INDEL_TYPE[which(apply(test_data[,c(6:14)], 1, sum) == 0)])
# DEL INS 
# 486  37 
# The vast majority of shared deletions at same spot are deletions

test_mat_1 <- test_data[, c(6:14)] > 20
table(test_data$INDEL_TYPE[which(apply(test_mat_1, 1, sum) == 0)])
#   DEL   INS 
# 17554 10010

test_mat_2 <- test_data[, c(6:14)] > 100
table(test_data$INDEL_TYPE[which(apply(test_mat_2, 1, sum) == 0)])
#   DEL   INS 
# 17831 15598
# If we start being more permissive with the distance between the positions,
#  then the difference between deletions and insertions goes down dramatically

test_mat_3 <- test_data[, c(6:14)] > 1
table(test_data$INDEL_TYPE[which(apply(test_mat_3, 1, sum) == 0)])
#   DEL   INS 
# 17376   175 
# If we chose 'same' and 'overlap', then deletions are a MUCH better choice
sum(test_data$INDEL_TYPE == 'DEL')
#[1] 21763
# ~80% of deletions are the same or have overlap across all samples

### What about for big indels
indel_size <- 1000
sub_data <- test_data[which(test_data$INDEL_SIZE >= indel_size), ]

sum(apply(sub_data[,c(6:14)], 1, sum) == 0)
# [1] 17

sub_mat_1 <- sub_data[, c(6:14)] > 1
sum(apply(sub_mat_1, 1, sum) == 0)
# [1] 2685
table(sub_data$INDEL_TYPE[which(apply(sub_mat_1, 1, sum) == 0)])
#  DEL  INS 
# 2657   28
sum(sub_data$INDEL_TYPE == 'DEL')
# [1] 3622
# ~73% of deletions 1k are larger are shared/overlap across all samples

# FOR NOW, WILL FOCUS ON DELETIONS 1K OR LARGER
indel_size <- 1000

bad_samps <- c('X13_4', 'X14_1')
bad_inds <- c()
for(bs in bad_samps){
  tmp_ind <- which(colnames(test_data) == bs)
  bad_inds <- c(bad_inds, tmp_ind)
}
keep_cols <- setdiff(seq(ncol(test_data)), bad_inds)

test_2 <- test_data[intersect(which(test_data$INDEL_TYPE == 'DEL'), 
  which(test_data$INDEL_SIZE >= indel_size)), keep_cols]

b13_5v3 <- which(test_2$X13_3 <= 1)
b13_5v2 <- which(test_2$X13_2 <= 1)
b13_5v1 <- which(test_2$X13_1 <= 1)

length(intersect(b13_5v3, b13_5v2))
# [1] 3253

length(intersect(b13_5v3, b13_5v1))
# [1] 3404

length(intersect(b13_5v2, b13_5v1))
# [1] 3261

length(intersect(intersect(b13_5v3, b13_5v2), b13_5v1))
# [1] 3220
# 3220 1K+ deletions are shared by the 4 clone-13 branches
3220/nrow(test_2)
# [1] 0.8890116
# 89% of 1K+ deletions are shared by the 4 clone-13 branches

#subt_1 <- test_2[, c(6:ncol(test_2))] > 1
#sum(apply(subt_1[,c(1:3)], 1, sum) == 0)
##[1] 3220

subt_1 <- test_2[, c(6:ncol(test_2))] > 0
sum(apply(subt_1[,c(1:3)], 1, sum) == 0)
# [1] 217

col_vec <- c(1:ncol(subt_1))
share_sum_vec <- c()
for(i in seq(1000)){
  tmp_val <- sum(apply(subt_1[,sample(col_vec, size = 3)], 1, sum) == 0)
  share_sum_vec <- c(share_sum_vec, tmp_val)
}

summary(share_sum_vec)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     115     147     162     164     182     222 

sum(share_sum_vec > 217)
# [1] 25

sum(apply(subt_1, 1, sum) == 0)
# 32
# 32 deletions are shared across all libraries

all_lib_inds <- which(apply(subt_1, 1, sum) == 0)
subt_2 <- subt_1[-all_lib_inds,]

sum(apply(subt_2[,c(1:3)], 1, sum) == 0)
# [1] 185

share_sum_vec <- c()
for(i in seq(1000)){
  tmp_val <- sum(apply(subt_2[,sample(col_vec, size = 3)], 1, sum) == 0)
  share_sum_vec <- c(share_sum_vec, tmp_val)
}

sum(share_sum_vec > 185)
# [1] 22

apply(subt_2, 2, sum)
# X13_1 X13_2 X13_3 X14_2 X14_3 X14_4 X14_5 
#  3022  3060  2959  3037  3071  3108  3087

meta_file <- '/home/t4c1/WORK/grabowsk/data/poplar_branches/meta/poplar_branch_meta_v3.0.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t')

# Deletions only found in branch 5:
sum(apply(subt_2,1,sum) == 7)
# [1] 2201

# Deletions NOT found in 14, so only found in clone 13
length(setdiff(which(apply(subt_2[, c(4:7)], 1, sum) == 4), 
  which(apply(subt_2,1,sum) == 7)))
# [1] 260

sum(apply(subt_2[, c(1,3,4)], 1, sum) == 0)
# 1,3,4

sum(apply(subt_2[, c(4:7)], 1, sum) == 4)

# Need to figure out how to incorportate/remove genotype-specific markers...
quit(save = 'no')

