# Compare individual VCF files from runs through the pipeline

data_dir <- '/home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL/'

old_ver_short <- 'ref.PAXL.vcf'
new_ver_short <- 'Ptr145v1.PAXL_v2_r1.vcf'

tmp_old_vcf <- read.table(paste(data_dir, old_ver_short, sep = ''), sep = '\t',
  stringsAsFactors = F)

tmp_new_vcf <- read.table(paste(data_dir, new_ver_short, sep = ''), sep = '\t',
  stringsAsFactors = F)

old_info <- strsplit(tmp_old_vcf[ , 8], split = ';')
old_type <- unlist(lapply(old_info, 
  function(x) unlist(strsplit(x[[1]], split = '='))[2]))

tmp_old_vcf$type <- old_type
tmp_old_vcf$name <- paste(tmp_old_vcf[,1], tmp_old_vcf[,2], sep = '_')

new_info <- strsplit(tmp_new_vcf[ , 8], split = ';')
new_type <- unlist(lapply(new_info, 
  function(x) unlist(strsplit(x[[1]], split = '='))[2]))

tmp_new_vcf$type <- new_type
tmp_new_vcf$name <- paste(tmp_new_vcf[,1], tmp_new_vcf[,2], sep = '_')

# NUMBER OF EACH SV
table(tmp_old_vcf$type)
#  BND   DEL   INS   INV 
# 1066 26439 24022     7
table(tmp_new_vcf$type)
#   BND   DEL   INS   INV 
# 10364 25263 25358    21

# Almost 10x (10k) more BND than in old VCF; similar levels of other SVs

sum(tmp_new_vcf$name %in% tmp_old_vcf$name)
# [1] 33152
# 33k SVs at same position
new_in_old <- which(tmp_new_vcf$name %in% tmp_old_vcf$name)
new_vcf_old_pos <- sapply(tmp_new_vcf$name[new_in_old], function(x) 
  which(tmp_old_vcf$name == x))

new_in_old_nhits <- unlist(lapply(new_vcf_old_pos, length))
table(new_in_old_nhits)
#     1     2 
# 33028   124 

new_vcf_old_pos_1 <- unlist(lapply(new_vcf_old_pos, function(x) x[1]))
sum(tmp_new_vcf$type[new_in_old] != tmp_old_vcf$type[new_vcf_old_pos_1])
# 105

not_same_type_pos1 <- which(tmp_new_vcf$type[new_in_old] != 
  tmp_old_vcf$type[new_vcf_old_pos_1])

# tmp_new_vcf[new_same_pos[not_same_type[1]], ]
# tmp_old_vcf[new_vcf_old_pos_1[not_same_type][1], ]

length(intersect(which(new_in_old_nhits == 2), not_same_type_pos1))
# 30

not_same_2hits <- intersect(which(new_in_old_nhits == 2), not_same_type_pos1)

tmp_new_vcf$type[new_same_pos[not_same_2hits[1]]]
tmp_old_vcf$type[new_vcf_old_pos[[not_same_2hits[1]]][2]]

tmp_new_vcf$type[new_same_pos[not_same_2hits]]

not_same_2hit_list <- new_vcf_old_pos[not_same_2hits]
not_same_2hit_2ind <- unlist(lapply(not_same_2hit_list, function(x) x[2]))

sum(tmp_new_vcf$type[new_in_old[not_same_2hits]] != 
  tmp_old_vcf$type[not_same_2hit_2ind])
# 0
# the 30 'not_same_2hits' are all the same SV type when looking at the second
#   index

new_in_old_diff_type <- new_in_old[setdiff(not_same_type_pos1, not_same_2hits)]

length(new_in_old_diff_type)
# [1] 75
# so 105-30 = 75 SVs with same position are not the same type - need to 
#   look at those

# 1st is a INS in old and BND in new
# 2nd is BND in old and INS in new

diff_type_new_type <- tmp_new_vcf$type[new_in_old_diff_type]
diff_type_old_type <- tmp_old_vcf$type[unlist(new_vcf_old_pos[
  setdiff(not_same_type_pos1, not_same_2hits)])]
diff_type_df <- data.frame(old_type = diff_type_old_type, 
  new_type = diff_type_new_type, stringsAsFactors = F)

diff_type_bnd_inds <- union(which(diff_type_df$old_type == 'BND'), 
  which(diff_type_df$new_type == 'BND'))

length(diff_type_bnd_inds)
# [1] 45

diff_type_indel_inds <- setdiff(c(1:nrow(diff_type_df)), diff_type_bnd_inds)
# 30

tmp_new_vcf[new_in_old_diff_type[diff_type_indel_inds[3]], ]
tmp_old_vcf[which(tmp_old_vcf$name == tmp_new_vcf$name[
  new_in_old_diff_type[diff_type_indel_inds[3]]]),]

tmp_new_vcf[new_in_old_diff_type[diff_type_indel_inds], 8]
# 27 of 30 are tagged as TANDEM duplicates in tmp_new_vcf - probably would be
#   filtered out of final results

# Summary of comparisons between v2.0.1 and v2.1.1
# 33152 of 61016 (in new vcf) SVs are at the same position in both files
# of those, all but 75 are the same type of SV
# of those 75, 45 are different because the SV is a BND in one file but
#  not the other
# of the remaining 30, 27 are labeled as TANDEM repeats and may be removed
#  following filtering


old_same_pos <- which(tmp_old_vcf$name %in% tmp_new_vcf$name)

sum(tmp_old_vcf$type[old_same_pos] == tmp_new_vcf$type[new_same_pos])

old_bnd_inds <- which(old_type == 'BND')

new_bnd_inds <- which(new_type == 'BND')
# BND data is sort of a disaster
# The new round of results essentially has 10k more BND's, which seems odd...
# Ways to filter:
#   Want 2 adjacent BNDs within a certain distance 
      # otherwise, i'm not sure what's going on
#   Want adjacent BNDs to have correct orientation: first BND have 
      # the new sequence on the 5' end, second BND have new sequence on 3' end
#   Adjacent BNDs should typically have mates from the same chromosome
      # I'm not sure that represents true translocations, but my suspision
      #   is that many of the not-garbage BNDs are actually insertions of
      #   repetitive sequences

chr1_bnds <- intersect(which(tmp_new_vcf[,1] == 'Chr01'), new_bnd_inds)

bnd_3prime <- intersect(grep('^A|^T|^C|^G', tmp_new_vcf[,5]), new_bnd_inds)
bnd_5prime <- setdiff(new_bnd_inds, bnd_3prime)

bnd_3prime_cor_orient <- bnd_3prime[which((bnd_3prime + 1) %in% bnd_5prime)]

adj_dist

bnd_3prime_dist <- -1 * (tmp_new_vcf[bnd_3prime, 2] - 
  tmp_new_vcf[(bnd_3prime + 1),2])
bnd_3prime_dist_cut <- bnd_3prime[which(bnd_3prime_dist < 50 & 
  bnd_3prime_dist >= 0)]

bnd_3prime_cor_so <- intersect(bnd_3prime_dist_cut, bnd_3prime_cor_orient)

bnd_3prime_ginfo <- strsplit(tmp_new_vcf[bnd_3prime,10], split = ':')
bnd_3prime_ac <- lapply(bnd_3prime_ginfo, function(x) strsplit(x[[2]], 
  split = ','))
bnd_3prime_ac_mat <- matrix(data = as.numeric(unlist(bnd_3prime_ac)), 
  ncol = 2, byrow = T)

bnd_

bnd_3prime_min_count <- apply(bnd_3prime_ac_mat, 1, min)
bnd_3prime_max_count <- apply(bnd_3prime_ac_mat, 1, max)
bnd_3prime_depth <- as.numeric(unlist(lapply(bnd_3prime_ginfo, 
  function(x) x[[3]])))
bnd_3prime_min_pen <- bnd_3prime_min_count/bnd_3prime_depth

bnd_3prime_count_pass <- bnd_3prime[which(bnd_3prime_min_count > 2 & 
  bnd_3prime_min_pen > 0.1)]

bnd_3prime_cor_so_count <- intersect(bnd_3prime_cor_so, bnd_3prime_count_pass)

new_inv_inds <- which(new_type == 'INV')
