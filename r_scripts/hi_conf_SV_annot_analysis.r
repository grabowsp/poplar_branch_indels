# Analysis of the hi-confidence SVs that might be somatic mutations

annot_file <- '/home/t4c1/WORK/sujan/Ptricocarpa/PtrichocarpaStettler14_532_v1.1/annotation/PtrichocarpaStettler14_532_v1.1.gene.gff3'

annot <- read.table(annot_file, header = F, stringsAsFactors = F)

annot_gene_inds <- which(annot[,3] == 'gene')

####
gene_in_SV <- function(annot_df, chr_name, sv_start, sv_length){
  ########
  chr_inds <- which(annot_df[, 1] == chr_name)
  st_inds <- intersect(which(annot_df[, 4] >= sv_start), chr_inds)
  sv_end <- sv_start + sv_length
  en_inds <- intersect(which(annot_df[, 5] <= sv_end), chr_inds)
  overlap_inds <- intersect(st_inds, en_inds)
  return(overlap_inds)
}

gene_5overlap_SV <- function(annot_df, chr_name, sv_start, sv_length){
  ########
  chr_inds <- which(annot_df[, 1] == chr_name)
  st_5_inds <- intersect(which(annot_df[,4] < sv_start), chr_inds)
  en_5_inds <- intersect(which(annot_df[,5] >= sv_start), chr_inds)
  overlap_inds <- intersect(st_5_inds, en_5_inds)
  return(overlap_inds)
}

gene_3overlap_SV <- function(annot_df, chr_name, sv_start, sv_length){
  ######
  chr_inds <- which(annot_df[, 1] == chr_name)
  sv_end <- sv_start + sv_length
  st_3_inds <- intersect(which(annot_df[,4] <= sv_end), chr_inds)
  en_3_inds <- intersect(which(annot_df[,5] > sv_end), chr_inds)
  overlap_inds <- intersect(st_3_inds, en_3_inds)
  return(overlap_inds)
}

################
# Chr01_38823408_INS_67
ins_chr1_5_inds <- gene_5overlap_SV(annot_df = annot, chr_name = 'Chr01',
  sv_start = 38823408, sv_length = 0)
# integer(0)

ins_chr1_3_inds <- gene_3overlap_SV(annot_df = annot, chr_name = 'Chr01',
  sv_start = 38823408, sv_length = 0)
# integer(0)

chr01_inds <- which(annot[,1] == 'Chr01')

st_inds_1 <- intersect(which(annot[,4] <= 38823408), chr01_inds)
ed_inds_1 <- intersect(which(annot[,5] >= (38823408)), chr01_inds)

intersect(st_inds_1, ed_inds_1)
# integer(0)

max(st_inds_1)
# closest to PtStettler14.01G320000.v1.1

min(ed_inds_1)
# closest to PtStettler14.01G319900.1.v1.1

# seems to be between PtStettler14.01G319900.1.v1.1 and 
##  PtStettler14.01G320000.v1.1

# annot[grep('PtStettler14.01G319900', annot[,9]),]
# gene spas from 38817250 to 38823387
# insertion is 21bp downstream from end of the gene

# annot[grep('PtStettler14.01G320000', annot[,9]),]
# gene spans from 38823442 to 38827166 = 3724bp
# insertion is 34bp upstream of the end of the gene (gene is in 'negative'
#   orientiation

# Takehome: Insertion is in intergenic space between PtStettler14.01G319900
#  and PtStettler14.01G320000

################
# Chr05_17223380_DEL_9729
del_chr5_5_inds <- gene_5overlap_SV(annot_df = annot, chr_name = 'Chr05',
  sv_start = 17223380, sv_length = 9729)
intersect(del_chr5_5_inds, annot_gene_inds)
# [1] 179905
# PtStettler14.05G146600
annot[grep('PtStettler14.05G146600', annot[,9]),]
# spans from 17219228 to 17224228; 5000bp
# 848bp are within the deletion; these are all with the 1687bp 3' UTR of the
#   gene; so deletion is within the 3' UTR of the gene

del_chr5_IN_inds <- gene_in_SV(annot_df = annot, chr_name = 'Chr05',
  sv_start = 17223380, sv_length = 9729)
intersect(del_chr5_IN_inds, annot_gene_inds)
# integer(0)

del_chr5_3_inds <- gene_3overlap_SV(annot_df = annot, chr_name = 'Chr05',
  sv_start = 17223380, sv_length = 9729)
# integer(0)

chr05_inds <- which(annot[,1] == 'Chr05')

st_inds_2 <- intersect(which(annot[,4] <= 17223380), chr05_inds)
ed_inds_2 <- intersect(which(annot[,5] >= (17223380 + 9729)), chr05_inds)

intersect(st_inds_2, ed_inds_2)
# integer(0)

max(st_inds_2)
#[1] 179911
# PtStettler14.05G146600.1.v1.1

min(ed_inds_2)
#[1] 179912
# PtStettler14.05G146700

annot[grep('PtStettler14.05G146600', annot[,9]),]
# ends at 17224228; DEL start 848bp downstream of 3' UTR

annot[grep('PtStettler14.05G146700', annot[,9]),]
# starts at 17236555, DEL ends 3446bp upstream of 5' UTR

# Take-home: The deletion includes a portion of the 3' UTR of 
#  PtStettler14.05G146600 and the intergenic region between that gene and
#  PtStettler14.05G146700

###############
# Chr06_25962893_DUP_28528
dup_chr6_5_inds <- gene_5overlap_SV(annot_df = annot, chr_name = 'Chr06', 
  sv_start = 25962893, sv_length = 28528)
intersect(dup_chr6_5_inds, annot_gene_inds)
# [1] 226698
# PtStettler14.06G236700 overlaps the 5' end of the Duplication
annot[grep('PtStettler14.06G236700', annot[,9]), ]
# spans 25962085 to 25967594; 5509bp long gene
# 4701bp is within the DUP; this is CDS.2 and beyond for both mRNA versions;
### but is NOT a full version of the gene

dup_chr6_IN_inds <- gene_in_SV(annot_df = annot, chr_name = 'Chr06', 
  sv_start = 25962893, sv_length = 28528)
intersect(dup_chr6_IN_inds, annot_gene_inds)
# [1] 226716
# PtStettler14.06G236800 is completely contained within the duplication

dup_chr6_3_inds <- gene_3overlap_SV(annot_df = annot, chr_name = 'Chr06',
  sv_start = 25962893, sv_length = 28528)
# integer(0)


####
chr06_inds <- which(annot[,1] == 'Chr06')

st_inds_3 <- intersect(which(annot[,5] >= 25962893), chr06_inds)
ed_inds_3 <- intersect(which(annot[,4] <= (25962893 + 28528)), chr06_inds)

intersect(st_inds_3, ed_inds_3)
# 23 entries

max(st_inds_3)
# [1] 226710

min(ed_inds_3)
# [1] 226725

# Take-home: Duplication includes the 4701 3' region of PtStettler14.06G236700
## and all of PtStettler14.06G236800

################
# Chr09_9139125_DEL_7970
del_chr9_5_inds <- gene_5overlap_SV(annot_df = annot, chr_name = 'Chr09',
  sv_start = 9139125, sv_length = 7970)
# integer(0)

del_chr9_IN_inds <- gene_in_SV(annot_df = annot, chr_name = 'Chr09',
  sv_start = 9139125, sv_length = 7970)
# integer(0)

del_chr9_3_inds <- gene_3overlap_SV(annot_df = annot, chr_name = 'Chr09',
  sv_start = 9139125, sv_length = 7970)
# integer(0)

chr09_inds <- which(annot[,1] == 'Chr09')

st_inds_tmp <- intersect(which(annot[,4] <= 9139125), chr09_inds)
ed_inds_tmp <- intersect(which(annot[,5] >= (9139125 + 7970)), chr09_inds)

max(st_inds_tmp)
#[1] 293134
# PtStettler14.09G091900
annot[grep('PtStettler14.09G091900', annot[,9]), ]
# 1537 upstream of beginning of deletion

min(ed_inds_tmp)
#[1] 293135 
# PtStettler14.09G092000
# 1358bp downsteam of end of deletion

# Take-home: Deletion is in intergenic region between PtStettler14.09G091900 
#  and PtStettler14.09G092000

############
# Chr14_8492556_DUP_7858
dup_chr14_5_inds <- gene_5overlap_SV(annot_df = annot, chr_name = 'Chr14',
  sv_start = 8492556, sv_length = 7858)
intersect(dup_chr14_5_inds, annot_gene_inds)
# [1] 409795
# PtStettler14.14G106000
annot[grep('PtStettler14.14G106000', annot[,9]), ]
# Most of PtStettler14.14G106000 is within the duplication except for 
#  a portion of the 3' UTR of one of the transcripts and a short coding
#  sequencen

dup_chr14_IN_inds <- gene_in_SV(annot_df = annot, chr_name = 'Chr14',
  sv_start = 8492556, sv_length = 7858)
intersect(dup_chr14_IN_inds, annot_gene_inds)
# all the indices are from PtStettler14.14G106000

dup_chr14_3_inds <- gene_3overlap_SV(annot_df = annot, chr_name = 'Chr14',
  sv_start = 8492556, sv_length = 7858)
# integer(0)

# Take-home: Duplication includes most of PtStettler14.14G106000 except for
#  the 3' UTR and a putative coding region 3' of the 3' UTR.



