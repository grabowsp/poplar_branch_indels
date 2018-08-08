#!/bin/bash

FILE_DIR=/home/t4c1/WORK/grabowsk/data/poplar_branches/struc_vcfs

cd $FILE_DIR

for br in 13 14;
  do for i in {1..5};
    do for indel in ins del;
      do vcftools --vcf 'branch_'$br'_'$i'_chromosomes_only.recode.vcf' \
         --positions 'branch_'$br'_'$i'_'$indel'_positions.txt' --recode \
         --recode-INFO-all --out 'branch_'$br'_'$i'_chr_'$indel;
    done;
  done;
done
 

