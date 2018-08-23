#!/bin/bash

FILE_DIR=/home/t4c1/WORK/grabowsk/data/poplar_branches/struc_vcfs

cd $FILE_DIR

for branch in 13 14;
  do for i in {1..5};
    do vcftools --vcf 'branch_'$branch'_'$i'_struct.vcf' --chr Chr01 \
      --chr Chr02 --chr Chr03 --chr Chr04 --chr Chr05 --chr Chr06 --chr Chr07\
      --chr Chr08 --chr Chr09 --chr Chr10 --chr Chr11 --chr Chr12 --chr Chr13\
      --chr Chr14 --chr Chr15 --chr Chr16 --chr Chr17 --chr Chr18 --chr Chr19\
      --recode --recode-INFO-all \
      --out 'branch_'$branch'_'$i'_chromosomes_only';
  done;
done
