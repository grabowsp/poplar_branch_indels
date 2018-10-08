#!/bin/bash

FILE_DIR=/home/t4c1/WORK/grabowsk/data/poplar_branches/struc_vcfs

cd $FILE_DIR

for i in {1..5};
  do for j in 13 14;
    do vcftools --vcf 'branch_'$j'_'$i'_chromosomes_only.recode.vcf' --depth\
          --out 'branch_'$j'_'$i;
    done;
  done


