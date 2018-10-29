#!/bin/bash

FILE_DIR=/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller

cd $FILE_DIR

for branch in `ls *.vcf`;
    do vcftools --vcf $branch --chr Chr01 \
      --chr Chr02 --chr Chr03 --chr Chr04 --chr Chr05 --chr Chr06 --chr Chr07\
      --chr Chr08 --chr Chr09 --chr Chr10 --chr Chr11 --chr Chr12 --chr Chr13\
      --chr Chr14 --chr Chr15 --chr Chr16 --chr Chr17 --chr Chr18 --chr Chr19\
      --recode --recode-INFO-all \
      --out $branch'_chromOnly.vcf';
done

mv *chromOnly.vcf.recode.vcf ./chrom_only


