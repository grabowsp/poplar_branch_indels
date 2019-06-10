#!/bin/bash

#$ -q "all.q"
#$ -pe smp 1
#$ -M pgrabowski@hudsonalpha.org
#$ -m ae
#$ -V
#$ -cwd
#$ -N ngmlr_r03_perm_stats


source /home/raid2/LINUXOPT/miniconda2/bin/activate /home/grabowsky/.conda/envs/NGS_analysis

Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_Indel_permissive_analysis.r /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.ngmlr.ppsv_v2.2_1.full.call.r03.vcf_data_info.rds
