#!/bin/bash

#$ -q "all.q"
#$ -pe smp 1
#$ -M pgrabowski@hudsonalpha.org
#$ -m ae
#$ -V
#$ -cwd
#$ -N pbmm2_r02_stats

/home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/generate_pbsv_SV_stats.sh /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PtStettler14.pbmm2.ppsv_v2.2_1.full.call.r02.vcf_data_info.rds


