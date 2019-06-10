#!/bin/bash

#$ -q "all.q"
#$ -pe smp 1
#$ -M pgrabowski@hudsonalpha.org
#$ -m ae
#$ -V
#$ -cwd
#$ -N stats_8lib_r02

/home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/generate_pbsv_SV_stats.sh /home/grabowsky/tools/workflows/poplar_branch_indels/pbsv_res_info/pbsv_v2.2.0_ngmlr_8libs_r02.res_config


