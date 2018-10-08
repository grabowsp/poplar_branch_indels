#!/bin/bash
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/multisamp_indel_calling
samtools mpileup -f /home/smrtlink/references/Ptrichocarpa_444_v3_0/sequence/Ptrichocarpa_444_v3_0.fasta -o pop_branch.mpileup /home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams/branch13.1_sorted.bam /home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams/branch13.2_sorted.bam /home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams/branch13.3_sorted.bam /home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams/branch13.5_sorted.bam /home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams/branch14.2_sorted.bam /home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams/branch14.3_sorted.bam /home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams/branch14.4_sorted.bam /home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams/branch14.5_sorted.bam