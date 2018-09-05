#!/bin/bash
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams
samtools sort /home/smrtlink/userdata/jobs_root.local/001/001961/tasks/pbsvtools.tasks.gather_align-1/alignments.bam -o branch13.1_sorted.bam
samtools sort /home/smrtlink/userdata/jobs_root.local/001/001965/tasks/pbsvtools.tasks.gather_align-1/alignments.bam -o branch13.2_sorted.bam
samtools sort /home/smrtlink/userdata/jobs_root.local/001/001969/tasks/pbsvtools.tasks.gather_align-1/alignments.bam -o branch13.3_sorted.bam
samtools sort /home/smrtlink/userdata/jobs_root.local/002/002040/tasks/pbsvtools.tasks.gather_align-1/alignments.bam -o branch13.5_sorted.bam
samtools sort /home/smrtlink/userdata/jobs_root.local/002/002243/tasks/pbsvtools.tasks.gather_align-1/alignments.bam -o branch14.2_sorted.bam
samtools sort /home/smrtlink/userdata/jobs_root.local/002/002269/tasks/pbsvtools.tasks.gather_align-1/alignments.bam -o branch14.3_sorted.bam
samtools sort /home/smrtlink/userdata/jobs_root.local/002/002245/tasks/pbsvtools.tasks.gather_align-1/alignments.bam -o branch14.4_sorted.bam
samtools sort /home/smrtlink/userdata/jobs_root.local/002/002060/tasks/pbsvtools.tasks.gather_align-1/alignments.bam -o branch14.5_sorted.bam
