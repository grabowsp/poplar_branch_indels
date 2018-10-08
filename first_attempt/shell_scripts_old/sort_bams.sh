#!/bin/bash
cd /home/t4c1/WORK/grabowsk/data/poplar_branches/sorted_bams
/home/raid2/LINUXOPT/samtools-1.9/bin/samtools sort -T branch13.1tmpsort -m 2G -o branch13.1_sorted.bam /home/smrtlink/userdata/jobs_root.local/001/001961/tasks/pbsvtools.tasks.gather_align-1/alignments.bam
/home/raid2/LINUXOPT/samtools-1.9/bin/samtools sort -T branch13.2tmpsort -m 2G -o branch13.2_sorted.bam /home/smrtlink/userdata/jobs_root.local/001/001965/tasks/pbsvtools.tasks.gather_align-1/alignments.bam
/home/raid2/LINUXOPT/samtools-1.9/bin/samtools sort -T branch13.3tmpsort -m 2G -o branch13.3_sorted.bam /home/smrtlink/userdata/jobs_root.local/001/001969/tasks/pbsvtools.tasks.gather_align-1/alignments.bam
/home/raid2/LINUXOPT/samtools-1.9/bin/samtools sort -T branch13.5tmpsort -m 2G -o branch13.5_sorted.bam /home/smrtlink/userdata/jobs_root.local/002/002040/tasks/pbsvtools.tasks.gather_align-1/alignments.bam
/home/raid2/LINUXOPT/samtools-1.9/bin/samtools sort -T branch14.2tmpsort -m 2G -o branch14.2_sorted.bam /home/smrtlink/userdata/jobs_root.local/002/002243/tasks/pbsvtools.tasks.gather_align-1/alignments.bam
/home/raid2/LINUXOPT/samtools-1.9/bin/samtools sort -T branch14.3tmpsort -m 2G -o branch14.3_sorted.bam /home/smrtlink/userdata/jobs_root.local/002/002269/tasks/pbsvtools.tasks.gather_align-1/alignments.bam
/home/raid2/LINUXOPT/samtools-1.9/bin/samtools sort -T branch14.4tmpsort -m 2G -o branch14.4_sorted.bam /home/smrtlink/userdata/jobs_root.local/002/002245/tasks/pbsvtools.tasks.gather_align-1/alignments.bam
/home/raid2/LINUXOPT/samtools-1.9/bin/samtools sort -T branch14.5tmpsort -m 2G -o branch14.5_sorted.bam /home/smrtlink/userdata/jobs_root.local/002/002060/tasks/pbsvtools.tasks.gather_align-1/alignments.bam
