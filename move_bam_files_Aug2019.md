# Info about moving files to free up space

## Overview
* Several of the .bam files are really big and need to move them to free up \
space on the current disk
### Plan
* Move files to different disk - submit script or else will take forever
* Set up soft-links

## ngmlr BAMs
### New location for files
* `/home/grabowsky_scratch/poplar_branch_files/`
### Transfer files to new SCRATCH space
* use rsync to transfer files
* `/home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/move_bam_files_submit.sh`
### Remove bam files from old directory
```
cd /home/f1p1/tmp/poplar_branches/lib_mapping/
rm ./*.bam
```
### Generate soft-links in old directory
```
bash
cd /home/grabowsky_scratch/poplar_branch_files/
for i in *.bam;
do ln -s /home/grabowsky_scratch/poplar_branch_files/$i /home/f1p1/tmp/poplar_branches/lib_mapping/$i;
done
```

## pbmm2 BAMS
### New location for files
* `/home/grabowsky_scratch/poplar_branch_files/pbmm2_bams`
### Transfer files to new SCRATCH space
* `/home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts/move_pbmm2_bam_files_submit.sh`
### Remove bam files from old directory
```
cd /home/f1p2/tmp/Poplar14.5_pbsv/pbsv_v2.1.1_run
rm ./*.bam
```
### Generate soft-links in old directory
```
bash
cd /home/grabowsky_scratch/poplar_branch_files/pbmm2_bams
for i in *.bam;
do ln -s /home/grabowsky_scratch/poplar_branch_files/pbmm2_bams/$i /home/f1p2/tmp/Poplar14.5_pbsv/pbsv_v2.1.1_run/$i;
done
```


