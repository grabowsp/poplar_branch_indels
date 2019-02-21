# Test of repeatibility of the pbsv pipeline for detecting SVs in the poplar \
branches

## Overview
* Run recommended pbsv pipeline, as layed out on the pbsv Github page on \
Jan. 29 2019
* Run 3 times on 3 samples to look at repeatibility
  * PAXL
  * PAXN
  * PAYK

## Generate Tandem Repeat Annotation .bed
### Overview
* Github page recommends including a tandem repeat annotation .bed during \
the discovery phase to "increase sensitivity and recall"
### Info about generating .bed file
* `./generate_tandem_repeat_bed.md`
### .bed file
* `/home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes/poplar_14_5_v1_tandemrepeat.bed`

## Merge .subreads.bam files to make single file for mapping
### Overview
* I tried to incorporate separate `...subreads.bam` files into a single \
`pbmm2 align` command but couldn't get it to work
* It's not clear how to properly map the `...subreads.bam` files individually \
and then merge the outputs
* So will merge the `...subreads.bam` files into a single `.bam` file for \
mapping
* will use bamtools `merge` function for that
  * requires making a file of file names of the `...subreads.bam` files
### PAXL
* Just testing that `bamtools merge` will work
  * Take-home: it workds, but it takes a long time to run, so should submit \
the job on the cluster - see PAXN
```
bash
source activate NGS_analysis
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
ls m*subreads.bam > PAXL_subreads_files.fofn
bamtools merge -list PAXL_subreads_files.fofn -out PAXL.all.subreads.bam

```
### PAXN
* testing that running it via submision works properly
  * Take-home: it works
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXN
ls m*subreads.bam > PAXN_subreads_files.fofn
qsub PAXN.mergesubreads.sh
```
### Generate file of subread file names for each library
* already manually generated them for PAXL and PAXN
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do cd ./$LIB
ls m*subreads.bam > $LIB'_subreads_files.fofn'
cd ..; done
```
### Generate submit files
* already merged the files for PAXL and PAXN
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do sed 's/PAXN/'"$LIB"'/g' ./PAXN/PAXN.mergesubreads.sh \
> ./$LIB/$LIB.mergesubreads.sh; done
```
### Submit merge jobs
* already merged the files for PAXL and PAXN
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do cd ./$LIB
qsub $LIB.mergesubreads.sh
cd ..; done
```

## Align PacBio reads to refernce genome
### Overview
* use `pbmm2` to map `..subreads.bam` files to reference
  * had to merge the `subreads.bam` files to get this to work properly
### Example code
```
source /home/raid2/LINUXOPT/miniconda2/bin/activate /home/grabowsky/.conda/envs/NGS_analysis

pbmm2 align \
/home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta \
PAXL.all.subreads.bam ref.PAXL_v2_r1.bam --sort --sample 'PAXL' --median-filter
```
### PAXL
* Check that command words
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.align_v2_r1.sh
```
* Make submit scripts for reps 2 and 3
```
for i in _r2 _r3;
do 
sed 's/_r1/'"$i"'/g' PAXL.align_v2_r1.sh > PAXL.align_v2$i.sh;
done
```
* submit jobs
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.align_v2_r2.sh
qsub PAXL.align_v2_r3.sh
```
### PAXN and PAYK
* Generate submit scripts
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAXN PAYK;
  do sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.align_v2_r1.sh > \
  ./$LIB/$LIB.align_v2_r1.sh
  for i in _r2 _r3;
    do sed 's/_r1/'"$i"'/g' ./$LIB/$LIB.align_v2_r1.sh > \
    ./$LIB/$LIB.align_v2$i.sh;
  done;
done
```
* Submit jobs
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAXN PAYK;
  do cd ./$LIB
  for i in _r1 _r2 _r3;
#    do echo $LIB.align_v2$i.sh;
    do qsub $LIB.align_v2$i.sh;
  done;
  cd ..;
done
```
### Rest of the Libraries
* Generate submit scripts
  * Already generated scripts for PAXL, PAXN, and PAYK
  * Only doing one rep (for now) for remaining libraries
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
  do sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.align_v2_r1.sh > \
  ./$LIB/$LIB.align_v2_r1.sh;
done
```
* Submit scripts
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
  do cd ./$LIB
  qsub $LIB.align_v2_r1.sh
  cd ..;
done
```

## Discover signatures of structural variation
### Overview
* use `pbsv discover` to generate `...svsig.gz` file
* GitHub page highly recommends using a tandem repeat annotation .bed file
### Example code
```
pbsv discover --tandem-repeats \
/home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes/poplar_14_5_v1_tandemrepeat.bed \
ref.PAXL_v2_r1.bam Ptr145v1.PAXL_v2_r1.svsig.gz
```
### PAXL
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.disc_v2_r1.sh

```

## Call structural vairants and assign genotypes
### Overview
* use `pbsv call` to generate `...vcf` file
### Example Code

### PAXL
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.call_v2_r1.sh

```

* NEXT Steps:
  * check on all the mapping results
  * Check on PAXL SV calling/genotyping results
  * Generate commands for discovery for rest of the libraries
  * Submit discovery jobs
  * Check on discovery results
  * Generate commands for sv calling for individuals
  * Submit sv calling for individuals jobs
  * Generate commands fof sv calling using all samples
    * Generate 4 with different sample orders
  * Submit sv calling for all samples jobs
