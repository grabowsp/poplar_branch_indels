# Test of repeatibility of the pbsv pipeline for detecting SVs in the poplar \
branches

## Overview
* Run recommended pbsv pipeline, as layed out on the pbsv Github page on \
Jan. 29 2019

## Generate Tandem Repeat Annotation .bed
### Overview
* Github page recommends including a tandem repeat annotation .bed during \
the discovery phase to "increase sensitivity and recall"
* I'm starting with trying Tandem Repeat Finder (trf)
  * I'm using a wrapper called trfbig that runs trf on large sequences and \
can generate a .bed file.
    * Was generated for the USCS human genome browser data and is no longer \
supported, so there's limited documentation on it
  * I installed both trf and trfbig (usc-trfbig) into by `NGS_analysis` \
environment via the bioconda channel
### Test run
```
cd /home/f1p1/tmp/poplar_branches/ref_stuff/tmp
trfBig /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta trf_test_1.fa

trf /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta 2 7 7 80 10 50 2000 -d

```
* Note: the trfBig wrapper is not working properly.
* Next steps: 
  * try trf on a subset of the reference genome to see how long it takes to run
  * If need be, divide the reference into chromosomes and run independently \
on each chromosome
  * Figure out how to generate a .bed file from the output of trf
