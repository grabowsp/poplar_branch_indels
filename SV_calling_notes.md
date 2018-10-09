# Detecting and calling SV in poplar data using poplar ref v4.0
## Overview
### Goals
* Call SVs in branches via multiple methods
* Identify method(s) that call SVs with enough confidence to map accumulation
of somatic SV mutations going up the tree
* Validate some SVs using Illumina data
* Describe "what are the SV mutations?" (general goal of project)
* Characterize where are the mutations found and do they affect genes? (general
goal of project)
### Steps
1) Call SVs using HAGSC PacBio pipeline - each sample separate (default)
2) Call SVs using HAGSC PacBio pipeline - all samples together
3) Call SVs using ngmlr/SNIFFLES/SURVIVOR
4) Call SVs using VarScan
5) Generate statistics that indicate the quality of SV calling
6) Compare the results from each method
7) Generate pipeline for checking Illumina data
8) Validate 10 SVs in Illunina data (as test case)

## PacBio single sample pipeline
Mike automatically runs this pipeline
### VCF files

## PacBio combine-samples pipeline
Mike is trying this out now

## ngmlr/SNIFFLES/SURVIVOR
### Steps
1) ngmlr on each library separately using "SAM attributes enables"
2) sort .bams
3) SNIFFLES on each .bam separately
4) sort vcfs
5) SURVIVOR on all vcfs
6) rerun SNIFFLES to force-call all SVs based on SURVIVOR output
7) sort vcfs
8) Re-run SURVIVOR to get merged vcf
### Analysis Notes
`./ngmlr_sniffles_notes.md`



