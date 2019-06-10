# Analysis comparing genetic distance and ages of branches and nodes
* Note: I developed this for result that were biased by a sample-order bias \
in a previous version of pbsv; There were too few variable SVs to do \
this analyis with the updated results

## Overview
### Goals
* Model relationship between SVs and age of branches
* Estimate age of node of bottom branch on Tree 14
### Approach
* Generate different models to explain relationship between SVs and age
  * Model 1: Use Branch and Node Ages
  * Model 2: Used fixed age difference between Tree 13 and Tree 14
  * Model 3: Use only Node Ages
* Compare SV-based genetic distance and age-based distance
  * Mantel Tests
  * Linear models
### Analysis Files
* R script with analysis steps
  * `~/pb_dist_comp_analysis.r`
* Powerpoint explaining analysis and results
  * `~/ppt_files/Branch_age_vs_genetic_distance.pptx`
* Conda environment
  * `/home/grabowsky/.conda/envs/NGS_analysis`
 
