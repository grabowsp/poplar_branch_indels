# Effect of including all 10 Branches
## Overview
### Background Info
* Branches 13.4 and 14.1 have been recommended to be removed because they showed signes of disease which could affect genomic results
* Previous analysis shows that samples 13.5 and 14.5 were switched at some point in the data generation process
* Perhaps the elevated numbers of unique SVs in 13.2 and 14.3 are actually a result of sample switches with 13.4 and/or 14.1.
### Approach
* Calculate SV-based genetic distances and use to generate trees
* Calculate numbers of sample-specific and clone-specific SVs
### Analysis Files
* R script
  * `~/r_scripts/pb_allBranch_SV_analysis.r`
* PowerPoint with explanations of analysis and results
  * `~/ppt_files/Ten_Branch_SV_analysis.pptx`
* Conda environment
  * `/home/grabowsky/.conda/envs/NGS_analysis`

