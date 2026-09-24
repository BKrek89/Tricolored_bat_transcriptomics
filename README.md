This repository contains the code for all analyses in the manuscript "Preliminary evidence suggests the transcriptomic response of tricolored bats to white-nose syndrome is similar but not identical to that of little brown bats"

Part of the included code was run on a linux based HPC and part was completed in R.

The order of the scripts is as follows:

On the HPC

    1. read filtering and alignment in Bash

         a. move to main analysis in R
         b. move to QC
         c. mov to P. destructans analysis

Main analysis in R
 
    1. post_alignment_filtering.R
    2. PC_regression.R
    3. DESeq2_and_venndiagram.R


QC (after completing read filtering and alignment in Bash)

  On the HPC

    1. QC in bash

  In R on the HPC
  
    1. run_dupradar.R
    2. summarize.dups.R
    3. raw.dup.R

  In R
  
    1. qc_graphs.R
    2. dupradar.plots.R


P. destructans analysis

