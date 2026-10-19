This repository contains the code for all analyses in the manuscript "Preliminary evidence suggests the transcriptomic response of tricolored bats to white-nose syndrome is similar but not identical to that of little brown bats"

Part of the included code was run on a linux based HPC and part was completed in R.

Required Bioinformatics packages
   1. BARCODESPLITTER
   2. bbtools
   3. TRIMMOMATIC
   4. FastQScreen
   5. STAR
   6. samtools
   7. HTSeq
   8. Subread
   9. Picard
   10. genmap
   11. bedtools
    
Required R packages
   1. dplyr
   2. ggplot2
   3. ggbiplot
   4. tidyverse
   5. DESeq2
   6. dupRadar
   7. emmeans
   8. DHARMa
   9. writexl
   10. eulerr
   11. tidyr
   12. BiocManager
   13. purrr
   14. tibble





The order of the scripts is as follows:

On the HPC

    1. genome_prep
    2. read filtering and alignment in Bash

         a. move to main analysis in R
         b. move to QC
         c. move to P. destructans analysis

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

On the HPC

    1. Pd in Bash

In R

    1. pd.R
