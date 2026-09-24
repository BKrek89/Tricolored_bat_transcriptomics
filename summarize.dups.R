##script to run this on an hpc can be found in QC in bash
#AI assisted in writing this code


library(dupRadar)
library(dplyr)
library(purrr)
library(tibble)

results_dir <- "results"

sample_dirs <- list.dirs(
  results_dir,
  recursive = FALSE,
  full.names = TRUE
)

summary_table <- map_dfr(sample_dirs, function(sample_dir) {
  
  rds_file <- file.path(sample_dir, "DupMat.rds")
  
  if (!file.exists(rds_file)) {
    return(NULL)
  }
  
  dm <- readRDS(rds_file)
  
  fit <- duprateExpFit(dm)
  
  tibble(
    sample = basename(sample_dir),
    intercept = fit$intercept,
    slope = fit$slope
  )
})

summary_table

write.csv(
  summary_table,
  file.path(results_dir, "dupRadar_summary.csv"),
  row.names = FALSE
)
