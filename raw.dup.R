##script to run this on an hpc can be found in QC in bash
#AI assisted in writing this R code

out_dir <- "results"

#create sample file list
dup_files <- list.files(
  out_dir,
  pattern = "DupMat\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Found", length(dup_files), "DupMat files\n")

if (length(dup_files) == 0) {
  stop("No DupMat.tsv files found.")
}

results <- list()

#import dupradar files for each sample

for (file in dup_files) {
  
  sample <- basename(dirname(file))
  
  cat("Processing:", sample, "\n")
  dm <- read.delim(
    file,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  dm <- dm[
    dm$allCounts > 0 &
      !is.na(dm$dupRate) &
      is.finite(dm$dupRate),
  ]
  
#log transform counts
  
  log10_counts <- log10(dm$allCounts)
  
#Regression
  
  fit <- lm(
    dm$dupRate ~ log10_counts
  )
  
  coefficients <- coef(fit)
  
  intercept <- coefficients[1]
  slope <- coefficients[2]
  
  model_summary <- summary(fit)
  
  r_squared <- model_summary$r.squared
  
  slope_p <- model_summary$coefficients[
    "log10_counts",
    "Pr(>|t|)"
  ]
  dupRate_1read <- predict(
    fit,
    newdata = data.frame(
      log10_counts = log10(1)
    )
  )
  
#extract predicted duplication rate at 10 reads 
  dupRate_10reads <- predict(
    fit,
    newdata = data.frame(
      log10_counts = log10(10)
    )
  )
  

  
  results[[sample]] <- data.frame(
    sample = sample,
    n_genes = nrow(dm),
    intercept = intercept,
    slope = slope,
    dupRate_1read = as.numeric(dupRate_1read),
    dupRate_10reads = as.numeric(dupRate_10reads),
    R_squared = r_squared,
    slope_p = slope_p
  )
}

#make table and export

results_df <- do.call(
  rbind,
  results
)

rownames(results_df) <- NULL
output_file <- file.path(
  out_dir,
  "dupRadar_raw_count_regression_results.tsv"
)

write.table(
  results_df,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
