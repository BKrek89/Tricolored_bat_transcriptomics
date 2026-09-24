library(dupRadar)

#set up paths with config file

config_file <- commandArgs(trailingOnly = TRUE)[1]

if (is.na(config_file) || config_file == "") {
  stop("Usage: Rscript run_dupRadar.R config.txt")
}

config <- readLines(config_file)
config <- config[!grepl("^\\s*#", config)]
config <- config[nzchar(trimws(config))]

cfg <- list()

for (line in config) {
  parts <- strsplit(line, "=", fixed = TRUE)[[1]]
  
  key <- trimws(parts[1])
  value <- trimws(paste(parts[-1], collapse = "="))
  
  cfg[[key]] <- value
}

bam_dir  <- cfg$BAM_DIR
gtf      <- cfg$GTF
out_dir  <- cfg$OUT_DIR
tmp_dir  <- cfg$TMP_DIR

stranded <- as.integer(cfg$STRANDED)
paired   <- as.logical(cfg$PAIRED)
threads  <- as.integer(cfg$THREADS)


#create output directories

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)


#import bam files
bam_files <- list.files(
  bam_dir,
  pattern = "\\.bam$",
  full.names = TRUE
)

if (length(bam_files) == 0) {
  stop("No BAM files found in: ", bam_dir)
}

cat("Found", length(bam_files), "BAM files\n")

#get paths and create directories for each sample

for (bam in bam_files) {
  
  sample <- basename(bam)
  sample <- sub("\\.bam$", "", sample)
  
  cat("\n========================================\n")
  cat("Processing:", sample, "\n")
  cat("========================================\n")
  
  sample_out <- file.path(out_dir, sample)
  
  dir.create(
    sample_out,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  #Run dupRadar
  
  cat("Running analyzeDuprates...\n")
  
  dm <- analyzeDuprates(
    bam = bam,
    gtf = gtf,
    stranded = stranded,
    paired = paired,
    threads = threads,
    verbose = TRUE,
    tmpDir = tmp_dir
  )

  #Save files
  
  saveRDS(
    dm,
    file.path(sample_out, "DupMat.rds")
  )
  
  write.table(
    dm,
    file = file.path(sample_out, "DupMat.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  #density plot
  
  pdf(
    file.path(sample_out, "dupRadar_plots.pdf"),
    width = 8,
    height = 7
  )
  
  try(
    duprateExpDensPlot(dm),
    silent = TRUE
  )
  