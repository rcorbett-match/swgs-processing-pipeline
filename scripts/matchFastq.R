# Required Packages
library(dplyr)
library(tidyr)
library(glue)
library(stringr)

# Access command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 1) {
  stop("Please provide the required input.")
} else {
  samples_path <- args[1]
  reads_path <- args[2]
  paired <- args[3]
}

if (grepl("*", reads_path, fixed = TRUE)) {
    stop('`reads` must not contain wildcards.')
}

samples <- read.csv(samples_path)

# prep sample to fastq path df
if (paired == "true") {
  cols <- c("sample_id", "read1", "read2")
} else {
  cols <- c("sample_id", "read")
}
sample_fastq <- data.frame(matrix(ncol=length(cols),nrow=0, dimnames=list(NULL, cols))) %>%
  mutate_all(as.character)

for (i in 1:nrow(samples)) {
  sample_id <- samples$sample_id[i]
  reads_id <- samples$reads_id[i]

  los <- list.files(path = reads_path, full.names = TRUE, pattern = glue(".*({reads_id}).*"), recursive = TRUE) %>%
    normalizePath()
  row <- data.frame(t(c(sample_id, los)))
  if (ncol(row) != length(cols)) {
    warning(glue("{sample_id} did not match to the correct number of fastq file(s) (expects {length(cols)-1}, but matched {ncol(row) - 1})"))
    next
  }
  colnames(row) <- cols
  sample_fastq <- bind_rows(sample_fastq, row)
}

write.csv(sample_fastq, glue("sample_fastq.csv"), row.names = FALSE, quote = FALSE)