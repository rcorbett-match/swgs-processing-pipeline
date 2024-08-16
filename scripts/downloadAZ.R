# Required Packages
library(AzureStor)
library(processx)
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
  sas <- args[1]
  cont_url <- args[2]
  samples_path <- args[3]
  paired <- args[4]
  out_dir <- args[5]
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

# loop over each sample and download fastqs
for (i in 1:nrow(samples)) {
  sample_id <- samples$sample_id[i]
  az_url <- samples$az_url[i]

  dir <- str_remove(az_url, "(https://)?.+/.+/")

  AzureStor::call_azcopy("copy",
                          glue("{cont_url}/*?{sas}"),
                          out_dir,
                          "--recursive=true",
                          '--include-regex', ".*\\.(fastq\\.gz|fq\\.gz|fastq|fq)",
                          "--include-path", dir)

  file.rename(glue("{out_dir}/{dir}"), glue("{out_dir}/{sample_id}"))
  los <- list.files(path = glue("{out_dir}/{sample_id}"), full.names = TRUE) %>%
    normalizePath()
  row <- data.frame(t(c(sample_id, los)))
  if (ncol(row) != length(cols)) {
    warning(glue("{sample_id} did not match to the correct number of fastq file(s) (expects {length(cols)-1}, but matched {ncol(row) - 1})"))
    next
  }
  colnames(row) <- cols
  sample_fastq <- bind_rows(sample_fastq, row)
}

write.csv(sample_fastq, glue("{out_dir}/sample_fastq.csv"), row.names = FALSE, quote = FALSE)