suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(future)
  library(CGHcall)
})

# Access command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 1) {
  stop("Please provide the required input.")
} else {
  nthreads <- args[1]
  type <- args[2]
  dir_path <- args[3]
  genome <- args[4]
  out_path <- args[5]
}

binsize <- as.numeric(paste0(stringr::str_remove(basename(out_path), 'kb'), "000"))

future::plan("multisession", workers=as.integer(nthreads))

# Convert WX files to QDNAseq object
if (type == "wisecondorx") {
  cnobj <- utanos::BuildWxQdnaObject(dir_path, genome_used = genome, 
                                       bin_size = binsize, Xincluded = TRUE)
  
  # find solutions
  solutions <- utanos::FindRascalSolutions(cnobj)
  
  # scaling to ACN
  results <- utanos::CalculateACNs(cnobj, "mad", solutions, return_sols = TRUE)
  
  results[["all_solutions"]] <- solutions
  saveRDS(results, glue::glue("{out_path}/wx_acn_{binsize/1000}kb.rds"))
} else if (type == "qdnaseq") {
  for (i in c("noY","noXY")) {
    cnobj <- readRDS(glue::glue("{dir_path}/{basename(dir_path)}_{i}_copyNumbersSegmented.rds"))
    
    # find solutions
    solutions <- utanos::FindRascalSolutions(cnobj)
    
    # scaling to ACN
    results <- utanos::CalculateACNs(cnobj, "mad", solutions, return_sols = TRUE)
    
    results[["all_solutions"]] <- solutions
    saveRDS(results, glue::glue("{out_path}/qdnaseq_{i}_acn_{binsize/1000}kb.rds"))
  }
}
