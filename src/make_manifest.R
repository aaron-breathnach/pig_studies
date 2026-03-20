#!/usr/bin/Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

get_abs_filepaths <- function(inp_dir, sample_id) {
  
  r1 <- normalizePath(sprintf("%s/%s_R1_001.fastq.gz", inp_dir, sample_id))
  r2 <- normalizePath(sprintf("%s/%s_R2_001.fastq.gz", inp_dir, sample_id))
  
  cols <- c("sample-id", "forward-absolute-filepath", "reverse-absolute-filepath")
  
  tibble(sample_id, r1, r2) %>%
    setNames(cols)
  
}

make_manifest <- function(inp_dir, out_dir) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  sample_ids <- list.files(inp_dir, pattern = "R1", full.names = TRUE) %>%
    purrr::map(function(x) x %>% basename() %>% str_split("_R") %>% unlist() %>% nth(1)) %>%
    unlist()
  
  manifest <- purrr::map(sample_ids, function(x) get_abs_filepaths(inp_dir, x)) %>%
    bind_rows()
  
  filename <- sprintf("%s/manifest.tsv", out_dir)
  write_tsv(manifest, filename)
  
}

if (sys.nframe() == 0) {
  make_manifest(args[1], args[2])
}
