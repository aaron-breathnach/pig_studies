#!/usr/bin/Rscript

library(tidyverse)

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)

get_cutoff <- function(filename) {
  
  read_delim(filename)[4, -1] %>%
    pivot_longer(cols = everything()) %>%
    dplyr::rename(pos = 1, qual = 2) %>%
    mutate_all(as.numeric) %>%
    filter(qual < 30) %>%
    filter(pos == min(pos)) %>%
    pull(pos) %>%
    as.character()
  
}

get_cutoffs <- function(inp_dir, out_dir) {
  
  filenames <- list.files(
    inp_dir,
    pattern = "-seven-number-summaries.tsv", 
    full.names = TRUE
  )
  
  cutoffs <- unlist(purrr::map(filenames, function(x) get_cutoff(x)))
  
  filename <- sprintf("%s/trunc_len.txt", out_dir)
  
  writeLines(cutoffs, filename)
  
}

if (sys.nframe() == 0) {
  get_cutoffs(args[1], args[2])
}
