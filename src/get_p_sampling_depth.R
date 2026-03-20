#!/usr/bin/Rscript

library(tidyverse)

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)

get_p_sampling_depth <- function(inp, out_dir) {
  
  dat <- read_delim(inp) %>%
    filter(input != "numeric")
  
  for (colname in colnames(dat)[-1]) {
    dat[[colname]] <- as.numeric(dat[[colname]])
  }
  
  p_sampling_depth <- dat %>%
    filter(merged == min(merged)) %>%
    pull(merged) %>%
    as.character()
  
  filename <- sprintf("%s/p_sampling_depth.txt", out_dir)
  writeLines(p_sampling_depth, filename)
  
}

if (sys.nframe() == 0) {
  get_p_sampling_depth(args[1], args[2])
}
