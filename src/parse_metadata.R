library(tidyverse)

setwd("/Users/aaron.walsh/Desktop/pig_studies/")

xlsx <- "data/metadata_813_and_815.xlsx"

metadata <- readxl::read_excel(xlsx, skip = 2) %>%
  rename("sample-id" = 1)

colnames(metadata) <- colnames(metadata) %>%
  str_replace_all(" ", "_") %>%
  str_replace_all("__", "_") %>%
  str_to_lower() %>%
  str_replace("/", "_")

tsv <- str_replace(xlsx, ".xlsx", ".tsv")

write_tsv(metadata, tsv)
