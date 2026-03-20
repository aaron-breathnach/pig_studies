library(tidyverse)

source("src/get_genus_relab.R")
source("src/get_pal.R")
source("src/make_heatmap.R")
source("src/run_diversity_analysis.R")
source("src/run_lmer.R")

.run_analysis <- function(features, taxonomy, metadata, out_dir) {
  
  run_lmer(genus_relab, meta)
  
  run_diversity_analysis(genus_relab, meta, trial_no)
  
  
  make_heatmap(genus_relab, meta, trial_no)
  
}

features <- qiime2R::read_qza("data/qiime/table.qza")$data

taxonomy <- qiime2R::read_qza("data/qiime/taxonomy.qza")$data

run_analysis <- function(experiment_id, features, taxonomy) {
  
  if (experiment_id == 813) {
    cols <- c(1:4, 7:8)
  } else {
    cols <- c(1:3, 5, 7:8)
  }
  
  metadata <- read_delim("data/metadata.tsv", col_select = cols) %>%
    rename(sample_id = 1, sample_type = 3, animal_id = 4) %>%
    filter(experiment == experiment_id & !is.na(animal_id)) %>%
    mutate(sample_id = as.character(sample_id)) %>%
    mutate(day = day %>%
             str_replace("d_", "") %>%
             str_replace("d", "") %>%
             str_replace("_pw", "") %>%
             str_replace("48_hrs", "2") %>%
             str_pad(3, "left", "0") %>%
             str_c("D", .)
    ) %>%
    mutate(sample_type = gsub(" ", "_", sample_type)) %>%
    mutate(treatment = sprintf("Group %s", treatment)) %>%
    select(1:5)
  
  sample_types <- unique(metadata$sample_type)
  
  
  
}

.run_analysis <- function(experiment_id, sample_type, features, taxonomy, metadata) {
  
  out_dir <- sprintf("output/%s/%s", experiment_id, sample_type)
  
  metadata <- metadata[which(metadata$sample_type == sample_type),]
  
  run_diversity_analysis(features, taxonomy, metadata, out_dir)
  
  make_heatmap(genus_relab, meta, trial_no)

  metadata <- metadata %>%
    filter(sample_type == sample_types[1]) %>%
    mutate(sample_id = as.character(sample_id))
  
  sample_type <- str_replace_all(sample_type, " ", "_")
  
  out_dir <- sprintf("output/%s/%s", experiment_id, sample_type)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
