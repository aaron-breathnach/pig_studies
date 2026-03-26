library(tidyverse)

source("src/get_genus_relab.R")
source("src/get_pal.R")
source("src/make_heatmap.R")
source("src/run_diversity_analysis.R")
source("src/run_lmer.R")

.run_analysis <- function(experiment_id, sample_type, features, taxonomy, metadata) {
  
  print(sprintf("%s|%s", experiment_id, sample_type))
  
  out_dir <- sprintf("output/%s/%s", experiment_id, sample_type)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  metadata <- metadata[which(metadata$sample_type == sample_type),]
  
  run_diversity_analysis(features, taxonomy, metadata, out_dir)
  
  run_lmer(features, taxonomy, metadata, out_dir)
  
}

run_analysis <- function(experiment_id, features, taxonomy, metadata) {
  
  metadata <- metadata %>%
    filter(experiment == experiment_id)
  
  # keep <- metadata %>%
  #   select(sample_type, day, treatment) %>%
  #   distinct() %>%
  #   group_by(sample_type, day) %>%
  #   tally() %>%
  #   ungroup() %>%
  #   filter(n > 1) %>%
  #   select(1, 2)
  # 
  # metadata <- metadata %>%
  #   inner_join(keep, by = c("sample_type", "day"))
  
  sample_types <- unique(metadata$sample_type)
  
  lapply(sample_types, function(sample_type) {
    .run_analysis(experiment_id, sample_type, features, taxonomy, metadata)
  })
  
}

wrapper <- function() {
  
  features <- qiime2R::read_qza("data/qiime/table.qza")$data
  
  taxonomy <- qiime2R::read_qza("data/qiime/taxonomy.qza")$data
  
  experiment_ids <- c(813, 815)
  
  cols <- c(1:5, 5, 7:8)
  
  days <- c(
    "D006",
    "D021",
    "D049",
    "D002",
    "D027",
    "Before farrow",
    "Lactation"
  )
  
  metadata <- read_delim("data/metadata.tsv", col_select = all_of(cols)) %>%
    mutate(animal_id = case_when(
      !is.na(offspring_id) ~ offspring_id,
      !is.na(sow_id) ~ sow_id
    )) %>%
    rename(sample_id = 1, sample_type = 3) %>%
    mutate(sample_id = as.character(sample_id)) %>%
    mutate(sample_type = gsub(" ", "_", sample_type)) %>%
    mutate(treatment = sprintf("Group %s", treatment)) %>%
    select(sample_id, day, sample_type, animal_id, treatment, experiment) %>%
    mutate(day = day %>%
             str_replace("d5_", "") %>%
             str_replace("d_", "") %>%
             str_replace("d", "") %>%
             str_replace("_pw", "") %>%
             str_replace("48_hrs", "2") %>%
             str_pad(3, "left", "0") %>%
             str_replace("_", " ") %>%
             str_to_sentence()
    ) %>%
    mutate(day = ifelse(nchar(day) == 3, paste0("D", day), day)) %>%
    filter(experiment == 815 & day %in% days | experiment == 813)
  
  lapply(experiment_ids, function(experiment_id) {
    run_analysis(experiment_id, features, taxonomy, metadata)
  })
  
}
