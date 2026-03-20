tidy_genus_name <- function(x) {
  
  taxon_levels <- x %>%
    str_split("; ") %>%
    unlist()
  
  for (i in 6:1) {
    taxon <- taxon_levels[i]
    if (!grepl("__$", taxon) & !grepl("\\d", taxon)) {
      break
    }
  }
  
  taxon <- paste0(taxon_levels[i:6], collapse = "; ")
  
  return(taxon)
  
}

get_genus_relab <- function(features, taxonomy, metadata, abundance = 0, prevalence = 0) {
  
  tax <- taxonomy %>%
    dplyr::select(1, 2)
  
  genus_relab <- features[,metadata$sample_id] %>%
    as.data.frame() %>%
    rownames_to_column("Feature.ID") %>%
    inner_join(tax, by = "Feature.ID") %>%
    pivot_longer(cols = !c(Feature.ID, Taxon), names_to = "sample_id", values_to = "count") %>%
    filter(count > 0 & grepl("g__", Taxon)) %>%
    mutate(genus = gsub("\\; s__.*", "", Taxon)) %>%
    group_by(sample_id, genus) %>%
    summarise(count = sum(count)) %>%
    group_by(sample_id) %>%
    mutate(relab = 100 * count / sum(count)) %>%
    ungroup() %>%
    select(1, 2, 4)
  
  genera <- genus_relab %>%
    filter(relab > abundance) %>%
    select(genus, sample_id) %>%
    distinct() %>%
    group_by(genus) %>%
    tally() %>%
    filter(n > prevalence * nrow(metadata)) %>%
    pull(genus)
  
  genus_relab <- genus_relab %>%
    filter(genus %in% genera) %>%
    pivot_wider(names_from = genus, values_from = relab, values_fill = 0) %>%
    column_to_rownames("sample_id")
  
  colnames(genus_relab) <- colnames(genus_relab) %>%
    purrr::map(function(x) tidy_genus_name(x)) %>%
    unlist()
  
  return(genus_relab)
  
}
