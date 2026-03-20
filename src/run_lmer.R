.run_models <- function(taxon, dat) {
  
  print(taxon)
  
  inp <- dat %>%
    select(animal_id, treatment, day, all_of(taxon)) %>%
    dplyr::rename(relab = 4) %>%
    mutate(pres_abs = ifelse(relab > 0, 1, 0)) %>%
    mutate_if(is.character, as.factor)
  
  val <- inp %>%
    filter(relab > 0) %>%
    filter(relab == min(relab)) %>%
    pull(relab) %>%
    unique()
  
  inp$relab <- log(inp$relab + (0.5 * val))
  
  lmm <- lmerTest::lmer(
    relab ~ treatment * day + (1 | animal_id),
    data = inp
  ) %>%
    anova() %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(taxon = taxon) %>%
    select(8, 1, 7) %>%
    dplyr::rename(p.value = 3) %>%
    mutate(method = "LMM")
  
  # run_glmm <- 0.25 * nrow(inp) <= sum(inp$pres_abs) & sum(inp$pres_abs) <= 0.75 * nrow(inp)
  # 
  # if (run_glmm) {
  #   
  #   glmm <- glmmTMB::glmmTMB(
  #     pres_abs ~ treatment * day + (1 | animal_id),
  #     data = inp,
  #     family = binomial
  #   ) %>%
  #     car::Anova() %>%
  #     as.data.frame() %>%
  #     rownames_to_column("term") %>%
  #     mutate(taxon = taxon) %>%
  #     select(5, 1, 4) %>%
  #     rename(p.value = 3) %>%
  #     mutate(method = "GLMM")
  #   
  # } else {
  #   
  #   glmm <- NULL
  #   
  # }
  # 
  # rbind(lmm, glmm) %>%
  #   select(method, taxon, term, p.value)
  
  lmm
  
}

.plot_lmer <- function(feat, covariate, pval, dat, metadata, type = "relab") {
  
  df <- dat %>%
    select(all_of(feat)) %>%
    dplyr::rename(value = 1) %>%
    rownames_to_column("sample_id") %>%
    inner_join(metadata, by = "sample_id")
  
  if (type == "relab") {
    
    y_lab <- "log(abundance [%])"
    df$value <- log1p(df$value)
    title <- feat
    
  } else {
    
    y_lab <- "Shannon index"
    title <- ""
    
    if (covariate == "group") {
      df <- df %>% mutate(group = gsub(" .*", "", group))
    }
  }
  
  if (covariate == "group") {
    
    x_lab <- "Group"
    
    p <- ggplot(df, aes(x = group, y = value)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(fill = group),
                  pch = 21,
                  position = position_jitterdodge(jitter.width = 0.25),
                  show.legend = FALSE)
  } else {
    
    x_lab <- "Timepoint"
    
    df <- df %>%
      mutate(day = as.numeric(gsub(".* ", "", day)))
    
    if (covariate == "group:day") {
      
      p <- ggplot(df, aes(x = day, y = value)) +
        geom_smooth(aes(group = group, colour = group, fill = group),
                    method = "lm")
      
    } else {
      
      p <- ggplot(df, aes(x = day, y = value)) +
        geom_point(pch = 21, fill = "grey") +
        geom_smooth(method = "lm")
      
    }
  }
  
  subtitle <- ifelse(pval < 0.001,
                     "*q* < 0.001",
                     paste0("*q* = ", round(pval, 3)))
  
  p +
    theme_classic(base_size = 12.5) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = ggtext::element_markdown(),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +
    scale_colour_viridis_d() +
    scale_fill_viridis_d() +
    labs(title = title,
         subtitle = subtitle,
         x = x_lab,
         y = y_lab,
         colour = "Group",
         fill = "Group")
  
}

plot_lmer <- function(i, lmer, genus_relab, metadata) {
  
  .plot_lmer(lmer[i, "taxon"],
             lmer[i, "covariate"],
             lmer[i, "p"],
             genus_relab,
             metadata)
  
}

run_lmer <- function(features, taxonomy, metadata, trial_no) {
  
  genus_relab <- get_genus_relab(features, taxonomy, metadata, 0.01, 0.1)
  
  dat <- genus_relab %>%
    rownames_to_column("sample_id") %>%
    inner_join(metadata, by = "sample_id") 
  
  res <- colnames(genus_relab) %>%
    purrr::map(function(x) .run_models(x, dat)) %>%
    bind_rows() %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, "BH"))
  
  dir.create("output/tables", showWarnings = FALSE)
  write_tsv(lmer, sprintf("output/tables/lmer_results.%s.tsv", trial_no))
  
  lmer <- lmer %>%
    filter(p <= 0.25)
  
  plot_list <- 1:nrow(lmer) %>%
    purrr::map(function(x) plot_lmer(x, lmer, genus_relab, metadata))
  
  plot_name <- lmer %>%
    select(1, 2) %>%
    mutate(taxon = taxon %>%
             str_replace_all("\\; ", "\\.") %>%
             str_replace_all("[|]", "") %>%
             str_replace_all(" ", "_")) %>%
    mutate(covariate = str_replace(covariate, "\\:", "X")) %>%
    mutate(filename = sprintf("output/plots/pngs/%s/%s/%s.png", trial_no, covariate, taxon)) %>%
    pull(filename)
  
  plot_name %>%
    purrr::map(function(x) str_split(x, "/") %>% unlist() %>% head(-1) %>% paste0(collapse = "/")) %>%
    unlist() %>%
    unique() %>%
    purrr::map(function(x) dir.create(x, showWarnings = FALSE, recursive = TRUE))
  
  n <- length(plot_list)
  
  purrr::map(1:n, function(x) ggsave(plot_name[x], plot_list[[x]], width = 7.5, height = 5))
  
  filename <- sprintf("output/plots/lmer_plot.%s.pdf", str_to_lower(trial_no))
  pdf(filename, width = 7.5, height = 5)
  print(plot_list)
  dev.off()
  
}
