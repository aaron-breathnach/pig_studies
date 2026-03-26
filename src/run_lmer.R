.run_models <- function(taxon, dat) {
  
  print(taxon)
  
  inp <- dat %>%
    select(animal_id, treatment, day, all_of(taxon)) %>%
    dplyr::rename(relab = 4) %>%
    mutate_if(is.character, as.factor)
  
  min_val <- inp %>%
    filter(relab > 0) %>%
    filter(relab == min(relab)) %>%
    pull(relab)
  
  inp$relab <- log(inp$relab + 0.5 * min_val)
  
  if (length(unique(dat$day)) == 1) {
    
    mod <- lm(relab ~ treatment, data = inp)
    
    cols <- c(7, 1, 6)
    
    method <- "LM"
    
  } else {
    
    mod <- lmerTest::lmer(
      relab ~ treatment * day + (1 | animal_id),
      data = inp
    )
    
    cols <- c(8, 1, 7)
    
    method <- "LMM"
    
  }
  
  mod %>%
    anova() %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(taxon = taxon) %>%
    select(all_of(cols)) %>%
    dplyr::rename(p.value = 3) %>%
    mutate(method = method) %>%
    drop_na()
  
}

viz_sig_spp <- function(sp, dat, sig_res) {
  
  p_inp <- dat %>%
    select(treatment, day, all_of(sp)) %>%
    rename(relab = 3)
  
  min_val <- p_inp %>%
    filter(relab > 0) %>%
    filter(relab == min(relab)) %>%
    pull(relab)
  
  p_inp$relab <- log(p_inp$relab + 0.5 * min_val)
  
  subtitle <- sig_res %>%
    filter(taxon == sp) %>%
    pull(fdr) %>%
    rstatix::p_round() %>%
    sprintf("*q*-value: %s", .)
  
  n <- length(unique(p_inp$day))
  
  if (n > 1) {
    
    p <- ggplot(p_inp, aes(x = day, y = relab)) +
      geom_boxplot(aes(colour = treatment), outlier.shape = NA) +
      geom_jitter(
        aes(colour = treatment, fill = treatment),
        position = position_jitterdodge(jitter.width = 0.25),
        alpha = 0.5,
      )
    
  } else {
    
    p <- ggplot(p_inp, aes(x = treatment, y = relab)) +
      geom_boxplot(aes(colour = treatment), outlier.shape = NA) +
      geom_jitter(
        aes(colour = treatment, fill = treatment),
        alpha = 0.5,
      )
    
  }
  
  p +
    ggtitle(sp) +
    labs(
      x = xlab,
      y = "**log<sub>10</sub>(abundance)**",
      colour = "Group",
      fill = "Group",
      subtitle = subtitle
    ) +
    scale_colour_manual(values = get_pal("treatment")) +
    scale_fill_manual(values = get_pal("treatment")) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title.x = element_text(face = "bold"),
      axis.title.y = ggtext::element_markdown(),
      legend.title = element_text(face = "bold"),
      plot.subtitle = ggtext::element_markdown()
    )
  
}

run_lmer <- function(features, taxonomy, metadata, out_dir) {
  
  genus_relab <- get_genus_relab(features, taxonomy, metadata, 0.01, 0.1, clr = FALSE)
  
  dat <- genus_relab %>%
    rownames_to_column("sample_id") %>%
    inner_join(metadata, by = "sample_id")
  
  res <- colnames(genus_relab) %>%
    purrr::map(function(x) .run_models(x, dat)) %>%
    bind_rows() %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, "BH")) %>%
    ungroup()
  
  tsv <- sprintf("%s/stats.tsv", out_dir)
  
  write_tsv(res, tsv)
  
  sig_res <- res %>%
    filter(grepl(":", term) & fdr <= 0.05)
  
  if (nrow(sig_res > 0)) {
    
    plot_list <- lapply(sig_res$taxon, function(x) viz_sig_spp(x, dat, sig_res))
    
    filename <- sprintf("%s/sig_res.pdf", out_dir)
    pdf(filename, width = 7.5, height = 5)
    print(plot_list)
    dev.off()
    
  }
  
}
