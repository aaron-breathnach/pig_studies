alpha_diversity <- function(genus_relab, metadata) {
  
  alpha <- data.frame(vegan::diversity(genus_relab)) %>%
    rownames_to_column("sample_id") %>%
    rename(shannon_index = 2) %>%
    merge(metadata, by = "sample_id") %>%
    as_tibble()
  
  pval <- lmerTest::lmer(
    shannon_index ~ treatment * day + (1 | animal_id),
    alpha
  ) %>%
    anova() %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(grepl(":", term)) %>%
    pull(7)
  
  string <- "Interaction *p*-value = %s"
  
  if (pval >= 0.001) {
    title <- sprintf(
      string,
      round(pval, 3)
    )
  } else {
    title <- sprintf(
      string,
      scales::scientific(pval, 3)
    )
  }
  
  ggplot(alpha, aes(x = day, y = shannon_index)) +
    geom_boxplot(aes(colour = treatment), outlier.shape = NA) +
    geom_jitter(
      aes(colour = treatment, fill = treatment),
      position = position_jitterdodge(jitter.width = 0.25),
      alpha = 0.5
    ) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.title = ggtext::element_markdown()
    ) +
    ggtitle(title) +
    labs(
      x = "Day",
      y = "Shannon index",
      fill = "Treatment",
      colour = "Treatment"
    ) +
    scale_colour_manual(values = get_pal("treatment")) +
    scale_fill_manual(values = get_pal("treatment"))
  
}

make_adonis_plot <- function(veg_dis, metadata) {
  
  data <- metadata %>%
    column_to_rownames("sample_id") %>%
    {.[labels(veg_dis),]}
  
  p_inp <- vegan::adonis2(
    veg_dis ~ treatment + day + treatment * day,
    data = data,
    strata = data$animal_id,
    by = "terms"
  ) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    select(1, 4, 6) %>%
    drop_na() %>%
    dplyr::rename(p = 3) %>%
    mutate(term = term %>%
             str_replace("treatment", "Group") %>%
             str_replace("day", "Day")) %>%
    arrange(R2) %>%
    mutate(term = factor(term, levels = .$term))
  
  ggplot(p_inp, aes(x = R2, y = term)) +
    geom_bar(stat = "identity", colour = "black", fill = "#CB7A5C", width = 0.75) +
    geom_text(
      aes(label = sprintf('paste(italic("p")," = %s")', p)),
      hjust = -0.1,
      parse = TRUE
    ) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown()
    ) +
    labs(x = "***R*<sup>2</sup>**", y = "**Term**") +
    scale_x_continuous(expand = expansion(mult = c(0, 1/3)))
  
}

make_pcoa_plot <- function(points, variable, facet, xlab, ylab, legend_title) {
  
  if (variable == "treatment") {
    pal <- get_pal("treatment")
    legend_title <- "Group"
  } else {
    n <- length(unique(points$day))
    pal <- viridis::viridis(n)
    legend_title = "Day"
  }
  
  p_inp <- points %>%
    select(Axis.1, Axis.2, all_of(variable), all_of(facet)) %>%
    rename(variable = 3, facet = 4)
  
  ggplot(p_inp, aes(x = Axis.1, y = Axis.2)) +
    facet_wrap(~ facet, scales = "free") +
    geom_point(aes(fill = variable), pch = 21, size = 3) +
    theme_bw(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(x = xlab, y = ylab, fill = "Group") +
    scale_fill_manual(values = pal)
  
}

viz_int_pig_dis <- function(veg_dis, metadata) {
  
  remove <- metadata %>%
    group_by(animal_id, day) %>%
    tally() %>%
    filter(n > 1) %>%
    pull(animal_id) %>%
    as.character()
  
  data <- metadata %>%
    filter(!animal_id %in% remove) %>%
    mutate(day_num = as.numeric(as.factor(day))) %>%
    select(sample_id, animal_id, day, day_num, treatment)
  
  cols <- c("treatment", "animal_id", "day", "distance")
  
  p_inp <- veg_dis %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("sample1") %>%
    pivot_longer(!sample1, names_to = "sample2", values_to = "distance") %>%
    as_tibble() %>%
    inner_join(data, by = c("sample1" = "sample_id")) %>%
    inner_join(data, by = c("sample2" = "sample_id")) %>%
    filter(animal_id.x == animal_id.y) %>%
    filter(day_num.y - day_num.x == 1) %>%
    arrange(animal_id.x, day.y) %>%
    select(treatment.x, animal_id.x, day.y, distance) %>%
    setNames(cols)
  
  ggplot(p_inp, aes(x = day, y = distance)) +
    geom_boxplot(aes(colour = treatment), outlier.shape = NA) +
    geom_jitter(
      aes(colour = treatment, fill = treatment),
      position = position_jitterdodge(jitter.width = 0.25)
    ) +
    scale_colour_manual(values = get_pal("treatment")) +
    scale_fill_manual(values = get_pal("treatment")) +
    labs(x = "Day", y = "Bray-Curtis distance", colour = "Group", fill = "Group") +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

beta_diversity <- function(genus_relab, metadata) {
  
  veg_dis <- vegan::vegdist(genus_relab)
  
  ## PCoA
  pcoa <- ape::pcoa(veg_dis)
  
  rel_eig <- pcoa$values$Relative_eig[1:2]
  
  xlab <- paste0("PCoA1 (", round(100 * rel_eig[1], 2), "%)")
  ylab <- paste0("PCoA2 (", round(100 * rel_eig[2], 2), "%)")
  
  points <- pcoa$vectors[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    select(1:3) %>%
    inner_join(metadata, by = "sample_id")
  
  p1 <- make_pcoa_plot(points, "treatment", "day", xlab, ylab)
  p2 <- make_pcoa_plot(points, "day", "treatment", xlab, ylab)
  
  ## adonis bar chart
  p3 <- make_adonis_plot(veg_dis, metadata)
  
  ## Bray-Curtis pairwise distances
  p4 <- viz_int_pig_dis(veg_dis, metadata)
  
  ## combine the plots
  
  top <- patchwork::wrap_plots(
    p1,
    p2,
    heights = c(2, 1)
  )
  
  bottom <- patchwork::wrap_plots(
    p3,
    p4,
    nrow = 1
  )
  
  patchwork::wrap_plots(top, bottom, heights = c(3, 1)) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A")
  
}

run_diversity_analysis <- function(features, taxonomy, metadata, out_dir) {
  
  genus_relab <- get_genus_relab(features, taxonomy, metadata)
  
  filenames <- sprintf("%s/%s_diversity.png", out_dir, c("alpha", "beta"))
  
  p1 <- alpha_diversity(genus_relab, metadata)
  ggsave(filenames[1], p1, width = 7.5, height = 5, bg = "white")
  
  p2 <- beta_diversity(genus_relab, metadata)
  ggsave(filenames[2], p2, width = 10, height = 15, bg = "white")
  
}
