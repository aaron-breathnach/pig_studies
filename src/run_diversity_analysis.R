alpha_diversity <- function(genus_relab, metadata) {
  
  alpha <- data.frame(vegan::diversity(genus_relab)) %>%
    rownames_to_column("sample_id") %>%
    rename(shannon_index = 2) %>%
    merge(metadata, by = "sample_id") %>%
    as_tibble()
  
  string <- "Interaction *p*-value: %s"
  
  single_timepoint <- length(unique(alpha$day)) == 1
  
  if (single_timepoint) {
    
    pval <- lm(shannon_index ~ treatment, data = alpha) %>%
      anova() %>%
      drop_na() %>%
      pull(5)
    
    xlab <- "Group"
    p <- ggplot(alpha, aes(x = treatment, y = shannon_index))
    
  } else {
    
    pval <- lmerTest::lmer(
      shannon_index ~ treatment * day + (1 | animal_id),
      alpha
    ) %>%
      anova() %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      filter(grepl(":", term)) %>%
      pull(7)
    
    xlab <- "Day"
    p <- ggplot(alpha, aes(x = day, y = shannon_index))
    
  }
  
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
  
  p +
    geom_boxplot(aes(colour = treatment), outlier.shape = NA) +
    geom_jitter(
      aes(colour = treatment, fill = treatment),
      position = position_jitterdodge(jitter.width = 0.25),
      alpha = 0.5
    ) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      plot.title = ggtext::element_markdown()
    ) +
    ggtitle(title) +
    labs(
      x = xlab,
      y = "Shannon index"
    ) +
    scale_colour_manual(values = get_pal("treatment")) +
    scale_fill_manual(values = get_pal("treatment"))
  
}

make_adonis_plot <- function(genus_relab, metadata) {
  
  veg_dis <- vegan::vegdist(genus_relab)
  
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
    mutate(p = rstatix::p_round(p)) %>%
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
    pal <- get_pal("day")
    legend_title <- "Day"
  }
  
  p_inp <- points %>%
    select(Axis.1, Axis.2, all_of(variable), all_of(facet)) %>%
    rename(variable = 3, facet = 4)
  
  p <- ggplot(p_inp, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = variable), pch = 21, size = 3) +
    theme_bw(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(x = xlab, y = ylab, fill = legend_title) +
    scale_fill_manual(values = pal)
  
  n <- length(unique(points$day))
  
  if (n > 1) {
    p <- p +
      facet_wrap(~ facet, scales = "free", ncol = 2)
  }
  
  return(p)
  
}

viz_int_pig_dis <- function(genus_relab, metadata) {
  
  veg_dis <- vegan::vegdist(genus_relab)
  
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
  
  n <- length(unique(p_inp$day))
  
  if (n == 1) {
    
    p <- ggplot(p_inp, aes(x = treatment, y = distance)) +
      geom_boxplot(aes(colour = treatment), outlier.shape = NA) +
      geom_jitter(
        aes(colour = treatment, fill = treatment)
      )
    
    xlab <- "Group"
    
  } else {
    
    p <- ggplot(p_inp, aes(x = day, y = distance)) +
      geom_boxplot(aes(colour = treatment), outlier.shape = NA) +
      geom_jitter(
        aes(colour = treatment, fill = treatment),
        position = position_jitterdodge(jitter.width = 0.25)
      )
    
    xlab <- "Day"
    
  }
  
  p +
    scale_colour_manual(values = get_pal("treatment")) +
    scale_fill_manual(values = get_pal("treatment")) +
    labs(x = xlab, y = "Bray-Curtis distance") +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
}

make_pcoa_plots <- function(genus_relab, metadata) {
  
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
  
  n <- length(unique(points$day))
  
  if (n > 1) {
    
    p1 <- make_pcoa_plot(points, "treatment", "day", xlab, ylab)
    p2 <- make_pcoa_plot(points, "day", "treatment", xlab, ylab)
    
    if (n > 2) {
      
      p <- patchwork::wrap_plots(
        p1,
        p2,
        nrow = 2,
        heights = c(2, 1)
      )
      
    } else {
      
      p <- patchwork::wrap_plots(
        patchwork::plot_spacer(),
        p1,
        p2,
        patchwork::plot_spacer(),
        nrow = 4,
        heights = c(1, 2, 2, 1)
      )
      
    }
    
  } else {
    
    data <- metadata %>%
      column_to_rownames("sample_id")
    
    permanova <- vegan::adonis2(
      veg_dis ~ treatment,
      data = data
    ) %>%
      {.[1, 3:4]} %>%
      as.numeric()
    
    title <- sprintf(
      "PERMANOVA: *R*<sup>2</sup>=%s, *p*=%s",
      round(permanova[1], 2),
      round(permanova[2], 3)
    )
    
    p <- make_pcoa_plot(points, "treatment", "day", xlab, ylab) +
      theme(
        axis.line = element_line(),
        panel.border = element_blank(),
        plot.title = ggtext::element_markdown()
      ) +
      ggtitle(title)
    
  }
  
  return(p)
  
}

run_diversity_analysis <- function(features, taxonomy, metadata, out_dir) {
  
  genus_relab <- get_genus_relab(features, taxonomy, metadata)
  
  filename <- sprintf("%s/diversity.png", out_dir)
  
  p1 <- alpha_diversity(genus_relab, metadata)
  
  if (length(unique(metadata$day)) > 1) {
    
    p2 <- make_adonis_plot(genus_relab, metadata)
    p3 <- viz_int_pig_dis(genus_relab, metadata)
    p4_p5 <- make_pcoa_plots(genus_relab, metadata)
    p1_p2_p3 <- patchwork::wrap_plots(p1, p2, p3, nrow = 3)
    
    p <- patchwork::wrap_plots(
      p1_p2_p3,
      p4_p5,
      nrow = 1,
      widths = c(1, 2)
    ) +
      patchwork::plot_annotation(tag_levels = "A") +
      patchwork::plot_layout(guides = "collect")
    
    w <- 15
    h <- 12.5
    
  } else {
    
    p2 <- make_pcoa_plots(genus_relab, metadata)
    
    p <- patchwork::wrap_plots(p1, p2)
    
    w <- 10
    h <- 5
    
  }
  
  ggsave(filename, p, width = w, height = h, bg = "white")
  
}
