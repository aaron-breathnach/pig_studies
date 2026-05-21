library(tidyverse)

source("src/get_genus_relab.R")
source("src/get_pal.R")

#############################
## beta diversity analysis ##
#############################

run_adonis2 <- function(veg_dis, metadata) {
  
  meta <- metadata %>%
    column_to_rownames("sample_id")
  
  adonis2 <- vegan::adonis2(
    veg_dis ~ timepoint + hygiene + pigs_per_pen + feeders_per_pen,
    data = meta,
    permutations = 999,
    strata = meta$pig_id,
    by = "margin"
  ) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    drop_na() %>%
    mutate(
      term = term %>%
        str_replace_all("_", " ") %>%
        str_to_sentence()
    ) %>%
    select(term, 4, 6) %>%
    rename(r2 = 2, p = 3) %>%
    arrange(r2) %>%
    mutate(term = factor(term, levels = term))
  
  p <- ggplot(adonis2, aes(x = r2, y = term)) +
    geom_point(colour = "darkred") +
    geom_segment(aes(xend = 0, yend = term), colour = "darkred") +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(x = expression(bold(R^2))) +
    geom_text(aes(label = sprintf("p: %s", p)), hjust = -0.1) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.25)))
  
  list(table = adonis2, plot = p)
}

make_pcoa_plot <- function(pcoa_df, pcoa, var) {
  
  p_inp <- pcoa_df %>%
    select(Axis.1, Axis.2, all_of(var), timepoint) %>%
    rename(variable = 3)
  
  legend_title <- case_when(
    var == "hygiene" ~ "Hygiene",
    var == "feeders_per_pen" ~ "Feeders per pen",
    var == "pigs_per_pen" ~ "Pigs per pen",
    TRUE ~ var
  )
  
  xlab <- sprintf("PCoA 1 (%.2f%%)", pcoa$values$Relative_eig[1] * 100)
  ylab <- sprintf("PCoA 2 (%.2f%%)", pcoa$values$Relative_eig[2] * 100)
  
  ggplot(p_inp, aes(x = Axis.1, y = Axis.2)) +
    facet_wrap(~ timepoint, scales = "free") +
    geom_point(aes(fill = variable), pch = 21) +
    stat_ellipse(
      aes(colour = variable),
      level = 0.95,
      linetype = "dashed",
      show.legend = FALSE
    ) +
    theme_bw(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      x = xlab,
      y = ylab,
      colour = legend_title,
      fill = legend_title
    ) +
    scale_colour_manual(values = get_pal(var)) +
    scale_fill_manual(values = get_pal(var))
}

#####################################
## differential abundance analysis ##
#####################################

get_term_info <- function(term) {
  case_when(
    term == "hygiene" ~ "Hygiene",
    term == "feeders_per_pen" ~ "Feeders per pen",
    term == "pigs_per_pen" ~ "Pigs per pen",
    TRUE ~ term
  )
}

.viz_sig_gen <- function(gns, term, fdr, genus_relab, metadata, relab = TRUE, ann) {
  
  var <- term
  legend_title <- get_term_info(term)
  
  p_inp <- genus_relab %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    inner_join(metadata, by = "sample_id") %>%
    select(sample_id, timepoint, all_of(var), all_of(gns)) %>%
    rename(variable = 3, relab = 4)
  
  if (relab) {
    min_val <- p_inp %>%
      filter(relab > 0) %>%
      summarise(min_val = min(relab, na.rm = TRUE)) %>%
      pull(min_val)
    
    if (!is.finite(min_val)) {
      min_val <- 1e-06
    }
    
    p_inp <- p_inp %>%
      mutate(relab = log10(relab + 0.5 * min_val))
    
    title <- gns
    y_lab <- "**log<sub>10</sub>(abundance)**"
    
  } else {
    title <- NULL
    y_lab <- "Shannon index"
    ann <- ann %>%
      mutate(fdr = rstatix::p_format(p.value))
  }
  
  ann <- ann %>%
    mutate(relab = max(p_inp$relab, na.rm = TRUE) + 0.05 * diff(range(p_inp$relab, na.rm = TRUE)))
  
  ggplot(p_inp, aes(x = timepoint, y = relab)) +
    geom_boxplot(aes(colour = variable), outlier.shape = NA) +
    geom_jitter(
      aes(colour = variable, fill = variable),
      position = position_jitterdodge(jitter.width = 0.25),
      alpha = 0.5
    ) +
    geom_text(
      data = ann,
      aes(x = timepoint, y = relab, label = fdr),
      inherit.aes = FALSE
    ) +
    ggtitle(title) +
    labs(
      x = "Timepoint",
      y = y_lab,
      colour = legend_title,
      fill = legend_title,
      subtitle = fdr
    ) +
    scale_colour_manual(values = get_pal(var)) +
    scale_fill_manual(values = get_pal(var)) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title.x = element_text(face = "bold"),
      axis.title.y = ggtext::element_markdown(),
      legend.title = element_text(face = "bold"),
      plot.subtitle = ggtext::element_markdown()
    )
}

viz_sig_gen <- function(n, contrasts, genus_relab, metadata, ann) {
  
  gns <- contrasts[[n, "genus"]]
  term <- contrasts[[n, "term"]]
  cnt <- contrasts[[n, "contrast"]]
  
  ann <- ann %>%
    filter(genus == gns, term == !!term, contrast == cnt)
  
  if (gns == "shannon_index") {
    relab <- FALSE
    p <- "p.value"
  } else {
    relab <- TRUE
    p <- "fdr"
  }
  
  fdr <- contrasts[[n, p]] %>%
    rstatix::p_format() %>%
    sprintf("*q*-value: %s", .)
  
  .viz_sig_gen(gns, term, fdr, genus_relab, metadata, relab, ann)
}

run_mod <- function(genus, genus_relab, metadata) {
  
  print(genus)
  
  inp <- genus_relab %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    select(sample_id, all_of(genus)) %>%
    rename(relab = 2) %>%
    inner_join(metadata, by = "sample_id") %>%
    mutate_if(is.character, as.factor)
  
  min_val <- inp %>%
    filter(relab > 0) %>%
    summarise(min_val = min(relab, na.rm = TRUE)) %>%
    pull(min_val)
  
  if (!is.finite(min_val)) {
    min_val <- 1e-06
  }
  
  inp <- inp %>%
    mutate(relab = log10(relab + 0.5 * min_val))
  
  mod <- lmerTest::lmer(
    relab ~ timepoint + hygiene + pigs_per_pen + feeders_per_pen + (1 | pig_id),
    data = inp
  )
  
  # --- ANOVA ---
  anova_res <- mod %>%
    anova() %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(genus = genus)
  
  # --- PRIMARY CONTRASTS ---
  hygiene_contr <- emmeans::emmeans(mod, ~ hygiene) %>%
    emmeans::contrast(method = "pairwise") %>%
    as.data.frame() %>%
    mutate(term = "hygiene", genus = genus)
  
  pigs_contr <- emmeans::emmeans(mod, ~ pigs_per_pen) %>%
    emmeans::contrast(method = "pairwise") %>%
    as.data.frame() %>%
    mutate(term = "pigs_per_pen", genus = genus)
  
  feeders_contr <- emmeans::emmeans(mod, ~ feeders_per_pen) %>%
    emmeans::contrast(method = "pairwise") %>%
    as.data.frame() %>%
    mutate(term = "feeders_per_pen", genus = genus)
  
  contrasts_res <- bind_rows(hygiene_contr, pigs_contr, feeders_contr)
  
  # --- TIMEPOINT-SPECIFIC CONTRASTS ---
  # These models are only used to annotate the plots with within-timepoint contrasts.
  # The 'term' column is kept throughout so hygiene/feeders/pigs annotations cannot be mixed up.
  
  mod_hygiene <- lmerTest::lmer(
    relab ~ timepoint * hygiene + pigs_per_pen + feeders_per_pen + (1 | pig_id),
    data = inp
  )
  
  sec_cnt_hygiene <- emmeans::emmeans(mod_hygiene, ~ hygiene | timepoint) %>%
    emmeans::contrast("pairwise") %>%
    as_tibble() %>%
    mutate(term = "hygiene", genus = genus) %>%
    select(genus, term, everything())
  
  mod_pigs_per_pen <- lmerTest::lmer(
    relab ~ timepoint * pigs_per_pen + hygiene + feeders_per_pen + (1 | pig_id),
    data = inp
  )
  
  sec_cnt_pigs_per_pen <- emmeans::emmeans(mod_pigs_per_pen, ~ pigs_per_pen | timepoint) %>%
    emmeans::contrast("pairwise") %>%
    as_tibble() %>%
    mutate(term = "pigs_per_pen", genus = genus) %>%
    select(genus, term, everything())
  
  mod_feeders_per_pen <- lmerTest::lmer(
    relab ~ timepoint * feeders_per_pen + hygiene + pigs_per_pen + (1 | pig_id),
    data = inp
  )
  
  sec_cnt_feeders_per_pen <- emmeans::emmeans(mod_feeders_per_pen, ~ feeders_per_pen | timepoint) %>%
    emmeans::contrast("pairwise") %>%
    as_tibble() %>%
    mutate(term = "feeders_per_pen", genus = genus) %>%
    select(genus, term, everything())
  
  sec_contrasts_res <- bind_rows(
    sec_cnt_hygiene,
    sec_cnt_pigs_per_pen,
    sec_cnt_feeders_per_pen
  )
  
  list(
    anova = anova_res,
    contrasts = contrasts_res,
    sec_contrasts = sec_contrasts_res
  )
}

get_sig_gen <- function(anova_tab, contrasts_tab) {
  unique(
    c(
      filter(anova_tab, fdr <= 0.05, term != "timepoint")$genus,
      filter(contrasts_tab, fdr <= 0.05, term != "timepoint")$genus
    )
  )
}

make_heatmap <- function(genus_relab, metadata) {
  
  meta <- metadata %>%
    select(-pig_id) %>%
    rename(Timepoint = 2, Hygiene = 3, "Pigs per pen" = 4, "Feeders per pen" = 5) %>%
    column_to_rownames("sample_id")
  
  ann_col <- list(
    "Feeders per pen" = get_pal("feeders_per_pen"),
    "Pigs per pen" = get_pal("pigs_per_pen"),
    "Hygiene" = get_pal("hygiene"),
    "Timepoint" = get_pal("timepoint")
  )
  
  pheatmap::pheatmap(
    log10(t(genus_relab + 1e-06)),
    color = viridis::viridis(100, option = "D"),
    show_colnames = FALSE,
    annotation_col = meta,
    annotation_colors = ann_col,
    filename = "output/817/plots/heatmap.png",
    cellwidth = 1,
    cellheight = 15
  )
}

###########################################
## run alpha and beta diversity analysis ##
###########################################

run_alpha_diversity_analysis <- function(mat, metadata) {
  
  diversity <- vegan::diversity(mat) %>%
    as.data.frame() %>%
    rename(shannon_index = 1)
  
  div_results <- run_mod("shannon_index", diversity, metadata)
  
  cntrsts <- div_results$contrasts %>%
    filter(p.value <= 0.05)
  
  sec_con <- div_results$sec_contrasts
  
  if (nrow(cntrsts) > 0) {
    alpha_diversity_plots <- lapply(seq_len(nrow(cntrsts)), function(x) {
      viz_sig_gen(
        x,
        cntrsts,
        diversity,
        metadata,
        sec_con
      )
    }) %>%
      patchwork::wrap_plots(nrow = 2)
    
    ggsave(
      "output/817/plots/alpha_diversity.png",
      alpha_diversity_plots,
      width = 7.5,
      height = 7.5
    )
  }
}

run_beta_diversity_analysis <- function(mat, metadata) {
  
  veg_dis <- vegan::vegdist(mat, method = "bray")
  
  adonis2 <- run_adonis2(veg_dis, metadata)
  ggsave("output/817/plots/adonis2.png", adonis2$plot, width = 7.5, height = 5)
  
  pcoa <- ape::pcoa(veg_dis)
  
  pcoa_df <- pcoa$vectors[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    inner_join(metadata, by = "sample_id") %>%
    as_tibble()
  
  variables <- c("hygiene", "feeders_per_pen", "pigs_per_pen")
  
  lapply(variables, function(var) {
    p <- make_pcoa_plot(pcoa_df, pcoa, var)
    filename <- sprintf("output/817/plots/pcoa_%s.png", var)
    ggsave(filename, p, width = 7.5, height = 5)
  })
}

run_diversity_analysis <- function(mat, metadata) {
  run_alpha_diversity_analysis(mat, metadata)
  run_beta_diversity_analysis(mat, metadata)
}

##################
## get metadata ##
##################

get_metadata <- function() {
  
  cols <- c(
    "sample_id",
    "timepoint",
    "hygiene",
    "pigs_per_pen",
    "feeders_per_pen",
    "pig_id"
  )
  
  xlsx <- "data/Copy of Meta data 817.xlsx"
  
  readxl::read_excel(xlsx) %>%
    select(1, 11, 8, 9, 10, 5) %>%
    setNames(cols) %>%
    mutate(sample_id = gsub(".*\\.|.*_", "", sample_id)) %>%
    mutate(
      timepoint = ifelse(grepl("Transfer", timepoint), "D45", timepoint) %>%
        str_replace(" .*", "") %>%
        factor(levels = c("D0", "D14", "D45", "Finisher"))
    ) %>%
    mutate(hygiene = factor(hygiene, levels = c("Basic", "Optimal"))) %>%
    mutate(
      pigs_per_pen = str_to_sentence(english::english(pigs_per_pen)) %>%
        factor(levels = c("Seven", "Fourteen"))
    ) %>%
    mutate(
      feeders_per_pen = str_to_sentence(english::english(feeders_per_pen)) %>%
        factor(levels = c("One", "Two"))
    )
}

###############
## Trial 817 ##
###############

outdirs <- sprintf("output/817/%s", c("plots", "tables"))

lapply(outdirs, function(outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
})

metadata <- get_metadata()

features <- qiime2R::read_qza("data/qiime/table.qza")$data %>%
  suppressWarnings()

colnames(features) <- gsub(".*\\.|.*_", "", colnames(features))

taxonomy <- qiime2R::read_qza("data/qiime/taxonomy.qza")$data %>%
  suppressWarnings()

genus_relab <- get_genus_relab(features, taxonomy, metadata, 0.01, 0.1, clr = FALSE)

## run genus-level differential abundance analysis
results <- lapply(colnames(genus_relab), function(genus) {
  run_mod(genus, genus_relab, metadata)
})

anova_tab <- results %>%
  lapply(function(x) x$anova) %>%
  bind_rows() %>%
  group_by(term) %>%
  mutate(fdr = p.adjust(`Pr(>F)`, "BH")) %>%
  ungroup() %>%
  select(genus, everything()) %>%
  arrange(genus)

contrasts_tab <- results %>%
  lapply(function(x) x$contrasts) %>%
  bind_rows() %>%
  group_by(term, contrast) %>%
  mutate(fdr = p.adjust(p.value, "BH")) %>%
  ungroup() %>%
  select(genus, everything()) %>%
  arrange(genus)

sig_gen <- get_sig_gen(anova_tab, contrasts_tab)

sec_contrasts_tab <- results %>%
  lapply(function(x) x$sec_contrasts) %>%
  bind_rows() %>%
  group_by(term, contrast, timepoint) %>%
  mutate(fdr = p.adjust(p.value, "BH")) %>%
  mutate(fdr = rstatix::p_format(fdr)) %>%
  ungroup() %>%
  filter(genus %in% sig_gen) %>%
  select(timepoint, genus, term, contrast, fdr, everything())

# Write three worksheets to one Excel file
writexl::write_xlsx(
  list(
    "MixedEffectsModel" = anova_tab,
    "PrimaryContrasts" = contrasts_tab,
    "SecondaryContrasts" = sec_contrasts_tab
  ),
  path = "output/817/tables/model_results.xlsx"
)

sig_cntrst <- contrasts_tab %>%
  filter(fdr <= 0.05, term != "timepoint")

if (nrow(sig_cntrst) > 0) {
  
  plot_list <- lapply(seq_len(nrow(sig_cntrst)), function(n) {
    
    var <- sig_cntrst[[n, "term"]]
    
    out_dir <- sprintf("output/817/plots/sig_gen/%s", var)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    filename <- sig_cntrst[[n, "genus"]] %>%
      str_replace_all("; ", ".") %>%
      str_replace_all("/", "_") %>%
      sprintf("%s/%s.png", out_dir, .)
    
    p <- viz_sig_gen(n, sig_cntrst, genus_relab, metadata, sec_contrasts_tab)
    
    ggsave(filename, p, width = 7.5, height = 5)
    
    return(p)
  })
  
  pdf("output/817/plots/sig_genera.pdf", width = 7.5, height = 5)
  print(plot_list)
  dev.off()
}

make_heatmap(genus_relab[, sig_gen, drop = FALSE], metadata)

mat <- get_genus_relab(features, taxonomy, metadata, abundance = 0, prevalence = 0)

run_diversity_analysis(mat, metadata)
