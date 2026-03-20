get_breaks <- function(x) {
  seq(min(x), max(x), by = (max(x) - min(x)) / (length(x) - 1))
}

get_pch <- function(x) {
  lapply(x, function(x) ifelse(x < 0.05, "*", NA)) %>%
    unlist()
}

make_heatmap <- function(genus_relab, meta, trial_no) {
  
  ####################
  ## col annotation ##
  ####################
  
  ann_col <- meta %>%
    arrange(sample_id) %>%
    column_to_rownames("sample_id") %>%
    select(day, group) %>%
    mutate(day = sprintf("Day %s", str_pad(day, 2, "left", "0"))) %>%
    mutate(dss = ifelse(grepl("\\+", group), "DSS+", "DSS-")) %>%
    mutate(group = gsub(" .*", "", group)) %>%
    select(group, dss, day) %>%
    setNames(c("Group", "DSS", "Day"))
  
  ## group palette
  pal_group <- viridis::viridis(3)
  names(pal_group) <- sort(unique(ann_col$Group))
  ## day palette
  days <- sort(unique(ann_col$Day))
  pal_day <- RColorBrewer::brewer.pal(length(days), "Reds")
  names(pal_day) <- days
  ## dss palette
  pal_dss <- c("grey", "blue")
  names(pal_dss) <- sort(unique(ann_col$DSS))
  
  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = ann_col$Group,
    DSS = ann_col$DSS,
    Day = ann_col$Day,
    col = list(Group = pal_group, DSS = pal_dss, Day = pal_day),
    show_legend = FALSE
  )
  
  ####################
  ## row annotation ##
  ####################
  
  top_tax <- genus_relab %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(!sample_id) %>%
    group_by(name) %>%
    summarise(value = mean(value)) %>%
    ungroup() %>%
    top_n(25, value) %>%
    pull(name)
  
  ann_row <- read_delim(sprintf("output/tables/lmer_results.%s.tsv", trial_no)) %>%
    pivot_wider(names_from = "covariate", values_from = "p") %>%
    filter(taxon %in% top_tax) %>%
    column_to_rownames("taxon")
  
  col <- circlize::colorRamp2(
    breaks = get_breaks(ann_row),
    colors = viridis::viridis(length(ann_row), direction = -1, option = "magma")
  )
  
  right_ann <- ComplexHeatmap::HeatmapAnnotation(
    `Genus ~ Group` = ComplexHeatmap::anno_simple(
      ann_row[,1],
      pch = get_pch(ann_row[,1]),
      col = col
    ),
    `Genus ~ Day` = ComplexHeatmap::anno_simple(
      ann_row[,2],
      pch = get_pch(ann_row[,2]),
      col = col
    ),
    `Genus ~ Group:Day` = ComplexHeatmap::anno_simple(
      ann_row[,3],
      pch = get_pch(ann_row[,3]),
      col = col
    ),
    which = "row",
    show_legend = FALSE
  )
  
  #############
  ## legends ##
  #############
  
  legend_0 <- ComplexHeatmap::Legend(
    title = "log(abundance [%])",
    at = 1:4,
    col_fun = circlize::colorRamp2(
      breaks = 1:4,
      colors = viridis::viridis(4, option = "turbo")
    )
  )
  
  legend_1 <- ComplexHeatmap::Legend(
    labels = names(pal_group),
    title = "Group",
    legend_gp = grid::gpar(fill = as.character(pal_group))
  )
  
  legend_2 <- ComplexHeatmap::Legend(
    labels = names(pal_dss),
    title = "DSS",
    legend_gp = grid::gpar(fill = as.character(pal_dss))
  )
  
  legend_3 <- ComplexHeatmap::Legend(
    labels = names(pal_day),
    title = "Day",
    legend_gp = grid::gpar(fill = as.character(pal_day))
  )
  
  legend_4 <- ComplexHeatmap::Legend(
    title = "p-value",
    at = c(0.001, 0.01, 0.1, 1),
    col_fun = col
  )
  
  legend <- ComplexHeatmap::packLegend(legend_0,
                                       legend_1,
                                       legend_2,
                                       legend_3,
                                       legend_4)
  
  ######################
  ## make the heatmap ##
  ######################
  
  mat <- as.matrix(log1p(t(genus_relab))[top_tax, ])
  
  hmap <- ComplexHeatmap::Heatmap(
    mat,
    name = "log(abundance [%])",
    rect_gp = grid::gpar(col = "black", lwd = 0.5),
    top_annotation = top_ann,
    right_annotation = right_ann,
    show_column_names = FALSE,
    col = viridis::viridis(100, option = "turbo"),
    show_heatmap_legend = FALSE
  )
  
  filename <- sprintf("output/plots/heatmap.%s.png", trial_no)
  png(filename, width = 17.5, height = 8.125, units = "in", res = 300)
  ComplexHeatmap::draw(
    hmap,
    annotation_legend_list = list(legend),
    annotation_legend_side = "left",
    padding = unit(c(2.5, 2.5, 2.5, 50), "mm")
  )
  dev.off()
  
}
