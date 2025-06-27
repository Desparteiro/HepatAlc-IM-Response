perform_beta_diversity <- function(physeq_object, beta_metric, grouping_var, title, output_folder, distance_matrix = NULL) {
  
  order <- get(paste("order_by_", grouping_var, sep = ""))
  if (length(order) > 0) {
    num_groups <- length(order)
    if (num_groups == 1) {
      colors <- c("blue")
    } else if (num_groups == 2) {
      colors <- c("blue", "red")
    } else if (num_groups == 3) {
      colors <- c("darkgreen", "orange", "darkred")
    } else if (num_groups == 4) {
      colors <- c("darkgreen", "orange", "red", "darkred")
    } else {
      colors <- RColorBrewer::brewer.pal(num_groups, "Set1")
    }
  }
  
  physeq_object@sam_data[[grouping_var]] <- factor(physeq_object@sam_data[[grouping_var]], levels = order)
  
  non_na_samples <- !is.na(sample_data(physeq_object)[[grouping_var]])
  physeq_object <- prune_samples(non_na_samples, physeq_object)
  
  df_physeq <- data.frame(physeq_object@sam_data)
  df_physeq[[grouping_var]] <- factor(df_physeq[[grouping_var]], levels = order)
  df_physeq <- subset(df_physeq, df_physeq[[grouping_var]] %in% order)
  
  if (is.null(distance_matrix)) {
    if (beta_metric == "aitchison") {
      otu_matrix <- as.matrix(physeq_object@otu_table + 1)
      clr_matrix <- compositions::clr(otu_matrix)
      dist_matrix <- vegan::vegdist(clr_matrix, method = "euclidean")

      } else {
      dist_matrix <- phyloseq::distance(physeq_object, method = beta_metric)
    }
  } else {
    dist_matrix <- distance_matrix
  }

  formula <- as.formula(paste("dist_matrix ~", paste(grouping_var, collapse = " + ")))
  set.seed(123)
  
  permanova_result <- adonis2(formula, data = df_physeq, permutations = 9999, na.action = na.exclude)
  
  r_squared <- round(permanova_result$R2[1], 2)
  pvalue <- round(permanova_result$`Pr(>F)`[1], 3)
  
  if (pvalue < 0.001) {
    text_annotation <- paste0("PERMANOVA, R² = ", r_squared, ", p < ", 0.001)
  } else {
    text_annotation <- paste0("PERMANOVA, R² = ", r_squared, ", p = ", pvalue)
  }
  
  ordination <- ordinate(physeq_object, method = "PCoA", distance = dist_matrix)
  
  plot <- plot_ordination(physeq_object, ordination, color = grouping_var) +
    stat_ellipse(geom = "polygon", alpha = 0.03, linetype = 1) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.position = "bottom"    )+
  
    geom_point(size = 1)
  
  plot_build <- ggplot_build(plot)
  x_limits <- plot_build$layout$panel_scales_x[[1]]$range$range
  x_center <- mean(x_limits)
  
  plot <- plot + 
    annotate("text", x = x_center, y = Inf, hjust = 0.5, vjust = 1, 
             label = text_annotation, size = 6, fontface = "bold")
  
  permanova_table <- as.data.frame(permanova_result)
  permanova_table <- tableGrob(permanova_table, theme = my_theme)
  
  dir_path <- paste(output_folder, "beta_metrics/", comparison, "/", sep = "")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  pdf(file = paste0(dir_path, beta_metric, "_by_", grouping_var, ".pdf"))
  grid.arrange(permanova_table, plot, ncol = 1)
  dev.off()
  
  plot_name <- paste("plot", beta_metric, "by", grouping_var, sep = "_")
  assign(plot_name, plot, envir = .GlobalEnv)
  graphics.off()
}
