perform_correlation <- function(DATAFRAME,
                                VARIABLE_X,
                                VARIABLE_Y,
                                output_folder,
                                BREAKS_X = NULL,
                                BREAKS_Y = NULL,
                                LIMITS_X = NULL,
                                LIMITS_Y = NULL,
                                EXPAND_X = NULL,
                                EXPAND_Y = NULL,
                                LABEL_X = NULL,
                                LABEL_Y = NULL) {
  print(paste(
    "Performing correlation of ",
    VARIABLE_Y,
    " by ",
    VARIABLE_X,
    sep = ""
  ))
  DATAFRAME[[VARIABLE_X]] <- as.numeric(DATAFRAME[[VARIABLE_X]])
  DATAFRAME[[VARIABLE_Y]] <- as.numeric(DATAFRAME[[VARIABLE_Y]])
  DATAFRAME <- DATAFRAME[complete.cases(DATAFRAME[[VARIABLE_Y]], DATAFRAME[[VARIABLE_X]]), ]
  
  p_plot <- ggplot(DATAFRAME, aes_string(x = VARIABLE_X, y = VARIABLE_Y)) +
    geom_point() +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "black",
      linetype = "solid",
      size = 1.5
    ) +
    labs(x = VARIABLE_X, y = VARIABLE_Y) +
    scale_x_continuous(breaks = BREAKS_X, limits = LIMITS_X) +
    scale_y_continuous(breaks = BREAKS_Y, limits = LIMITS_Y) +
    theme_classic() +
    theme(text = element_text(size = 20))
  
  
  custom_round <- function(x) {
    if (x == 0) {
      return(0)
    }
    
    magnitude <- floor(log10(abs(x)))
    
    factor <- 10 ^ (-magnitude + 2)  
    
  
    x_scaled <- x * factor
    
    remainder <- x_scaled %% 5
    
    if (remainder < 2.5) {
      rounded_value <- (x_scaled - remainder) / factor
    } else {
      rounded_value <- (x_scaled + 5 - remainder) / factor
    }
    
    return(rounded_value)
  }
  
  if (is.null(BREAKS_Y)) {
    min_y <- custom_round(1.1 * min(DATAFRAME[[VARIABLE_Y]]))
    max_y <- custom_round(1.1 * max(DATAFRAME[[VARIABLE_Y]]))
    BREAKS_Y <- seq(min_y, max_y, by = 0.25 * (max_y - min_y))
  }
  if (is.null(LIMITS_Y)) {
    LIMITS_Y <- c(min_y, max_y)
  }
  if (is.null(EXPAND_Y)) {
    EXPAND_Y <- c(0, 0)
  }
  
  p_plot <- p_plot + scale_y_continuous(breaks = BREAKS_Y,
                                        limits = LIMITS_Y,
                                        expand = EXPAND_Y)
  
  if (is.null(BREAKS_X)) {
    min_X <- custom_round(1.1 * min(DATAFRAME[[VARIABLE_X]]))
    max_X <- custom_round(1.1 * max(DATAFRAME[[VARIABLE_X]]))
    BREAKS_X <- seq(min_X, max_X, by = 0.25 * (max_X - min_X))
    
  }
  if (is.null(LIMITS_X)) {
    LIMITS_X <- c(min_X, max_X)
  }
  if (is.null(EXPAND_X)) {
    EXPAND_X <- c(0, 0)
  }
  
  p_plot <- p_plot + scale_x_continuous(breaks = BREAKS_X,
                                        limits = LIMITS_X,
                                        expand = EXPAND_X)
  if (!is.null(LABEL_X)) {
    p_plot <-  p_plot + xlab(LABEL_X)
  }
  
  if (!is.null(LABEL_Y)) {
    p_plot <-  p_plot + ylab(LABEL_Y)
  }
  shapiro_x <- shapiro_test(DATAFRAME[[VARIABLE_X]])
  shapiro_y <- shapiro_test(DATAFRAME[[VARIABLE_Y]])
  
  if (any(c(shapiro_x$p.value, shapiro_y$p.value) < 0.05)) {
    corr_method <- "spearman"
  } else {
    corr_method <- "pearson"
  }
  corr_label <- toTitleCase(corr_method)
  corr_test <- cor.test(DATAFRAME[[VARIABLE_X]], DATAFRAME[[VARIABLE_Y]], method = corr_method)
  
  if (corr_method == "spearman") {
    label <- "rho"
  } else {
    label <- "r"
  }
  Value <- paste(
    label,
    " = ",
    round(corr_test$estimate, 2),
    ", p =",
    formatC(corr_test$p.value, format = "e", digits = 2)
  )
  
  p_plot <-  p_plot + annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1,
    vjust = 1,
    label = corr_label,
    size = 6,
    fontface = "bold"
  ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1,
      vjust = 2.5,
      label = Value,
      size = 6,
      fontface = "bold"
    )
  dir_path <- paste(output_folder, "Correlation/",
                    VARIABLE_Y,
                    "_by_",
                    VARIABLE_X,
                    sep = "")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  pdf(paste(dir_path, "/", VARIABLE_Y, "_by_", VARIABLE_X, ".pdf", sep = ""))
  print(p_plot)
  dev.off()
  graphics.off()
  plot_name <- paste("plot", VARIABLE_Y, "by", VARIABLE_X, sep = "_")
  assign(plot_name, p_plot, envir = .GlobalEnv)
  return(p_plot)
}