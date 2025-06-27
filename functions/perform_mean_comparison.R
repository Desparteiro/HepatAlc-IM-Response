perform_mean_comparison <- function(dataframe,
                                    name,
                                    response_var,
                                    grouping_var,
                                    order = NULL,
                                    colors = NULL,
                                    breaks = NULL,
                                    limits = NULL,
                                    expand = NULL,
                                    ypositions = NULL,
                                    paired = NULL,
                                    label_y = NULL,
                                    add_p_value = NULL,
                                    output_folder = NULL) {

    print(paste(
    "Performing the comparison of",
    response_var,
    "by",
    grouping_var,
    sep = " "
  ))
  
  format_p_value <- function(p_value, digits = 4) {
    if (p_value < 0.001) {
      return(sprintf("%.1e", p_value)) 
    } else {
      return(sprintf(paste0("%.", digits, "f"), p_value)) 
    }
  }
  
  variability <- dataframe %>%
    group_by_at(grouping_var) %>%
    summarise(sd = sd(!!sym(response_var), na.rm = TRUE)) %>%
    ungroup()
  
  zero_variability <- all(!is.na(variability$sd) &
                            variability$sd == 0)
  
  if (zero_variability) {
    print("Zero variability detected across all groups. Skipping further analysis.")
    return(NULL) 
  }
  
  
  plot_list <- list()
  if (!is.factor(dataframe[[grouping_var]])) {
    dataframe[[grouping_var]] <- as.factor(dataframe[[grouping_var]])
  }
  
  if (is.null(order) || length(order) == 0) {
    order <- unique(na.omit(dataframe[[grouping_var]]))
  }
  
  
  if (length(order) > 0) {
    dataframe[[grouping_var]] <- factor(dataframe[[grouping_var]], levels = order)
  }
  
  dataframe <- dataframe[!is.na(dataframe[[grouping_var]]), ]
  
  dataframe[[response_var]] <- as.numeric(as.character(dataframe[[response_var]]))
  
  dataframe <- dataframe[!is.na(dataframe[[response_var]]), ]
  
  dataframe[[grouping_var]] <- droplevels(dataframe[[grouping_var]])
  
  

  if (!(is.null(dataframe[[grouping_var]])) &
      nlevels(dataframe[[grouping_var]]) > 0) {
    p_plot <- ggboxplot(
      dataframe,
      y = response_var,
      x = grouping_var,
      color = grouping_var,
      fill = grouping_var,
      xlab = grouping_var,
      ylab = label_y,
      add = "jitter",
      bxp.errorbar = TRUE,
      outlier.shape = 19,
      bxp.errorbar.width = 0.15
    ) +
      guides(color = "none", fill = "none") +
      theme_classic() +
      theme(axis.title.x = element_blank()) +
      font("ylab", size = 20) +
      font("xy.text", size = 20) +
      font("legend.title", size = 16) +
      font("legend.text", size = 16) +
      rotate_x_text(45)
    
    
    if (is.null(colors)) {
      if (length(levels(dataframe[[grouping_var]])) == 1) {
        colors <- c("blue")
      } else if (length(levels(dataframe[[grouping_var]])) == 2) {
        colors <- c("#000099", "#F8766D")
      } else if (length(levels(dataframe[[grouping_var]])) == 3) {
        colors <- c("darkgreen", "orange", "darkred")
      } else if (length(levels(dataframe[[grouping_var]])) == 4) {
        colors <- c("darkgreen", "orange", "red", "darkred")
      }
    }
    
    p_plot <- p_plot +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = yarrr::transparent(colors, 0.8))
    
    if (length(levels(dataframe[[grouping_var]])) > 0) {
      summary <- dataframe %>%
        group_by_at(grouping_var) %>%
        summarise(
          mean = round(mean(!!sym(response_var)), 1),
          median = round(median(!!sym(response_var)), 1),
          sd = round(sd(!!sym(response_var)), 1),
          min = round(min(!!sym(response_var)), 1),
          max = round(max(!!sym(response_var)), 1)
        )
      summary_table <- tableGrob(as.data.frame(summary), theme = my_theme)
    } else {
      summary <- NULL
      summary_table <- NULL
    }
    
    num_levels <- length(unique(dataframe[[grouping_var]]))
    if (paired == TRUE) {
      individuals_to_remove <- unique(dataframe$Individual[is.na(dataframe[[response_var]])])
      dataframe <- dataframe[!dataframe$Individual %in% individuals_to_remove, ]
    }
    adhoc_success <- FALSE
    posthoc_success <- FALSE
    two_sample_success <- FALSE
    two_sample_test <- NULL
    adhoc_result <- NULL
    posthoc_result <- NULL
    posthoc_result_filtered <- NULL
    formula <- as.formula(paste("`", response_var, "` ~ ", grouping_var, sep = ""))
    
    shapiro <- tryCatch({
      byf.shapiro(formula, data = dataframe)
    }, error = function(e) {
      message("Could not calculate Shapiro test: ", e)
      NULL
    })
    normality <- !(is.null(shapiro) ||
                     any(shapiro$tab[2] < 0.05)) 
    
    levene <- tryCatch({
      levene_test(data = dataframe, formula = dataframe[[response_var]] ~ dataframe[[grouping_var]])
    }, error = function(e) {
      message("Could not calculate Levene test: ", e)
      NULL
    })
    homoscedasticity <- !(is.null(levene) ||
                            levene$p < 0.05) 
    if (num_levels < 2) {
      test <- NULL
    } else if (num_levels == 2) {
      if (normality == TRUE) {
        if (paired == TRUE) {
          test <- "Paired Student"
        } else if (homoscedasticity == TRUE) {
          test <- "Unpaired Student"
        } else {
          test <- "Welch Student"
        }
      } else {
        if (paired == TRUE) {
          test <- "Paired Wilcoxon"
        } else {
          test <- "Unpaired Wilcoxon"
        }
      }
    } else if (normality == TRUE) {
      if (homoscedasticity == TRUE) {
        test <- "ANOVA"
        posthoc_test <- "Tukey HSD"
      } else {
        test <- "Welch ANOVA"
        posthoc_test <- "Dunn's test"
      }
    } else {
      test <- "Kruskal-Wallis"
      posthoc_test <- "Dunn's test"
    }
    if (add_p_value == FALSE) {
      test <- NULL
    }
    if (is.null(test)) {
      print(
        paste(
          "No testing can be performed for the analysis of",
          response_var,
          "by",
          grouping_var,
          sep = " "
        )
      )
    } else {
      print(
        paste(
          test,
          "is the appropriate hypothesis testing for the analysis of",
          response_var,
          "by",
          grouping_var,
          sep = " "
        )
      )
    }
    
    if (!is.null(test) && test == "Paired Student") {
      two_sample_success <- tryCatch({
        two_sample_test <- t_test(data = dataframe,
                                  formula = as.formula(paste(
                                    "`", response_var, "` ~ ", grouping_var, sep = ""
                                  )),
                                  paired = TRUE)
      }, error = function(e) {
        message("Could not perform Paired Student's t-test: ", e)
        NULL
      })
      if (!is.null(two_sample_test)) {
        two_sample_test$p <- formatC(two_sample_test$p,
                                     format = "f",
                                     digits = 3)
      }
    }
    
    if (!is.null(test) && test == "Unpaired Student") {
      two_sample_success <- tryCatch({
        two_sample_test <- t_test(
          data = dataframe,
          formula = as.formula(paste(
            "`", response_var, "` ~ ", grouping_var, sep = ""
          )),
          paired = FALSE,
          var.equal = TRUE
        )
      }, error = function(e) {
        message("Could not perform Unpaired Student's t-test: ", e)
        NULL
      })
      if (!is.null(test) && !is.null(two_sample_test)) {
        two_sample_test$p <- formatC(two_sample_test$p,
                                     format = "f",
                                     digits = 3)
      }
    }
    
    if (!is.null(test) && test == "Welch Student") {
      two_sample_success <- tryCatch({
        two_sample_test <- t_test(
          data = dataframe,
          formula = as.formula(paste(
            "`", response_var, "` ~ ", grouping_var, sep = ""
          )),
          paired = FALSE,
          var.equal = FALSE
        )
      }, error = function(e) {
        message("Could not perform Welch Student's t-test: ", e)
        NULL
      })
      if (!is.null(two_sample_test)) {
        two_sample_test$p <- formatC(two_sample_test$p,
                                     format = "f",
                                     digits = 3)
      }
    }
    
    if (!is.null(test) && test == "Paired Wilcoxon") {
      two_sample_success <- tryCatch({
        two_sample_test <- wilcox_test(data = dataframe,
                                       formula = as.formula(paste(
                                         "`", response_var, "` ~ ", grouping_var, sep = ""
                                       )),
                                       paired = TRUE)
      }, error = function(e) {
        message("Could not perform Paired Wilcoxon: ", e)
        NULL
      })
      
      if (!is.null(two_sample_test)) {
        two_sample_test$p <- formatC(two_sample_test$p,
                                     format = "f",
                                     digits = 3)
      }
    }
    if (!is.null(test) && test == "Unpaired Wilcoxon") {
      two_sample_success <- tryCatch({
        two_sample_test <- wilcox_test(data = dataframe,
                                       formula = as.formula(paste(
                                         "`", response_var, "` ~ ", grouping_var, sep = ""
                                       )),
                                       paired = FALSE)
      }, error = function(e) {
        message("Could not perform Unpaired Wilcoxon: ", e)
        NULL
      })
      if (!is.null(two_sample_test)) {
        two_sample_test$p <- formatC(two_sample_test$p,
                                     format = "f",
                                     digits = 3)
      }
    }
    
    
    if (!is.null(test) && test == "ANOVA") {
      adhoc_success <- tryCatch({
        adhoc_result <- aov(data = dataframe, as.formula(paste(
          "`", response_var, "` ~ ", grouping_var, sep = ""
        )))
      }, error = function(e) {
        message("Could not perform ANOVA: ", e)
        NULL
      })
      
      posthoc_success <- tryCatch({
        posthoc_result <- dataframe %>% tukey_hsd(as.formula(paste(
          "`", response_var, "` ~ ", grouping_var, sep = ""
        )), p.adjust.method = "BH")
      }, error = function(e) {
        message("Could not perform Tukey's test: ", e)
        NULL
      })
    }
    
    if (!is.null(test) && test == "Welch ANOVA") {
      adhoc_success <- tryCatch({
        adhoc_result <- oneway.test(data = dataframe, as.formula(
          paste("`", response_var, "` ~ ", grouping_variable, sep = "")
        ))
      }, error = function(e) {
        message("Could not perform Welch ANOVA: ", e)
        NULL
      })
      posthoc_success <- tryCatch({
        posthoc_result <- dunn_test(data = dataframe,
                                    as.formula(paste(
                                      "`", response_var, "` ~ ", grouping_var, sep = ""
                                    )),
                                    p.adjust.method = "BH")
        
        posthoc_result
      }, error = function(e) {
        message("Could not perform Dunn's test: ", e$message)
        
        NULL
      })
      
      if (!is.null(posthoc_success)) {
        p_values <- posthoc_success$p.value
        
        if (length(p_values) > 0) {
          dataframe <- add_column(dataframe, p_values = p_values)
        } else {
          warning("No p-values were generated by Dunn's test.")
        }
      } else {
        warning("Dunn's test failed, skipping posthoc analysis.")
      }
      
    }
    
    if (!is.null(test) && test == "Kruskal-Wallis") {
      adhoc_success <- tryCatch({
        adhoc_result <- kruskal_test(data = dataframe, as.formula(paste(
          "`", response_var, "` ~ ", grouping_var, sep = ""
        )))
      }, error = function(e) {
        message("Could not perform Kruskal-Wallis test: ", e)
        NULL
      })
      posthoc_success <- tryCatch({
        posthoc_result <- dunn_test(data = dataframe,
                                    as.formula(paste(
                                      "`", response_var, "` ~ ", grouping_var, sep = ""
                                    )),
                                    p.adjust.method = "BH")
        
        if (is.null(posthoc_result) ||
            nrow(posthoc_result) == 0) {
          stop("Dunn's test returned no results or failed.")
        }
        
        posthoc_result
      }, error = function(e) {
        message("Could not perform Dunn's test: ", e$message)
        NULL
      })
      
      if (!is.null(posthoc_success)) {
        p_values <- posthoc_success$p.value
        
        if (length(p_values) == nrow(dataframe)) {
          dataframe <- add_column(dataframe, p_values = p_values)
        } else {
          warning(
            "Mismatch between p-values length and dataframe rows or no p-values were generated."
          )
        }
      } else {
        warning("Dunn's test failed, skipping posthoc analysis.")
        dataframe <- add_column(dataframe, p_values = rep(NA, nrow(dataframe)))
      }
      
    }
     if (!(is.null(two_sample_test))) {
      if (two_sample_test$p <= 1) {
        two_sample_test <- two_sample_test %>% add_xy_position(x = grouping_var,
                                                               step.increase = 0.1,
                                                               fun = "max")
        
        two_sample_test$p_signif <- ifelse(
          as.numeric(two_sample_test$p) < 0.0001,
          "****",
          ifelse(
            as.numeric(two_sample_test$p) < 0.001,
            "***",
            ifelse(
              as.numeric(two_sample_test$p) < 0.01,
              "**",
              ifelse(
                as.numeric(two_sample_test$p) < 0.05,
                "*",
                as.numeric(two_sample_test$p)
              )
            )
          )
        )
        if (!(is.null(two_sample_test))) {
          if (two_sample_test$p <= 1) {
            two_sample_test <- two_sample_test %>% add_xy_position(x = grouping_var,
                                                                   step.increase = 0.1,
                                                                   fun = "max")
            
            threshold <- 1e-4
            
            two_sample_test$p_signif <- ifelse(
              as.numeric(two_sample_test$p) < 0.0001,
              "****",
              ifelse(
                as.numeric(two_sample_test$p) < 0.001,
                "***",
                ifelse(
                  as.numeric(two_sample_test$p) < 0.01,
                  "**",
                  ifelse(
                    as.numeric(two_sample_test$p) < 0.05,
                    "*",
                    ifelse(
                      as.numeric(two_sample_test$p) < threshold,
                      format(
                        two_sample_test$p,
                        scientific = TRUE,
                        digits = 2
                      ),
                      round(as.numeric(two_sample_test$p), 3)
                    )
                  )
                )
              )
            )
            
            if (!is.null(ypositions)) {
              p_plot <- p_plot + stat_pvalue_manual(
                tip.length = 0.01,
                two_sample_test,
                label = "p",
                size = 6,
                y.position = ypositions
              )
            } else {
              p_plot <- p_plot + stat_pvalue_manual(
                tip.length = 0.01,
                two_sample_test,
                label = "p",
                size = 6
              )
            }
          }
        }
      }
    }

    if (!(is.null(posthoc_result))) {
      posthoc_result_filtered <- subset(posthoc_result, p.adj <= 0.1)
      posthoc_result_filtered <- posthoc_result_filtered %>% add_xy_position(x = grouping_var,
                                                                             step.increase = 0.1,
                                                                             fun = "max")
      
      posthoc_result_filtered$p.adj.signif <- sapply(posthoc_result_filtered$p.adj, format_p_value)
      
      if (is.null(ypositions)) {
        p_plot <- p_plot + stat_pvalue_manual(
          posthoc_result_filtered,
          label = "p.adj.signif",
          size = 4,
          tip.length = 0.01
        )
      } else {
        p_plot <- p_plot + stat_pvalue_manual(
          posthoc_result_filtered,
          label = "p.adj.signif",
          size = 4,
          y.position = ypositions,
          tip.length = 0.01
        )
      }
    }
    
    
    custom_round <- function(x) {
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

    max_y <- max(dataframe[[response_var]], na.rm = TRUE)
    if (!is.null(posthoc_result_filtered) &&
        nrow(posthoc_result_filtered) > 0) {
      max_y <- max(max_y,
                   max(posthoc_result_filtered$y.position, na.rm = TRUE))
    }
    if (!is.null(two_sample_test) && two_sample_test$p[1] < 0.1) {
      max_y <- max(max_y, max(two_sample_test$y.position, na.rm = TRUE))
    }
    
    if (is.null(breaks)) {
      breaks <- seq(0, 1.1 * max_y, length.out = 5)
    }
    if (is.null(limits)) {
      limits <- c(0, 1.1 * max_y)
    }
    if (is.null(expand)) {
      expand <- c(0.01, 0)
    }
    
    p_plot <- p_plot + scale_y_continuous(breaks = breaks,
                                          limits = limits,
                                          expand = expand)
    
    

    threshold <- 1e-4
    
    if (!(is.null(adhoc_result))) {
      if (test == "ANOVA") {
        text_annotation_1 <- "ANOVA"
        
        p_value <- summary(adhoc_result)[[1]][grouping_var, "Pr(>F)"]
        if (p_value < threshold) {
          p_value_formatted <- format(p_value, scientific = TRUE, digits = 2)
        } else {
          p_value_formatted <- round(p_value, 3)
        }
        
        text_annotation_2 <- paste0("F = ",
                                    round(summary(adhoc_result)[[1]][grouping_var, "F value"], 1),
                                    ", ",
                                    "p = ",
                                    p_value_formatted)
      }
      
      if (test == "Welch ANOVA") {
        text_annotation_1 <- "Welch ANOVA"
        
        p_value <- adhoc_result$p.value[1]
        if (p_value < threshold) {
          p_value_formatted <- format(p_value, scientific = TRUE, digits = 2)
        } else {
          p_value_formatted <- round(p_value, 3)
        }
        
        text_annotation_2 <- paste0("F = ",
                                    round(adhoc_result$statistic[1], 1),
                                    ", ",
                                    "p = ",
                                    p_value_formatted)
      }
      
      if (test == "Kruskal-Wallis") {
        text_annotation_1 <- "Kruskal-Wallis"
        
        p_value <- adhoc_result$p  
        
        if (is.na(p_value) || is.nan(p_value)) {
          p_value_formatted <- "NaN"
        } else if (p_value < threshold) {
          p_value_formatted <- format(p_value, scientific = TRUE, digits = 2)
        } else {
          p_value_formatted <- round(p_value, 3)
        }
        
        text_annotation_2 <- paste("H = ",
                                   round(adhoc_result$statistic, 1),
                                   ", ",
                                   "p = ",
                                   p_value_formatted)
      }
      
      
      
      num_categories <- length(levels(dataframe[[grouping_var]]))
      
      combined_annotation <- paste0(text_annotation_1, ", ", text_annotation_2)

      p_plot <- p_plot +
        annotate(
          "text",
          x = ((1:num_categories) %>% mean()),
          y = Inf,
          hjust = 0.5,
          vjust = 1,
          label = combined_annotation,
          size = 4,
          fontface = "bold"
        )
    }
    
    print(
      paste(
        "Performed with success the comparison of",
        response_var,
        "by",
        grouping_var,
        sep = " "
      )
    )
    dir_path <- paste(output_folder, name, "/", grouping_var, "/", sep = "")
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    pdf(paste(dir_path, response_var, "_by_", grouping_var, ".pdf", sep = ""))
    grid.arrange(p_plot)
    if (!(is.null(summary_table))) {
      grid.arrange(summary_table)
    }
    
    
    if (!is.null(shapiro)) {
      shapiro_table <- tableGrob(as.data.frame(shapiro$tab), theme = my_theme)
      grid.arrange(shapiro_table)
      grid.text("Shapiro normality test", x = 0.5, y = 0.9)
    }
    
    if (!is.null(levene)) {
      levene_table <- tableGrob((data.frame(
        W = levene$statistic[1], `p-value` = levene$p[1]
      )), theme = my_theme)
      grid.arrange(levene_table)
      grid.text("levene homoscedasticity test",
                x = 0.5,
                y = 0.9)
    }
    
    if (!is.null(two_sample_test) && two_sample_test$p <= 1) {
      two_sample_table <- tableGrob((
        data.frame(
          group1 = two_sample_test$group1,
          group2 = two_sample_test$group2,
          n1 = two_sample_test$n1,
          n2 = two_sample_test$n2,
          statistic = two_sample_test$statistic,
          p = two_sample_test$p
        )
      ), theme = my_theme)
      grid.arrange(two_sample_table)
      grid.text(test, x = 0.5, y = 0.9)
    }
    
    if (!is.null(adhoc_result)) {
      if (test == "Kruskal-Wallis") {
        adhoc_table <- tableGrob((as.data.frame(adhoc_result)), theme = my_theme)
      } else {
        adhoc_table <- tableGrob(as.data.frame(tidy(adhoc_result)), theme = my_theme)
      }
      grid.arrange(adhoc_table)
      grid.text(test, x = 0.5, y = 0.9)
    }
    
    if (!is.null(posthoc_result)) {
      posthoc_table <- tableGrob((as.data.frame(posthoc_result)), theme = my_theme)
      grid.arrange(posthoc_table)
      grid.text(posthoc_test, x = 0.5, y = 0.9)
    }
    plot_name <- paste("plot", response_var, "by", grouping_var, sep = "_")
    assign(plot_name, p_plot, envir = .GlobalEnv)
    dev.off()
    graphics.off()
    return(p_plot)
  } else {
    print(
      paste(
        "Could not perform the comparison of",
        response_var,
        "by",
        grouping_var,
        sep = " "
      )
    )
  }
}
