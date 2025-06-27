perform_full_stats <- function(dataset, 
                               grouping_variable, 
                               grouping_variable_levels,
                               columns_to_remove, 
                               rounding_number, 
                               ID, 
                               pairing, 
                               sample_name) { 
  dataset <- arrange(dataset, ID, grouping_variable) 
  dataset <- dataset %>% mutate_if(is.integer, as.numeric)
  for (column in columns_to_remove) {
    dataset[[column]] <- NULL
  }
  
  dataset[[grouping_variable]] <- factor(dataset[[grouping_variable]], levels = grouping_variable_levels)
  dataset <- subset(dataset, dataset[[grouping_variable]] %in% grouping_variable_levels)
  numeric_columns <- sapply(dataset, is.numeric)
  
  round_safe <- function(x, digits) {
    tryCatch({
      ifelse(is.na(x) | !is.numeric(x), NA, format(round(x, digits), nsmall = digits))
    }, error = function(e) NA)
  }
  
  clean_string <- function(x) {
    gsub("\\s*-\\s*", "-", gsub("\\s+", " ", gsub("\\( ", "(", gsub(" \\)", ")", x))))
  }
  

  numeric_columns <- sapply(dataset, is.numeric)
  stratified_stats_dataset <- dataset %>%
    group_by_at(vars(grouping_variable)) %>%
    summarise(across(where(is.numeric), list(
      mean = ~ round_safe(mean(., na.rm = TRUE), rounding_number),
      sd = ~ round_safe(sd(., na.rm = TRUE), rounding_number),
      median = ~ round_safe(median(., na.rm = TRUE), rounding_number),
      q1 = ~ round_safe(quantile(., 0.25, na.rm = TRUE), rounding_number),
      q3 = ~ round_safe(quantile(., 0.75, na.rm = TRUE), rounding_number)
    ), .names = "{.col}_{.fn}"))
  
  stratified_stats_dataset <- as.data.frame(t(stratified_stats_dataset))
  colnames(stratified_stats_dataset) <- levels(dataset[[grouping_variable]])
  stratified_stats_dataset <- stratified_stats_dataset[-1, ]
  
  numeric_column_names <- names(numeric_columns[numeric_columns == TRUE])
  empty_data <- data.frame(matrix(
    ncol = length(levels(dataset[[grouping_variable]])),
    nrow = length(numeric_column_names),
    dimnames = list(numeric_column_names, levels(dataset[[grouping_variable]]))
  ))
  
  new_stratified_stats_dataset <- dplyr::bind_rows(stratified_stats_dataset, empty_data)
  new_stratified_stats_dataset <- as.data.frame(as.matrix(new_stratified_stats_dataset))
  
  for (i in seq_len(nrow(new_stratified_stats_dataset))) {
    if (is.na(new_stratified_stats_dataset[i, 1])) {
      for (j in seq_len(ncol(new_stratified_stats_dataset))) {
        variable_name <- row.names(new_stratified_stats_dataset)[i]
        if (variable_name %in% colnames(dataset)) {
          mean_value <- stratified_stats_dataset[paste0(variable_name, "_mean"), j]
          sd_value <- stratified_stats_dataset[paste0(variable_name, "_sd"), j]
          median_value <- stratified_stats_dataset[paste0(variable_name, "_median"), j]
          q1_value <- stratified_stats_dataset[paste0(variable_name, "_q1"), j]
          q3_value <- stratified_stats_dataset[paste0(variable_name, "_q3"), j]
          
          N <- dataset %>%
            filter(!!sym(grouping_variable) == levels(dataset[[grouping_variable]])[j]) %>%
            summarise(N = sum(!is.na(!!sym(variable_name)))) %>%
            pull(N)
          
           new_stratified_stats_dataset[i, j] <- clean_string(paste0(
             mean_value, " ± ", sd_value, " [n = ", N, "] / ",
             median_value, " (", q1_value, "-", q3_value, ") [n = ", N, "]"
          ))
        }
      }
    }
  }
  
  summary_stats_dataset <- dataset %>%
    summarise(across(where(is.numeric), list(
      mean = ~ round(mean(., na.rm = TRUE), rounding_number),
      sd = ~ round(sd(., na.rm = TRUE), rounding_number),
      median = ~ round(median(., na.rm = TRUE), rounding_number),
      q1 = ~ round(quantile(., 0.25, na.rm = TRUE), rounding_number),
      q3 = ~ round(quantile(., 0.75, na.rm = TRUE), rounding_number),
      n = ~ sum(!is.na(.))
    ), .names = "{.col}_{.fn}"))
  
  summary_stats_dataset <- as.data.frame(t(summary_stats_dataset))
  
  new_summary_stats_dataset <- data.frame(Statistic = character(nrow(summary_stats_dataset) / 6), stringsAsFactors = FALSE)
  
  output_index <- 1
  
  mean_value <- NULL
  sd_value <- NULL
  median_value <- NULL
  q1_value <- NULL
  q3_value <- NULL
  N <- NULL
  
  for (i in seq_len(nrow(summary_stats_dataset))) {
    stat_type <- rownames(summary_stats_dataset)[i]
    
    value <- summary_stats_dataset[i, "V1"] 
    
    if (grepl("_mean$", stat_type)) {
      mean_value <- value
    } else if (grepl("_sd$", stat_type)) {
      sd_value <- value
    } else if (grepl("_median$", stat_type)) {
      median_value <- value
    } else if (grepl("_q1$", stat_type)) {
      q1_value <- value
    } else if (grepl("_q3$", stat_type)) {
      q3_value <- value
    } else if (grepl("_n$", stat_type)) {
      N <- value
      
      new_summary_stats_dataset[output_index, "Statistic"] <- paste0(
        mean_value, " ± ", sd_value, " [n = ", N, "] / ", 
        median_value, " (", q1_value, "-", q3_value, ") [n = ", N, "]"
      )
      output_index <- output_index + 1  
    }
  }
  
  print(paste("Descriptive statistics calculated for: ", grouping_variable))  
  
  new_stratified_stats_dataset <- cbind(new_summary_stats_dataset, new_stratified_stats_dataset)
  
  full_stats <- new_stratified_stats_dataset[-grep("_mean|_sd|_median|_q1|_q3", rownames(new_stratified_stats_dataset)), ]
  
  for (var in numeric_column_names) {
    loop_dataset <- subset(dataset, dataset[[grouping_variable]] %in% grouping_variable_levels)
    
    loop_dataset <- subset(loop_dataset, !is.na(loop_dataset[[var]]))
    loop_dataset[[grouping_variable]] <- droplevels(loop_dataset[[grouping_variable]])
    
    formula <- as.formula(paste(var, "~", grouping_variable))
    for (level in levels(loop_dataset[[grouping_variable]])) {
      shapiro_results <- NULL
      column_name <- paste("Shapiro_p_value", level, sep = "_")
      shapiro_results <- tryCatch(
        {
          byf.shapiro(formula, loop_dataset)
        },
        error = function(e) {
          message(
            "Failure of Shapiro for the level ",
            level,
            " of ",
            grouping_variable,
            " of ",
            var,
            " : ",
            e
          )
          NULL
        }
      )
      if (is.null(shapiro_results)) {
        full_stats[var, column_name] <- NA
      } else {
        full_stats[var, column_name] <-
          as.numeric(shapiro_results$tab[level, "p-value"])
      }
    }
    print(paste("Shapiro is calculated for: ", grouping_variable))
    
    levene_results <- tryCatch(
      {
        car::leveneTest(formula, loop_dataset)
      },
      error = function(e) {
        message(
          "Failure of Levene's test for the variable ",
          var,
          " : ",
          e
        )
        NULL
      }
    )
    if (is.null(levene_results) || is.nan(levene_results[1, "Pr(>F)"])) {
      full_stats[var, "Levene_p_value"] <- NA
    } else {
      full_stats[var, "Levene_p_value"] <-
        as.numeric(levene_results[1, "Pr(>F)"])
    }
    print(paste("Levene is calculated for: ", grouping_variable)) 
    
    
    if (length(levels(loop_dataset[[grouping_variable]])) == 2) {
      if (any(full_stats[var, grep("Shapiro_p_value", names(full_stats))] < 0.05) |
          any(is.na(full_stats[var, grep("Shapiro_p_value", names(full_stats))]))) {
        if (pairing == TRUE) {
          na_rows <- is.na(loop_dataset[[var]])
          IDs <- loop_dataset[[ID]][na_rows]
          loop_dataset <- loop_dataset[!loop_dataset[[ID]] %in% IDs, ]
          loop_dataset <- loop_dataset[order(loop_dataset[[ID]], loop_dataset[[grouping_variable]]), ]
          wilcox_results <- tryCatch(
            {
              wilcox.test(formula, data = loop_dataset, paired = TRUE)
            },
            error = function(e) {
              message(
                "Failure of paired Wilcoxon for the comparison by ",
                grouping_variable,
                " of ",
                var,
                " : ",
                e
              )
              NULL
            }
          )
          if (is.null(wilcox_results)) {
            full_stats[var, "paired_Wilcoxon_p_value"] <- NA
          } else {
            full_stats[var, "paired_Wilcoxon_p_value"] <-
              as.numeric(wilcox_results$p.value)
          }
        } else {
          wilcox_results <- tryCatch(
            {
              
              wilcox.test(formula, data = loop_dataset, exact = TRUE)
            },
            error = function(e) {
              message(
                "Failure of unpaired Wilcoxon for the comparison by ",
                grouping_variable,
                " of ",
                var,
                " : ",
                e
              )
              NULL
            }
          )
          if (is.null(wilcox_results) || is.na(wilcox_results$p.value))  {
            full_stats[var, "unpaired_Wilcoxon_p_value"] <- NA
          } else {
            full_stats[var, "unpaired_Wilcoxon_p_value"] <-
              as.numeric(wilcox_results$p.value)
          }
        }
      } else {
        if (pairing == TRUE &
            all(!is.na(full_stats[var, grep("Shapiro_p_value", names(full_stats))]))) {
          na_rows <- is.na(loop_dataset[[var]])
          IDs <- loop_dataset[[ID]][na_rows]
          loop_dataset <- loop_dataset[!loop_dataset[[ID]] %in% IDs, ]
          loop_dataset <- loop_dataset[order(loop_dataset[[ID]], loop_dataset[[grouping_variable]], decreasing = TRUE), ]
          
          t_results <- tryCatch(
            {
              t.test(
                formula,
                data = loop_dataset,
                paired = TRUE,
                var.equal = TRUE
              )
            },
            error = function(e) {
              message(
                "Failure of paired t-test for the comparison by ",
                grouping_variable,
                " of ",
                var,
                " : ",
                e
              )
              NULL
            }
          )
          if (is.null(t_results)) {
            full_stats[var, "paired_t_test_p_value"] <- NA
          } else {
            full_stats[var, "paired_t_test_p_value"] <-
              as.numeric(t_results$p.value)
          }
        }
        if (pairing == FALSE) {
          if (full_stats[var, "Levene_p_value"] < 0.05 &
              all(full_stats[var, grep("Shapiro_p_value", names(full_stats))] > 0.05)) {
            
            t_results <- tryCatch(
              {
                t.test(
                  formula,
                  data = loop_dataset,
                  var.equal = FALSE
                )
              },
              error = function(e) {
                message(
                  "Failure of Welch t-test for the comparison by ",
                  grouping_variable,
                  " of ",
                  var,
                  " : ",
                  e
                )
                NULL
              }
            )
            if (is.null(t_results)) {
              full_stats[var, "Welch_t_test_p_value"] <- NA
            } else {
              full_stats[var, "Welch_t_test_p_value"] <- as.numeric(t_results$p.value)
            }
          }
          
          if (full_stats[var, "Levene_p_value"] > 0.05 &
              all(full_stats[var, grep("Shapiro_p_value", names(full_stats))] > 0.05)) {
            
            t_results <- tryCatch(
              {
                t.test(
                  formula,
                  data = loop_dataset,
                  var.equal = TRUE
                )
              },
              error = function(e) {
                message(
                  "Failure of unpaired t-test for the comparison by ",
                  grouping_variable,
                  " of ",
                  var,
                  " : ",
                  e
                )
                NULL
              }
            )
            if (is.null(t_results)) {
              full_stats[var, "unpaired_t_test_p_value"] <- NA
            } else {
              full_stats[var, "unpaired_t_test_p_value"] <- as.numeric(t_results$p.value)
            }
          }
        }
      }
    } else {
      if (any(full_stats[var, grep("Shapiro_p_value", names(full_stats))] < 0.05) |
          any(is.na(full_stats[var, grep("Shapiro_p_value", names(full_stats))]))) {
        adhoc_result <- tryCatch(
          {
            kruskal_test(data = dataset, as.formula(paste("", var, " ~ ", grouping_variable, sep = "")))
          },
          error = function(e) {
            message(
              "Failure of Kruskal's test for the variable ",
              var,
              " : ",
              e
            )
            NULL
          }
        )
        if (is.null(adhoc_result) || is.nan(adhoc_result$p)) {
          full_stats[var, "kruskal test p-value"] <- NA
        } else {
          full_stats[var, "kruskal test p-value"] <- as.numeric(adhoc_result$p)
          full_stats[var, "kruskal test statistic"] <- as.numeric(adhoc_result$statistic)
        }
        
        posthoc_result <- tryCatch(
          {
            dunn_test(data = dataset, as.formula(paste("", var, " ~ ", grouping_variable, sep = "")), p.adjust.method = "BH")
          },
          error = function(e) {
            NULL
          }
        )
        if (!is.null(posthoc_result)) {
          p_values <- posthoc_result$p.adj
          p_values <- sapply(p_values, as.numeric)
          groups <- paste(posthoc_result$group1, posthoc_result$group2, sep = "-")
          new_data <- data.frame(group = groups, p_value = p_values)
          formatted_strings <- paste(groups, ": p =", p_values)
          final_string <- paste(formatted_strings, collapse = " ; ")
          full_stats[var, "Dunn's p-values"] <- final_string
          
          
          cld_results <- rcompanion::cldList(
            formula = p_value ~ group,
            data = new_data,
            threshold = 0.05,
            remove.space = FALSE,
            remove.zero = FALSE
          )
          
          for (level in cld_results$Group) {
            letter <- cld_results[cld_results$Group == level, "MonoLetter"]
            full_stats[var, level] <- paste(full_stats[var, level], letter, sep = "")
          }
        }
      } else {
        if (full_stats[var, "Levene_p_value"] < 0.05) {
          adhoc_result <- oneway.test(data = dataset, as.formula(paste("", var, " ~ ", grouping_variable, sep = "")), )
          if (is.null(adhoc_result)) {
            full_stats[var, "One-way p-value"] <- NA
          } else {
            full_stats[var, "One-way p-value"] <- as.numeric(adhoc_result$p.value)
            full_stats[var, "One-way F-value"] <- as.numeric(adhoc_result$statistic)
          }
        } else {
          adhoc_result <- aov(data = dataset, as.formula(paste("", var, " ~ ", grouping_variable, sep = "")), )
          if (is.null(adhoc_result)) {
            full_stats[var, "ANOVA p-value"] <- NA
          } else {
            full_stats[var, "ANOVA p-value"] <- as.numeric(summary(adhoc_result)[[1]][grouping_variable, "Pr(>F)"])
            full_stats[var, "ANOVA F-value"] <- summary(adhoc_result)[[1]][grouping_variable, "F value"]
          }
        }
        
        print(paste("Adhoc testing is calculated for: ", grouping_variable))  
        
        posthoc_result <- dataset %>% tukey_hsd(as.formula(paste("", var, " ~ ", grouping_variable, sep = "")), p.adjust.method = "BH")
        if (!is.null(posthoc_result)) {
          p_values <- posthoc_result$p.adj
          p_values <- sapply(p_values, as.numeric)
          groups <- paste(posthoc_result$group1, posthoc_result$group2, sep = "-")
          new_data <- data.frame(group = groups, p_value = p_values)
          formatted_strings <- paste(groups, ": p =", p_values)
          final_string <- paste(formatted_strings, collapse = " ; ")
          full_stats[var, "Tukey's p-values"] <- final_string
          cld_results <- rcompanion::cldList(
            formula = p_value ~ group,
            data = new_data,
            threshold = 0.05,
            remove.space = FALSE,
            remove.zero = FALSE
          )
          
          for (level in cld_results$Group) {
            letter <- cld_results[cld_results$Group == level, "MonoLetter"]
            full_stats[var, level] <- paste(full_stats[var, level], letter, sep = "")
          }
        }
      }
    }
  }
  
  print(paste("Posthoc testing is calculated for: ", grouping_variable))   
  
  dataset_excluded <- dataset[, !(names(dataset) %in% c(ID, sample_name, "sample","Sample"))]
  character_columns <- sapply(dataset_excluded, is.character)
  
  
  for (var in names(character_columns[character_columns == TRUE])) {
    count <- table(dataset_excluded[[var]], exclude = NA)
    for (i in 1:length(count)) {
      frequency <- count[i] / sum(count, na.rm = TRUE)
      
      formatted_level <- paste0(count[i], " (", round(frequency * 100, rounding_number), "%)")
      full_stats[paste(var, " (", names(count[i]), ")", sep = ""), "Total"] <- formatted_level
    }
    
    levels_variable_A <- levels(dataset_excluded[[grouping_variable]])
    for (level in levels_variable_A) {
      subset_data <- dataset_excluded[dataset_excluded[[grouping_variable]] == level, ]
      
      count <- table(subset_data[[var]], exclude = NA)
      for (i in 1:length(count)) {
        frequency <- count[i] / sum(count, na.rm = TRUE)
        formatted_level <- paste0(count[i], " (", round(frequency * 100, rounding_number), "%)")
        full_stats[paste(var, " (", names(count[i]), ")", sep = ""), level] <- formatted_level
      }
    }
    contingency_table <- table(dataset_excluded[[var]], dataset_excluded[[grouping_variable]])
    
    if (nrow(contingency_table) > 1 & ncol(contingency_table) > 1) {
      fisher_p_value <- tryCatch(
        {
          fisher.test(contingency_table, simulate.p.value = TRUE, B = 1e5)$p.value
        },
        error = function(e) {
          message("Fisher test failed for ", var, " : ", e$message)
          NA
        }
      )
      
      full_stats[paste(var, "test"), "Fisher_p_value"] <- as.numeric(fisher_p_value)
      
      chi_square_p_value <- tryCatch(
        {
          chisq.test(contingency_table)$p.value
        },
        error = function(e) {
          message("Chi-square test failed for ", var, " : ", e$message)
          NA
        }
      )
      
      full_stats[paste(var, "test"), "Chi_square_p_value"] <- as.numeric(chi_square_p_value)
      
    } else {
      full_stats[paste(var, "test"), "Fisher_p_value"] <- NA
      full_stats[paste(var, "test"), "Chi_square_p_value"] <- NA
    }
  }
  
  
  for (level in levels(dataset_excluded[[grouping_variable]])) {
    string <- paste(level, " : n = ", sum(dataset[[grouping_variable]] %in% level), sep = "")
    colname_to_change <- colnames(full_stats)[colnames(full_stats) == level]
    if (length(colname_to_change) > 0) {
      colnames(full_stats)[colnames(full_stats) == colname_to_change] <- string
    }
  }
  
  colnames(full_stats)[1] <- paste("Total : n = ", nrow(dataset), sep = "")

  print(paste("Formation of the output dataset for :", grouping_variable)) 
  return(full_stats) 
}


