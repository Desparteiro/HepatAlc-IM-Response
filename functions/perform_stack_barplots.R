perform_stacked_barplots <- function(dataframe,
                                    variables,
                                    grouping_var,
                                    individual,
                                    colors,
                                    title,
                                    hide_x_text,
                                    output_folder,
                                    color_palette,
                                    label_y,
                                    custom_labels = NULL,
                                    reorder_stacks = FALSE,
                                    legend_title = NULL,
                                    abundance_type = "relative",
                                    y_axis_limits = NULL,
                                    y_axis_breaks = NULL) {
  order_var_name <- paste("order_by_", grouping_var, sep = "")
  if (exists(order_var_name, envir = .GlobalEnv)) {
    order <- get(order_var_name, envir = .GlobalEnv)
  } else {
    order <- NULL  
  }
  
  selected_data <- dataframe[, c(individual, grouping_var, variables), drop = FALSE]
  
  selected_data <- selected_data[complete.cases(selected_data), ]
  if (nrow(selected_data) == 0) {
    message(paste("Warning: No valid data for", title, "- skipping plot."))
    return(NULL)
  }
  
  long_data <- pivot_longer(
    selected_data,
    cols = variables,
    names_to = "Variable",
    values_to = "Value"
  )
  long_data$Value <- as.numeric(long_data$Value)
  
  if (all(long_data$Value == 0 | is.na(long_data$Value))) {
    message(paste(
      "Warning: All values are zero or missing for",
      title,
      "- skipping plot."
    ))
    return(NULL)
  }
  
  is_already_normalized <- function(data, id_col, value_col, tol = 1e-3) {
    check <- data %>%
      group_by(!!sym(id_col)) %>%
      summarise(total = sum(!!sym(value_col), na.rm = TRUE)) %>%
      pull(total)
    all(abs(check - 100) < tol)
  }
  
  if (abundance_type == "relative") {
    needs_normalization <- !is_already_normalized(long_data, individual, "Value")
    if (needs_normalization) {
      long_data <- long_data %>%
        group_by(!!sym(individual)) %>%
        mutate(Value = (Value / sum(Value, na.rm = TRUE)) * 100)
    }
  }
  
  
  if (length(order) > 0) {
    if (length(order) == 1) {
      colors <- c("blue")
    } else if (length(order) == 2) {
      colors <- c("blue", "red")
    } else if (length(order) == 3) {
      colors <- c("darkgreen", "orange", "darkred")
    } else if (length(order) == 4) {
      colors <- c("darkgreen", "orange", "red", "darkred")
    }
  }
  
  print(paste("Performing stack barplot by", grouping_var))
  
  if (!is.null(grouping_var)) {
    dataframe[[grouping_var]] <- factor(dataframe[[grouping_var]], levels = order)
  }
  for (variable in variables) {
    dataframe[[variable]] <- as.numeric(dataframe[[variable]])
  }
  
  if (!is.null(grouping_var)) {
    dataframe <- subset(dataframe, dataframe[[grouping_var]] %in% order)
  }
  if (nrow(dataframe) == 0) {
    warning(paste("No data available for grouping variable:", grouping_var))
    return(NULL)
  }
  
  dataframe <- dataframe[, c(individual, grouping_var, variables)]
  
  individual_data <- dataframe %>%
    pivot_longer(
      cols = all_of(variables),
      names_to = "variables",
      values_to = "value"
    )
  
  if (nrow(individual_data) == 0) {
    warning(
      paste(
        "No data to plot for variable:",
        variables,
        "in grouping variable:",
        grouping_var
      )
    )
    return(NULL)
  }
  
  individual_data$variables <- factor(individual_data$variables, levels = variables)
  
  average_data <- individual_data %>%
    group_by(!!sym(grouping_var), variables) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  average_data$variables <- factor(average_data$variables, levels = variables)
  
  if (is.null(y_axis_breaks)) {
    max_value <- max(individual_data$value, na.rm = TRUE)
    step_size <- max_value / 4  
    y_axis_breaks <- seq(0, max_value, by = step_size)
  }
  

  
  plot_average <- ggplot(average_data, aes(
    x = !!sym(grouping_var),
    y = mean,
    fill = variables
  )) +
    geom_bar(stat = "identity", colour = "black") +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    font("ylab", size = 20) +
    font("xy.text", size = 20) +
    font("legend.title", size = 14) +
    font("legend.text", size = 14) +
    rotate_x_text(45) +
    labs(x = grouping_var,
         y = label_y,
         fill = legend_title %||% title) +
    scale_y_continuous(
      limits = if (!is.null(y_axis_limits))
        y_axis_limits
      else
        c(0, max_value),
      breaks = y_axis_breaks,
      expand = c(0.01, 0.01)
    ) +
    scale_fill_brewer(palette = "Set3")
  
  dir_path <- paste(output_folder, "stacked_grouped/", comparison, "/", sep = "")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  pdf(paste(dir_path, title, "_by_", grouping_var, ".pdf", sep = ""))
  grid.arrange(plot_average, ncol = 1)
  plot_name <- paste("plot_average_stacked", title, "by", grouping_var, sep = "_")
  assign(plot_name, plot_average, envir = .GlobalEnv)
  dev.off()
  
  y_label <- ifelse(abundance_type == "relative",
                    "Relative concentration (%)",
                    label_y)
  
  plot_individual <- ggplot(individual_data, aes(
    x = !!sym(individual),
    y = value,
    fill = variables
  )) +
    geom_bar(stat = "identity", colour = "black") +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      strip.text.x = element_text(size = 25),
      strip.text.y = element_text(size = 25),
      strip.background = element_rect(
        color = "white",
        fill = "white",
        size = 1.5,
        linetype = "solid"
      ),
      axis.text.x = if (hide_x_text)
        element_blank()
      else
        element_text(angle = 45, hjust = 1),
      axis.ticks.x = if (hide_x_text)
        element_blank()
      else
        element_line()
    ) +
    labs(x = grouping_var,
         y = y_label,
         fill = legend_title %||% title) +
    scale_fill_manual(values = colors,
                      labels = custom_labels %||% levels(individual_data$variables)) +
    scale_y_continuous(
      limits = if (!is.null(y_axis_limits))
        y_axis_limits
      else
        c(0, max_value),
      breaks = y_axis_breaks,
      expand = c(0.01, 0.01)
    )
  
  if (!is.null(grouping_var)) {
    plot_individual <- plot_individual +
      facet_wrap(as.formula(paste("~", grouping_var)), scales = "free_x")
  }
  
  if (is.null(color_palette)) {
    plot_individual <- plot_individual + scale_fill_brewer(palette = "Set3")
  } else {
    plot_individual <- plot_individual + scale_fill_manual(values = color_palette)
  }
  
  dir_path <- paste(output_folder,
                    "stacked_individual/",
                    grouping_var,
                    "/",
                    sep = "")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  output_file <- paste(dir_path, title, "_by_", grouping_var, ".pdf", sep = "")
  pdf(output_file)
  
  print(plot_individual)  
  dev.off()
  
  plot_name <- paste("plot_individual_stacked", title, "by", grouping_var, sep = "_")
  assign(plot_name, plot_individual, envir = .GlobalEnv)
    graphics.off()
}
