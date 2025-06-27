perform_target_comparison <- function(phyloseq_object, target, secondary_target = NULL, output_folder, formatted_target, grouping_var = NULL) {
  current_object <- phyloseq_object
  
  if (!is.null(grouping_var)) {
    non_na_samples <- !is.na(sample_data(current_object)[[grouping_var]])
    current_object <- prune_samples(non_na_samples, current_object)
  }
  
  rel_current_object <- microbiome::transform(current_object, "compositional")
  
  taxa_names_list <- taxa_names(rel_current_object)
  target_taxa <- taxa_names_list[taxa_names_list == target]
  secondary_target_taxa <- NULL
  
  if (!is.null(secondary_target)) {
    secondary_target_taxa <- taxa_names_list[taxa_names_list == secondary_target]
    if (length(secondary_target_taxa) == 0) {
      stop(paste("Secondary target taxa", secondary_target, "not found in the taxa names"))
    }
  }
  
  if (length(target_taxa) == 0) {
    stop(paste("Target taxa", target, "not found in the taxa names"))
  }
  
  rel_current_object <- prune_taxa(c(target_taxa, secondary_target_taxa), rel_current_object)
  otu_table_df <- t(as.data.frame(otu_table(rel_current_object)))
  
  sam_data <- data.frame(sample_data(rel_current_object))
  dataframe <- cbind(otu_table_df, sam_data)
  
  if (target %in% colnames(dataframe)) {
    dataframe <- dataframe %>%
      mutate(across(all_of(target), ~ . * 100))
  } else {
    stop(paste("Target", target, "not found in dataframe columns"))
  }
  
  if (!is.null(secondary_target)) {
    if (secondary_target %in% colnames(dataframe)) {
      dataframe <- dataframe %>%
        mutate(across(all_of(secondary_target), ~ . * 100))
    } else {
      stop(paste("Secondary target", secondary_target, "not found in dataframe columns"))
    }
  }
  
  if (is.null(grouping_var)) {
    stop("Grouping variable must be provided and cannot be NULL.")
    
  } else if (is.numeric(dataframe[[grouping_var]])) {
    if (!is.null(secondary_target)) {
      plot <- perform_correlation(
        dataframe,
        VARIABLE_X = secondary_target,
        VARIABLE_Y = target,
        BREAKS_X = NULL,
        BREAKS_Y = NULL,
        LIMITS_X = NULL,
        LIMITS_Y = NULL,
        EXPAND_X = NULL,
        EXPAND_Y = NULL,
        LABEL_X = secondary_target,
        LABEL_Y = formatted_target
      )
    } else {
      stop("For correlation, both target and secondary target must be provided when grouping_var is numeric.")
    }
  } else if (is.factor(dataframe[[grouping_var]]) || is.character(dataframe[[grouping_var]])) {
    plot <- perform_mean_comparison(
      dataframe,
      response_var = target,
      grouping_var = grouping_var,
      output_folder = output_folder,
      name = formatted_target,
      order = get(paste0("order_by_", grouping_var)),
      breaks = NULL,
      limits = NULL,
      expand = NULL,
      ypositions = NULL,
      paired = FALSE,
      add_p_value = TRUE
    )
  } else {
    stop("Grouping variable must be either categorical (factor/character) or numeric.")
  }
  
  return(plot)
}
