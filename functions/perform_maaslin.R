perform_maaslin <-  function(phyloseq_object,
                             comparison,
                             metadata,
                             fixed_effects,
                             random_effects,
                             references,
                             output_folder) {
  df_object <- as.data.frame(otu_table(get(phyloseq_object)))
  df_object[] <- lapply(df_object, as.numeric)
  df_object <-
    df_object[!grepl("GGB|SGB|UNKNOWN", row.names(df_object), perl = TRUE), ]
  grouping_var <- comparison
  order <- get(paste("order_by_", grouping_var, sep = ""))
  metadata_maaslin <- subset(metadata, metadata[[grouping_var]] %in% order)
  metadata_maaslin[[grouping_var]] <- factor(metadata_maaslin[[grouping_var]], levels = order)
  taxa_level <- str_match(phyloseq_object, "phyloseq_(\\w+)")[, 2]
  
  fixed_collapsed <- paste(fixed_effects, collapse = "_")
  random_collapsed <- paste(random_effects, collapse = "_")
  
  dir_path <- paste(output_folder,
                    "maaslin2/",
                    fixed_collapsed,
                    random_collapsed,
                    "/",
                    sep = "")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  Sys.sleep(0.5)
  
  output_results <- paste(dir_path, taxa_level, sep = "")
  if (grepl("clr", phyloseq_object)) {
    norm_method <- "NONE"
    transform_method <- "NONE"
  } else {
    norm_method <- "TSS"
    transform_method <- "LOG"
  }
  Sys.sleep(0.5)
  
  maaslin <- Maaslin2(
    input_data = df_object,
    input_metadata = metadata_maaslin,
    output = output_results,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    min_prevalence = 0.1,
    normalization = norm_method,
    transform = transform_method,
    correction = "BH",
    max_significance = 0.05,
    reference = references
  )
  maaslin_result <-
    paste(output_results, "/significant_results.tsv", sep = "")
}