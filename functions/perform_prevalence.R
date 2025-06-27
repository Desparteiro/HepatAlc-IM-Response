perform_prevalence <- function(phyloseq_object,
                               comparison,
                               output_folder) {
  taxa_level <- str_match(phyloseq_object, "phyloseq_(\\w+)")[, 2]

  dir_path <- paste(output_folder, "prevalence/", comparison, "/", sep = "")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE, )
  

  print(paste(
    "Performing the prevalence analysis of",
    taxa_level,
    "by",
    comparison,
    sep = " "
  ))
  
  df_object <- as.data.frame(otu_table(get(phyloseq_object)))
  df_object[] <- lapply(df_object, as.numeric)

  order <- get(paste("order_by_", comparison, sep = ""))
  print(paste("Order obtained for comparison:", toString(order)))  
  
  Presence_table <- t(apply(df_object, 2, function(x) {
    ifelse(x > 0, "Presence", "Absence")
  }))
  Presence_table <- as.data.frame(Presence_table)

  grouping_metadata <- sam_data(get(phyloseq_object))[[comparison]]
  Presence_table <- cbind(Presence_table, Groupe = grouping_metadata)
  colnames(Presence_table)[ncol(Presence_table)] <- comparison

  Presence_table <- subset(Presence_table, grouping_metadata %in% order)

  prevalence_df <- Presence_table %>%
    group_by_at(vars(!!comparison)) %>%
    summarise(across(where(is.character), ~ mean(. == "Presence") * 100))

  cols_to_keep <- sapply(prevalence_df[, -1], function(x) {
    length(unique(x)) > 1
  })
  prevalence_df <- prevalence_df[, c(TRUE, cols_to_keep)]

  
  cols <- names(prevalence_df)[-1]
  pvalues <- list()
  
  for (col in cols) {
    contingency_table <- table(Presence_table[[comparison]], Presence_table[[col]])
    print(paste("Contingency table for:", col))  
    pvalues[[col]] <- chisq.test(contingency_table)$p.value
  }
  
  adjusted_pvalues <- p.adjust(unlist(pvalues), method = "BH")
  print("Adjusted p-values calculated.")  
  
  pvalue_df <- data.frame(
    "metadata" = names(adjusted_pvalues),
    "pvalue unadjusted" = unlist(pvalues),
    "pvalue" = adjusted_pvalues
  )
  
  print("P-value data frame created.")  
  
  prevalence_df <- as.data.frame(t(prevalence_df))
  colnames(prevalence_df) <- as.character(prevalence_df[1, ])
  prevalence_df <- prevalence_df[-1, ]
  Prevalence_table <- bind_cols(prevalence_df, pvalue_df)
  
  Prevalence_table <- Prevalence_table %>%
    dplyr::select("metadata", order, "pvalue.unadjusted", "pvalue")
  Prevalence_table <- Prevalence_table[order(Prevalence_table$pvalue.unadjusted), ]
  print("Prevalence table created and sorted.") 
  
  new_Prevalence_table <- Prevalence_table
  new_Prevalence_table <- new_Prevalence_table[apply(new_Prevalence_table[, -which(names(new_Prevalence_table) == "metadata")], 1, function(x) {
    sum(x > 0) > 0
  }), ]
  print("Filtered new prevalence table for non-zero rows.")  
  
  for (level in order) {
    level_count <- sum(Presence_table[[comparison]] == level)
    print(paste("Processing level:", level, "with count:", level_count)) 
    
    colnames(Prevalence_table)[colnames(Prevalence_table) == level] <-
      paste(level, " (n = ", level_count, ")", sep = "")
    
    path <- paste(dir_path,
                  taxa_level,
                  ".csv",
                  sep = "")
    print(paste("Writing prevalence table to:", path)) 
    
    tryCatch({
      write.table(Prevalence_table, path, dec = ".", sep = ";", quote = FALSE, row.names = FALSE)
    }, error = function(e) {
      message("ERROR during write.table: ", e$message)
    })
    
  }
  
  rel_ab_ps <- microbiome::transform(get(phyloseq_object), "compositional")  
  rel_ab_df <- as.data.frame(otu_table(rel_ab_ps))
  print("Relative abundance data frame created.")  
  
  threshold <- 0.00001
  Prevalence_table$pvalue <- NULL
  Prevalence_table <- Prevalence_table %>%
    mutate_at(vars(-metadata), as.numeric)
  
  single_condition_species <- rownames(Prevalence_table)[apply(Prevalence_table[, -which(names(Prevalence_table) == "metadata")], 1, function(x) {
    sum(x > 0) == 1
  })]
  
  print("Identified single condition species.")  
  
  group_levels <- unique(Presence_table[[comparison]])
  group_dfs <- list()
  for (group in group_levels) {
    print(paste("Processing group:", group))  
    group_data <- Presence_table[Presence_table[[comparison]] == group, ]
    
    group_rel_ab_df <- rel_ab_df[, Presence_table[[comparison]] == group]
    group_rel_ab_df <- group_rel_ab_df[rownames(group_rel_ab_df) %in% single_condition_species, ]
    group_rel_ab_df <- group_rel_ab_df[rowSums(group_rel_ab_df > 0) > 0, ]
    
    subset_single_condition_species <- row.names(group_rel_ab_df)
    if (!is.null(subset_single_condition_species) &&
        length(subset_single_condition_species) > 0) {
      for (species in subset_single_condition_species) {
        print(paste("Calculating for species:", species)) 
        
        non_zero_values <- group_rel_ab_df[species, ][group_rel_ab_df[species, ] != 0]
        rel_abundance <- 100 * mean(non_zero_values)
        
        X <- sum(group_data[, species] == "Presence")
        Y <- sum(group_data[, species] == "Presence") + sum(group_data[, species] == "Absence")
        prevalence_fraction <- paste(" ", X, "/", Y)
        
        prevalence_percent <- as.numeric(new_Prevalence_table[species, group])
        
        single_condition_df <- data.frame(
          "Species" = species,
          "Group" = group,
          "Relative Abundance (%)" = rel_abundance,
          "Prevalence (%)" = prevalence_percent,
          "Prevalence_Fraction" = prevalence_fraction
        )
        
        if (prevalence_percent != 0) {
          group_dfs[[group]] <- rbind(group_dfs[[group]], single_condition_df)
          print(paste("Added data for species:", species)) 
        }
      }
      if (!(is.null(group_dfs[[group]]))) {
        names(group_dfs[[group]]) <- c(
          "Taxa",
          group,
          "Relative Abundance (%)",
          "Prevalence (%)",
          "Detection frequency"
        )
        group_dfs[[group]] <- group_dfs[[group]][order(-group_dfs[[group]]$`Prevalence (%)`), ]
        path <- paste(
          "output_plots/prevalence/",
          comparison,
          "/",
          taxa_level,
          "_",
          group,
          "_only.csv",
          sep = ""
        )
        print(paste("Writing group data to:", path))  
        write.table(
          group_dfs[[group]],
          path,
          dec = ".",
          sep = ";",
          quote = FALSE,
          row.names = FALSE
        )
      }
    }
    
    plot_list <- list()
    filtered_data <- Prevalence_table[Prevalence_table$pvalue.unadjusted < 0.10, ]
    long_data <- reshape2::melt(filtered_data, id.vars = c("metadata", "pvalue.unadjusted"))
    
    long_data$variable <- sub(" \\(.*", "", long_data$variable)
    print("Long data format prepared for plotting.")  
    
    group_levels <- unique(Presence_table[[comparison]])
    if (length(group_levels) == 1) {
      colors <- c("blue")
    } else if (length(group_levels) == 2) {
      colors <- c("#000099", "#F8766D")
    } else if (length(group_levels) == 3) {
      colors <- c("darkgreen", "orange", "darkred")
    } else if (length(group_levels) == 4) {
      colors <- c("darkgreen", "orange", "red", "darkred")
    }
    
    long_data$metadata <- sub("^[a-z]__", "", long_data$metadata)
    long_data$metadata <- gsub("_", " ", long_data$metadata)
    
    for (species in unique(long_data$metadata)) {
      print(paste("Creating plot for species:", species))  
      
      species_data <- subset(long_data, metadata == species)
      species_data$variable <- factor(species_data$variable, levels = order)
      
      plot <- ggplot(species_data, aes(x = variable, y = value, fill = variable)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(
          x = "Condition",
          y = paste0(species, "\nprevalence (%)"),
          fill = "Condition"
        ) +
        scale_y_continuous(limits = c(0, 100)) +  
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_fill_manual(values = colors) +
        scale_color_manual(values = colors) +
        guides(color = "none", fill = "none") +
        theme_classic() +
        theme(axis.title.x = element_blank()) +
        font("ylab", size = 16) +
        font("xy.text", size = 16) +
        font("legend.title", size = 16) +
        font("legend.text", size = 16) +
        rotate_x_text(45)
      
      plot_list[[species]] <- plot
      
      Cairo::CairoPNG(filename = paste0(dir_path, species, ".png"))
      print(plot)
      dev.off()
      print(paste("Plot saved for species:", species))  
    }
    
    return(plot_list)
  }
}