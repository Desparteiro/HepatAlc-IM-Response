perform_metaphlan_to_phyloseq <- function(mpa,
                                          metadata = NULL,
                                          version = 4,
                                          verbose = TRUE,
                                          tax_lvl = "Species") {
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  
  mpa <- as.data.frame(mpa)
  abundance_cols <- 3:ncol(mpa)
  for (i in abundance_cols) {
    mpa[, i] <- as.numeric(as.character(mpa[, i]))
  }
  tax_levels <- c("Kingdom",
                  "Phylum",
                  "Class",
                  "Order",
                  "Family",
                  "Genus",
                  "Species",
                  "SGB")
  tax_lvl_int <- match(tax_lvl, tax_levels)
  
  if (is.na(tax_lvl_int))
    stop("Invalid taxonomic level. Must be one of: ",
         paste(tax_levels, collapse = ", "))
  
  tax_lengths <- sapply(strsplit(as.character(mpa$clade_name), "|", fixed = TRUE), length)
  
  if (tax_lvl == "Species") {
    data_rows <- which(tax_lengths >= tax_lvl_int)
  } else {
    data_rows <- which(tax_lengths == tax_lvl_int + 1)
  }
  
  otu_cleaned <- mpa[data_rows, 3:ncol(mpa)]
  clade_names <- mpa[data_rows, 1]
  
  if (verbose) {
    cat("Found",
        length(data_rows),
        "taxa at or below",
        tax_lvl,
        "level\n")
  }
  
  if (!is.null(metadata)) {
    inters_names <- intersect(colnames(otu_cleaned), rownames(metadata))
    if (verbose) {
      lost_mpa <- setdiff(colnames(otu_cleaned), inters_names)
      lost_meta <- setdiff(rownames(metadata), inters_names)
      if (length(lost_mpa) > 0)
        cat("MetaPhlAn samples lost:",
            paste(lost_mpa, collapse = " "),
            "\n")
      if (length(lost_meta) > 0)
        cat("Metadata samples lost:",
            paste(lost_meta, collapse = " "),
            "\n")
    }
    otu_cleaned <- otu_cleaned[, inters_names, drop = FALSE]
  }
  
  if (tax_lvl == "Species") {
    taxonomy_tab <- clade_names %>%
      as.data.frame() %>%
      set_names("clade_name") %>%
      separate(
        col = "clade_name",
        into = tax_levels[1:tax_lvl_int],
        sep = "\\|",
        fill = "right",
        remove = TRUE,
        extra = "drop"
      ) %>%
      as.matrix()
  } else {
    taxonomy_tab <- clade_names %>%
      as.data.frame() %>%
      set_names("clade_name") %>%
      separate(
        col = "clade_name",
        into = tax_levels[1:(tax_lvl_int + 1)],
        sep = "\\|",
        fill = "right",
        remove = TRUE
      ) %>%
      as.matrix()
  }
  
  tax_column <- taxonomy_tab[, tax_lvl_int]
  
  unclassified_idx <- which(
    is.na(tax_column) |
      tax_column == "" |
      grepl("unclassified", tax_column, ignore.case = TRUE) |
      grepl("unknown", tax_column, ignore.case = TRUE)
  )
  
  keep_idx <- setdiff(seq_len(nrow(taxonomy_tab)), unclassified_idx)
  
  otu_known <- otu_cleaned[keep_idx, , drop = FALSE]
  tax_known <- taxonomy_tab[keep_idx, , drop = FALSE]
  
  otu_unknown <- otu_cleaned[unclassified_idx, , drop = FALSE]
  
  if (nrow(otu_unknown) > 0) {
    merged_unknown <- as.data.frame(matrix(colSums(otu_unknown), nrow = 1))
    colnames(merged_unknown) <- colnames(otu_cleaned)
    
    unknown_label <- paste0(tolower(substr(tax_lvl, 1, 1)), "__Unknown")
    
    taxonomy_unknown <- rep("Unknown", tax_lvl_int)
    taxonomy_unknown[tax_lvl_int] <- unknown_label
    taxonomy_unknown <- matrix(taxonomy_unknown, nrow = 1)
    colnames(taxonomy_unknown) <- tax_levels[1:tax_lvl_int]
    
    rownames(merged_unknown) <- unknown_label
    
    otu_final <- rbind(otu_known, merged_unknown)
    
    if (tax_lvl == "Species") {
      taxonomy_final <- rbind(tax_known, taxonomy_unknown)
    } else {
      tax_known <- tax_known[, -tax_lvl_int - 1]
      taxonomy_final <- rbind(tax_known, taxonomy_unknown)
    }
    
  } else {
    otu_final <- otu_known
    if (tax_lvl == "Species") {
      taxonomy_final <- tax_known
    } else {
      taxonomy_final <- tax_known[, -(tax_lvl_int + 1)]
    }
  }
  
  duplicate_taxa <- unique(taxonomy_final[, tax_lvl_int])[sapply(unique(taxonomy_final[, tax_lvl_int]), function(x) {
    sum(taxonomy_final[, tax_lvl_int] == x) > 1
  })]
  
  if (length(duplicate_taxa) > 0 && verbose) {
    cat("Found",
        length(duplicate_taxa),
        "duplicate taxa that will be merged:\n")
    cat(paste(duplicate_taxa, collapse = ", "), "\n")
  }
  
  if (length(duplicate_taxa) > 0) {
    new_taxonomy <- matrix(nrow = 0, ncol = ncol(taxonomy_final))
    colnames(new_taxonomy) <- colnames(taxonomy_final)
    new_otu <- matrix(nrow = 0, ncol = ncol(otu_final))
    colnames(new_otu) <- colnames(otu_final)
    
    processed_rows <- c()
    
    for (taxon in unique(taxonomy_final[, tax_lvl_int])) {
      taxon_rows <- which(taxonomy_final[, tax_lvl_int] == taxon)
      
      if (length(taxon_rows) > 1) {
        merged_tax <- taxonomy_final[taxon_rows[1], , drop = FALSE]
        
        
        merged_otu <- matrix(colSums(otu_final[taxon_rows, , drop = FALSE]), nrow = 1)
        colnames(merged_otu) <- colnames(otu_final)
        
        rownames(merged_tax) <- taxon
        rownames(merged_otu) <- taxon
        
      } else {
        merged_tax <- taxonomy_final[taxon_rows, , drop = FALSE]
        merged_otu <- otu_final[taxon_rows, , drop = FALSE]
        
        rownames(merged_tax) <- taxon
        rownames(merged_otu) <- taxon
      }
      
      new_taxonomy <- rbind(new_taxonomy, merged_tax)
      new_otu <- rbind(new_otu, merged_otu)
      
      processed_rows <- c(processed_rows, taxon_rows)
    }
    
    taxonomy_final <- new_taxonomy
    otu_final <- new_otu
    
    if (verbose) {
      cat("After merging duplicates:",
          nrow(taxonomy_final),
          "unique taxa remain\n")
    }
  }
  
  remaining_duplicates <- unique(taxonomy_final[, tax_lvl_int])[sapply(unique(taxonomy_final[, tax_lvl_int]), function(x) {
    sum(taxonomy_final[, tax_lvl_int] == x) > 1
  })]
  
  if (length(remaining_duplicates) > 0) {
    warning("Some duplicates still remain: ",
            paste(remaining_duplicates, collapse = ", "))
  }
  
  if (nrow(otu_final) > 0) {
    profiles <- phyloseq(
      otu_table(as.matrix(otu_final), taxa_are_rows = TRUE),
      tax_table(taxonomy_final),
      errorIfNULL = FALSE
    )
    
    if (verbose) {
      cat(
        "Created phyloseq object with",
        ntaxa(profiles),
        "taxa and",
        nsamples(profiles),
        "samples\n"
      )
    }
    
    return(profiles)
    
  } else {
    stop("No valid taxa found after filtering")
  }
}