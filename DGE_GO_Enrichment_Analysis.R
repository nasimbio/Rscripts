
# Required Libraries
library(Seurat)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)

perform_complete_analysis <- function(seurat_obj, treatment, timepoint1, timepoint2, output_dir, logfc_threshold = 0.5, pval_threshold = 0.05, top_n_genes = 50) {
  message("Starting the analysis pipeline...")
  
  # Step 1: Subset Seurat object
  subset_obj <- subset(seurat_obj, subset = Treatment == treatment) 
  message("Subset dimensions: ", paste(dim(subset_obj), collapse = " x "))
  message("Seurat object subset completed.")
  subset_obj@meta.data$Timepoint <- factor(subset_obj@meta.data$Timepoint, levels = c(timepoint1, timepoint2))
  
  # Step 2: Differential Gene Expression Analysis
  dge_results_list <- list()
  for (annotation in unique(subset_obj@meta.data$Res_0.4_Annotation)) {
    message(paste("Processing annotation:", annotation))
    annotation_subset <- subset(subset_obj, subset = Res_0.4_Annotation == annotation)
    Idents(annotation_subset) <- annotation_subset@meta.data$Timepoint
    
    # Check if there is enough number of cells
    cells_timepoint1 <- subset(annotation_subset, subset = Timepoint == timepoint1)
    cells_timepoint2 <- subset(annotation_subset, subset = Timepoint == timepoint2)
    
    # Skip analysis if any treatment group has fewer than 3 cells
    if (ncol(cells_timepoint1) < 3 | ncol(cells_timepoint2) < 3) {
      message(paste("Skipping annotation", annotation, "because one of the groups has fewer than 3 cells"))
      next
    }
    
    # Perform differential expression analysis
    dge_results <- FindMarkers(annotation_subset, ident.1 = timepoint1, ident.2 = timepoint2, test.use = "wilcox")
    dge_results_list[[annotation]] <- dge_results
  }
  message("DGE analysis completed.")
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Step 3: Save DGE results
  dge_dir <- file.path(output_dir, "DEGs")
  dir.create(dge_dir, recursive = TRUE, showWarnings = FALSE)
  for (annotation in names(dge_results_list)) {
    file_path <- file.path(dge_dir, paste0(annotation, ".csv"))
    write.csv(dge_results_list[[annotation]], file_path, row.names = TRUE)
  }
  message("DGE results saved.")
  
  # Step 4: Create and save volcano plots
  volcano_dir <- file.path(output_dir, "Volcano_Plots")
  dir.create(volcano_dir, showWarnings = FALSE)
  for (annotation in names(dge_results_list)) {
    dge_data <- dge_results_list[[annotation]]
    
    # Check if dge_data is empty
    if (length(dge_data) == 0) {
      message(paste("No DEGs found for annotation:", annotation, "- skipping volcano plot."))
      next
    }
    
    # Ensure p_val and avg_log2FC are numeric
    dge_data$avg_log2FC <- as.numeric(as.character(dge_data$avg_log2FC))
    dge_data$p_val <- as.numeric(as.character(dge_data$p_val))
    dge_data$p_val_adj <- as.numeric(as.character(dge_data$p_val_adj))
    dge_data$gene <- rownames(dge_data)
    
    # Create and save volcano plot
    volcano_file <- file.path(volcano_dir, paste0(annotation, "_volcano_plot.pdf"))
    EnhancedVolcano(
      dge_data,
      lab = dge_data$gene,
      x = 'avg_log2FC',
      y = 'p_val',
      pCutoff = pval_threshold,
      FCcutoff = logfc_threshold,
      cutoffLineType = 'twodash',
      title = paste(annotation, "(", timepoint1, "vs.", timepoint2, ")"),
      cutoffLineWidth = 0.8,
      pointSize = 4.0,
      labSize = 6.0,
      colAlpha = 1,
      legendLabels = c('Not sig.', 'Log (base 2) FC', 'p-value', 'p-value & Log(base 2) FC'),
      legendPosition = 'right',
      legendLabSize = 16,
      legendIconSize = 5.0
    )
    ggsave(filename = volcano_file, device = "pdf", width = 15, height = 10)
  }
  message("Volcano plots created and saved.")
  
  # Step 5: Generate Heatmap for Subset of Significant DEGs Based on -log10(p-value)
  heatmap_dir <- file.path(output_dir, "Heatmap_Plots")
  dir.create(heatmap_dir, showWarnings = FALSE)
  for (annotation in names(dge_results_list)) {
    dge_data <- dge_results_list[[annotation]]
    if (nrow(dge_data) > 0) {
      dge_data$neg_log10_pval <- -log10(dge_data$p_val)
      sig_genes <- rownames(dge_data[dge_data$neg_log10_pval > -log10(pval_threshold) & abs(dge_data$avg_log2FC) > logfc_threshold, ])
      if (length(sig_genes) > 0) {
        top_genes <- sig_genes[1:min(top_n_genes, length(sig_genes))]  # Limit to top N genes
        message("Generating heatmap for significant DE genes...")
        valid_genes <- top_genes[top_genes %in% rownames(GetAssayData(annotation_subset, slot = "data"))]
        
        # Generate Heatmap
        heatmap_file <- file.path(heatmap_dir, paste0(annotation, "_heatmap_plot.png"))
        if (length(valid_genes) > 0) {
          heatmap_plot <- DoHeatmap(annotation_subset, features = valid_genes, group.by = "Timepoint") +
            ggtitle("Significant DE Genes Heatmap")
          ggsave(filename = heatmap_file, plot = heatmap_plot, device = "png", width = 14, height = 10, dpi = 500)
          message("Heatmap saved to: ", heatmap_file)
        } else {
          message("No valid significant genes found in scaled data for heatmap.")
        }
      } else {
        message("No significant DEGs identified. Skipping heatmap generation.")
      }
    } else {
      message("No DEGs identified. Skipping heatmap generation.")
    }
  }
  
  message("Analysis pipeline completed.")
}





# call the function
perform_complete_analysis(
  seurat_obj = No_TCR_Subclusters_filtered, 
  treatment = "B+T2", 
  timepoint1="Day21",
  timepoint2="Day14",
  output_dir = "/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Inter_Timepoints/BT2/Day21_vs_Day14",
  logfc_threshold = 0.5, 
  pval_threshold = 0.05,
  top_n_genes = 50
)


#GO
perform_go_enrichment <- function(dge_dir, output_dir, logfc_threshold = 0.5, pval_threshold = 0.05) {
  message("Starting GO enrichment analysis...")
  
  # Step 1: Create output directory for GO enrichment if it doesn't exist
  go_dir <- file.path(output_dir, "GO_Enrichment")
  dir.create(go_dir, showWarnings = FALSE)
  
  # Step 2: Read each CSV file in the DGE results folder
  dge_files <- list.files(dge_dir, pattern = ".csv", full.names = TRUE)
  
  # Process each DGE file
  for (dge_file in dge_files) {
    message(paste("Processing DGE file:", dge_file))
    
    # Step 3: Read the DGE results
    dge_data <- read.csv(dge_file, row.names = 1)
    
    # Define upregulated and downregulated genes based on logFC and p-value thresholds
    dge_data$log10_pval <- -log10(dge_data$p_val)
    dge_data$expression <- "not significant"
    dge_data$expression[dge_data$log10_pval > -log10(pval_threshold) & dge_data$avg_log2FC > logfc_threshold] <- "upregulated"
    dge_data$expression[dge_data$log10_pval > -log10(pval_threshold) & dge_data$avg_log2FC < -logfc_threshold] <- "downregulated"
    
    # Subset upregulated and downregulated genes
    up_genes <- subset(dge_data, expression == "upregulated")
    up_genes$gene <- rownames(up_genes)
    down_genes <- subset(dge_data, expression == "downregulated")
    down_genes$gene <- rownames(down_genes)
    
    # Debugging: Print gene counts
    message(paste("Number of upregulated genes:", length(up_genes)))
    message(paste("Number of downregulated genes:", length(down_genes)))
    
    # Step 4: Remove missing or invalid gene symbols (empty or NA)
    up_genes_valid <- up_genes$gene[!is.na(up_genes$gene) & up_genes$gene != ""]
    down_genes_valid <- down_genes$gene[!is.na(down_genes$gene) & down_genes$gene != ""]
    
    # Debugging: Ensure there are no missing genes
    if (length(up_genes_valid) == 0) {
      message("No valid upregulated genes found.")
      next
    }
    if (length(down_genes_valid) == 0) {
      message("No valid downregulated genes found.")
      next
    }
    
    # Step 5: Check for valid gene symbols in the org.Mm.eg.db database
    valid_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")
    missing_up_genes <- setdiff(up_genes_valid, valid_symbols)
    missing_down_genes <- setdiff(down_genes_valid, valid_symbols)
    
    if (length(missing_up_genes) > 0) {
      message("Missing upregulated genes: ", paste(missing_up_genes, collapse = ", "))
    }
    if (length(missing_down_genes) > 0) {
      message("Missing downregulated genes: ", paste(missing_down_genes, collapse = ", "))
    }
    
    # Step 6: Map valid gene symbols to ENTREZ IDs
    up_genes_mapped <- bitr(
      up_genes_valid,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    )
    down_genes_mapped <- bitr(
      down_genes_valid,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    )
    
    # Step 7: Ensure proper mapping (no missing values)
    if (nrow(up_genes_mapped) == 0 || any(is.na(up_genes_mapped$ENTREZID))) {
      message("No valid Entrez IDs for upregulated genes.")
      next
    }
    if (nrow(down_genes_mapped) == 0 || any(is.na(down_genes_mapped$ENTREZID))) {
      message("No valid Entrez IDs for downregulated genes.")
      next
    }
    
    # Step 8: Run GO enrichment for upregulated genes
    if (nrow(up_genes_mapped) > 0) {
      message(paste("Running GO enrichment for upregulated genes in", dge_file))
      up_enrichment <- enrichGO(
        gene = up_genes_mapped$ENTREZID,
        OrgDb = org.Mm.eg.db,
        ont = "ALL",
        keyType = "ENTREZID",
        pAdjustMethod = 'BH',
        pvalueCutoff = pval_threshold
      )
      
      # Save the results and plot if enrichment is successful
      if (!is.null(up_enrichment) && nrow(up_enrichment@result) > 0) {
        up_file <- file.path(go_dir, paste0(tools::file_path_sans_ext(basename(dge_file)), "_upregulated_GO.csv"))
        write.csv(up_enrichment@result, up_file, row.names = FALSE)
        
        # Save the GO enrichment plot
        up_plot_file <- file.path(go_dir, paste0(tools::file_path_sans_ext(basename(dge_file)), "_upregulated_GO_plot.pdf"))
        ggsave(
          filename = up_plot_file,
          plot = enrichplot::dotplot(up_enrichment, split = "ONTOLOGY", showCategory = 6) +
            facet_grid(ONTOLOGY ~ ., scale = "free") +
            ggtitle("Upregulated GO Enrichment") +
            theme(text = element_text(size = 9),
                  axis.text.x = element_text(size = 9), 
                  axis.text.y = element_text(size = 9),
                  axis.title.y = element_text(size = 9),
                  axis.title.x = element_text(size = 9)),
          width = 10, height = 10
        )
      }
    }
    
    # Step 9: Run GO enrichment for downregulated genes
    if (nrow(down_genes_mapped) > 0) {
      message(paste("Running GO enrichment for downregulated genes in", dge_file))
      down_enrichment <- enrichGO(
        gene = down_genes_mapped$ENTREZID,
        OrgDb = org.Mm.eg.db,
        ont = "ALL",
        keyType = "ENTREZID",
        pAdjustMethod = 'BH',
        pvalueCutoff = pval_threshold
      )
      
      # Save the results and plot if enrichment is successful
      if (!is.null(down_enrichment) && nrow(down_enrichment@result) > 0) {
        down_file <- file.path(go_dir, paste0(tools::file_path_sans_ext(basename(dge_file)), "_downregulated_GO.csv"))
        write.csv(down_enrichment@result, down_file, row.names = FALSE)
        
        # Save the GO enrichment plot
        down_plot_file <- file.path(go_dir, paste0(tools::file_path_sans_ext(basename(dge_file)), "_downregulated_GO_plot.pdf"))
        ggsave(
          filename = down_plot_file,
          plot = enrichplot::dotplot(down_enrichment, split = "ONTOLOGY", showCategory = 6) +
            facet_grid(ONTOLOGY ~ ., scale = "free") +
            ggtitle("Downregulated GO Enrichment") +
            theme(text = element_text(size = 9),
                  axis.text.x = element_text(size = 9), 
                  axis.text.y = element_text(size = 9),
                  axis.title.y = element_text(size = 9),
                  axis.title.x = element_text(size = 9)),
          width = 10, height = 10
        )
      }
    }
  }
  message("GO enrichment results saved.")
}



# Call the GO enrichment function
perform_go_enrichment(dge_dir = "/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Inter_Timepoints/BT2/Day21_vs_Day14/DEGs", output_dir = "/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Inter_Timepoints/BT2/Day21_vs_Day14")

