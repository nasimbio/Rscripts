
library(clusterProfiler)
library(org.Mm.eg.db)        # Replace with org.Hs.eg.db for human
library(dplyr)
library(readr)
library(ggplot2)
library(forcats)

#  Step 1: Load the DE result file 
df <- read_csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Inter_Timepoints/C/Day21_vs_Day14/DEGs/combined.csv")  # Replace with your file

#  Step 2: Map SYMBOL to ENTREZ using bitr 
gene_map <- bitr(df$gene,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

df <- df %>%
  inner_join(gene_map, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(avg_log2FC), !is.na(ENTREZID))

#  Step 3: Run GSEA for Each Cell Type 
cell_types <- unique(df$Annotation)
gsea_results_list <- list()
combined_gsea_all <- list()

for (ct in cell_types) {
  message("Running GSEA for: ", ct)
  
  df_ct <- df %>% filter(Annotation == ct)
  
  gene_ranks <- df_ct$avg_log2FC
  names(gene_ranks) <- df_ct$ENTREZID
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  
  gsea_result <- tryCatch({
    gseGO(geneList = gene_ranks,
          OrgDb = org.Mm.eg.db,
          keyType = "ENTREZID",
          ont = "ALL",
          pvalueCutoff = 0.05,
          verbose = FALSE,
          nPermSimple = 10000,
          eps = 0)
  }, error = function(e) {
    warning(paste("GSEA failed for", ct, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result)) > 0) {
    gsea_df <- as.data.frame(gsea_result)
    gsea_df$cell_type <- ct
    
    # Save individual CSV
    write.csv(gsea_df, paste0("/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Inter_Timepoints/C/Day21_vs_Day14/GSEA/GSEA_", ct, "_C_Day21_vs_Day14.csv"), row.names = FALSE)
    
    gsea_results_list[[ct]] <- gsea_result
    combined_gsea_all[[ct]] <- gsea_df
  } else {
    message("⚠️ No enriched terms found for ", ct)
  }
}

#  Step 4: Combine and Save All Results 
combined_gsea_df <- bind_rows(combined_gsea_all)%>%
  filter(!is.na(NES), !is.na(p.adjust))
write_csv(combined_gsea_df, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Inter_Timepoints/C/Day21_vs_Day14/GSEA/Combined_GSEA_C_Day21_vs_Day14.csv")

# Step 5: Select Top 5 GO Terms per Cell Type 
top_pathways <- combined_gsea_df %>%
  group_by(cell_type) %>%
  slice_min(p.adjust, n = 10) %>%
  ungroup()

#  Step 6: Plot with Color = p.adjust 
plot <- ggplot(top_pathways,
               aes(x = fct_reorder(Description, NES),
                   y = NES,
                   fill = p.adjust)) +
  geom_col(show.legend = TRUE, position = position_dodge(width = 0.8)) +
  #scale_fill_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  coord_flip() +
  facet_wrap(~cell_type, scales = "free_y") +
  labs(title = "Top Enriched GO Terms per Cell Type (GSEA)",
       x = "GO Term",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")+
  theme(strip.text = element_text(face = "bold"))+
  scale_fill_gradient(low = "red", high = "blue")

#  Step 7: Save the Plot 
ggsave("/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Inter_Timepoints/C/Day21_vs_Day14/GSEA/Top_GSEA_Terms_C_Day21_vs_Day14.pdf", plot, width = 23, height = 10)






