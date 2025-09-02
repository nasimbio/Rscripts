

#get both positive and negative log2FC
library(dplyr)
DGE_Day14_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Intra_Timepoints/Day14/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

top30_genes_day14_BvsC <- DGE_Day14_BvsC %>%
  group_by(Annotation) %>%
  slice_max(avg_log2FC, n = 15) %>%  # ✅ Top 15 Upregulated (Highest log2FC)
  bind_rows(DGE_Day14_BvsC %>%
              group_by(Annotation) %>%
              slice_min(avg_log2FC, n = 15))  # ✅ Top 15 Downregulated (Lowest log2FC)


DGE_Day14_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Intra_Timepoints/Day14/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

top30_genes_day14_BT1vsC <- DGE_Day14_BT1vsC %>%
  group_by(Annotation) %>%
  slice_max(avg_log2FC, n = 15) %>%  # ✅ Top 15 Upregulated (Highest log2FC)
  bind_rows(DGE_Day14_BT1vsC %>%
              group_by(Annotation) %>%
              slice_min(avg_log2FC, n = 15))  # ✅ Top 15 Downregulated (Lowest log2FC)

DGE_Day14_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Intra_Timepoints/Day14/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

top30_genes_day14_BT2vsC <- DGE_Day14_BT2vsC %>%
  group_by(Annotation) %>%
  slice_max(avg_log2FC, n = 15) %>%  # ✅ Top 15 Upregulated (Highest log2FC)
  bind_rows(DGE_Day14_BT2vsC %>%
              group_by(Annotation) %>%
              slice_min(avg_log2FC, n = 15))  # ✅ Top 15 Downregulated (Lowest log2FC)

#top 30 combined
top30_combined <- do.call(rbind, list(top30_genes_day14_BvsC, top30_genes_day14_BT1vsC, top30_genes_day14_BT2vsC))
colnames(top30_combined)[colnames(top30_combined) == "Treatment"] <- "Comparison"

top30_combined_filter <- top30_combined[!duplicated(top30_combined$gene), ]%>% pull(gene)



library(dplyr)
DGE_Day21_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Intra_Timepoints/Day21/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

top30_genes_day21_BvsC <- DGE_Day21_BvsC %>%
  group_by(Annotation) %>%
  slice_max(avg_log2FC, n = 15) %>%  # ✅ Top 15 Upregulated (Highest log2FC)
  bind_rows(DGE_Day21_BvsC %>%
              group_by(Annotation) %>%
              slice_min(avg_log2FC, n = 15))  # ✅ Top 15 Downregulated (Lowest log2FC)

# 
# 
DGE_Day21_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Intra_Timepoints/Day21/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

top30_genes_day21_BT1vsC <- DGE_Day21_BT1vsC %>%
  group_by(Annotation) %>%
  slice_max(avg_log2FC, n = 15) %>%  # ✅ Top 15 Upregulated (Highest log2FC)
  bind_rows(DGE_Day21_BT1vsC %>%
              group_by(Annotation) %>%
              slice_min(avg_log2FC, n = 15))  # ✅ Top 15 Downregulated (Lowest log2FC)


DGE_Day21_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/DGE_Intra_Timepoints/Day21/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

top30_genes_day21_BT2vsC <- DGE_Day21_BT2vsC %>%
  group_by(Annotation) %>%
  slice_max(avg_log2FC, n = 15) %>%  # ✅ Top 15 Upregulated (Highest log2FC)
  bind_rows(DGE_Day21_BT2vsC %>%
              group_by(Annotation) %>%
              slice_min(avg_log2FC, n = 15))  # ✅ Top 15 Downregulated (Lowest log2FC)




#top 30 combine

top30_combined <- do.call(rbind, list(top30_genes_day21_BvsC, top30_genes_day21_BT1vsC, top30_genes_day21_BT2vsC))
colnames(top30_combined)[colnames(top30_combined) == "Treatment"] <- "Comparison"



#Pheatmap
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/no_TCR_subclustering/Pheatmap_Day21_avg_log2FC.pdf'), width = 15, height = 20)

library(pheatmap)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Remove duplicates
df_unique <- top30_combined[!duplicated(top30_combined$gene), ]
df_selected <- df_unique[, c("gene", "Comparison", "avg_log2FC", "Annotation")]

# **Step 4: Reshape Data (Pivot Table to Wide Format)**
expr_matrix <- df_selected %>%
  pivot_wider(names_from = Comparison, values_from = avg_log2FC, values_fill = list(avg_log2FC = 0)) %>%
  column_to_rownames("gene")

# **Step 4.1: Keep Only Treatment Columns (Remove Annotation & Extra Columns)**
expr_matrix <- expr_matrix[, c("B_vs_C", "BT1_vs_C", "BT2_vs_C")]  # ✅ Keep only relevant treatments

# **Step 4.2: Convert Data to a Numeric Matrix**
#expr_matrix <- as.matrix(expr_matrix)  # ✅ Ensure numeric matrix
#expr_matrix <- apply(expr_matrix, 2, as.numeric)  # ✅ Force numeric conversion

# Step 5: Ensure Genes Are Ordered by Annotation Before Clustering
df_selected <- df_selected %>% arrange(Annotation)
expr_matrix <- expr_matrix[match(df_selected$gene, rownames(expr_matrix)), ]

# Step 6: Define Annotations
row_annotation <- df_selected %>%
  distinct(gene, Annotation) %>%
  column_to_rownames("gene")  # Annotation info for genes

# **Step 6.1: Define Column Annotations with Only Selected Treatments**
col_annotation <- data.frame(Comparison = colnames(expr_matrix))  # ✅ Use only B_vs_C, BT1_vs_C, BT2_vs_C
rownames(col_annotation) <- colnames(expr_matrix)  # ✅ Ensure correct rownames

# Step 7: Define Colors
annotation_colors <- list(
  Annotation = setNames(colorRampPalette(brewer.pal(8, "Set3"))(length(unique(df_selected$Annotation))),
                        unique(df_selected$Annotation)),
  Comparison = c("B_vs_C" = "#1f78b4", "BT1_vs_C" = "#33a02c", "BT2_vs_C" = "#e31a1c")
)


log2FC_min <- -10  # Set minimum log2FC (adjust based on your data)
log2FC_max <- 10   # Set maximum log2FC

#Define color breaks (ensuring 0 is always white)
breaks_list <- seq(log2FC_min, log2FC_max, length.out = 100)

pheatmap(
  expr_matrix,
  name = "Log2 Fold Change",
  cluster_rows = TRUE,  
  cluster_cols = TRUE,  
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  scale = "none",
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 3,
  fontsize_col = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),  # Fixed color palette
  breaks = breaks_list,  # ✅ Ensures the same color mapping for Day 14 & Day 21
  main = "Heatmap for Day21 (avg_log2FC of Top30 Genes from DGE in Each Cell Type)",
  split = row_annotation$Annotation
)

# Step 8: Generate Heatmap
# pheatmap(
#   expr_matrix,
#   name = "Log2 Fold Change",
#   cluster_rows = TRUE,  # Cluster genes inside annotation groups
#   cluster_cols = TRUE,  # Cluster treatments
#   clustering_distance_cols = "correlation",  # ✅ Change to correlation-based distance
#   clustering_method = "average",
#   scale = "none",  # Use actual avg_log2FC values
#   annotation_row = row_annotation,  # Show annotation on rows
#   annotation_col = col_annotation,  # Show treatment annotation
#   annotation_colors = annotation_colors,
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   fontsize_row = 3,  # Adjust gene font size
#   fontsize_col = 10,  # Adjust treatment font size
#   color = colorRampPalette(c("blue", "white", "red"))(100),
#   main = "Heat map for D14 (avg_log2FC of Top30 Genes from DGE in Each Cell Type)",
#   split = row_annotation$Annotation
# )

dev.off()

pheatmap(
  expr_matrix,
  name = "Log2 Fold Change",
  cluster_rows = FALSE,  
  cluster_cols = TRUE,  
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  scale = "none",
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 3,
  fontsize_col = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  gaps_row = which(diff(as.numeric(factor(row_annotation$Annotation))) != 0),  # ✅ Add visible gaps for Annotation groups
  main = "Heat map for D14 (avg_log2FC of Top30 Genes from DGE in Each Cell Type)",
  split = row_annotation$Annotation
)

#heatmap for comparison across selected DGE results for TRM-like, vs C 

library(dplyr)
DGE_Day14_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

combined <- do.call(rbind, list(DGE_Day14_BvsC, DGE_Day14_BT1vsC, DGE_Day14_BT2vsC, DGE_Day21_BvsC, DGE_Day21_BT1vsC, DGE_Day21_BT2vsC))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"

# combined_test <- bind_rows(DGE_Day14_BvsC, DGE_Day14_BT1vsC, DGE_Day14_BT2vsC, DGE_Day21_BvsC, DGE_Day21_BT1vsC, DGE_Day21_BT2vsC)

combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Lag3", "Tigit", "Havcr2", "Tox", "Cxcr6", "Entpd1", "Cxcl13", "Prf1", "Gzma", "Gzmb", "Gzmh", "Gmmxk", "Nkg7", "Itgae", "Itga1", "Id2", "Id3")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("B_vs_C" = "#33a02c", "BT1_vs_C" = "#6a3d9a", "BT2_vs_C" = "#ff7f00")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_TRM_like_vs_C.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of Terminally exhausted-like, cytotoxic-like, TRM-like TEX Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_TRM_like_vs_C.csv" )


#heatmap for comparison across selected DGE results for TRM-like, vs B 
library(dplyr)
DGE_Day14_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day14_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_B/DEGs/combined.csv", row.names = NULL)


combined <- do.call(rbind, list(DGE_Day14_BT1vsB, DGE_Day14_BT2vsB, DGE_Day21_BT1vsB, DGE_Day21_BT2vsB))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Lag3", "Tigit", "Havcr2", "Tox", "Cxcr6", "Entpd1", "Cxcl13", "Prf1", "Gzma", "Gzmb", "Gzmh", "Gmmxk", "Nkg7", "Itgae", "Itga1", "Id2", "Id3")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("BT1_vs_B" = "#33a02c", "BT2_vs_B" = "#6a3d9a")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_TRM_like_vs_B.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of Terminally exhausted-like, cytotoxic-like, TRM-like TEX Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_TRM_like_vs_B.csv" )


#heatmap for comparison across selected DGE results for Effector-like, vs C 

library(dplyr)
DGE_Day14_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

combined <- do.call(rbind, list(DGE_Day14_BvsC, DGE_Day14_BT1vsC, DGE_Day14_BT2vsC, DGE_Day21_BvsC, DGE_Day21_BT1vsC, DGE_Day21_BT2vsC))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"

# combined_test <- bind_rows(DGE_Day14_BvsC, DGE_Day14_BT1vsC, DGE_Day14_BT2vsC, DGE_Day21_BvsC, DGE_Day21_BT1vsC, DGE_Day21_BT2vsC)

combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Ifng", "Cd69", "Ccl3", "Ccl4", "Tnf", "Nr4a2", "Nr4a1", "Tnfsf9", "Icos")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("B_vs_C" = "#33a02c", "BT1_vs_C" = "#6a3d9a", "BT2_vs_C" = "#ff7f00")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Effector_like_vs_C.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of Effector-like Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Effector_like_vs_C.csv" )


#heatmap for comparison across selected DGE results for Effector-like, vs B

library(dplyr)
DGE_Day14_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day14_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_B/DEGs/combined.csv", row.names = NULL)


combined <- do.call(rbind, list(DGE_Day14_BT1vsB, DGE_Day14_BT2vsB, DGE_Day21_BT1vsB, DGE_Day21_BT2vsB))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Ifng", "Cd69", "Ccl3", "Ccl4", "Tnf", "Nr4a2", "Nr4a1", "Tnfsf9", "Icos")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("BT1_vs_B" = "#33a02c", "BT2_vs_B" = "#6a3d9a")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Effector_like_vs_B.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of Effector-like Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Effector_like_vs_B.csv" )



#heatmap for comparison across selected DGE results for GZMK+ effector memory-like, vs C 

library(dplyr)
DGE_Day14_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

combined <- do.call(rbind, list(DGE_Day14_BvsC, DGE_Day14_BT1vsC, DGE_Day14_BT2vsC, DGE_Day21_BvsC, DGE_Day21_BT1vsC, DGE_Day21_BT2vsC))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Gzmk", "Gzmb", "Il7r", "Cxcr4", "Cxcr3", "Klrg1")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("B_vs_C" = "#33a02c", "BT1_vs_C" = "#6a3d9a", "BT2_vs_C" = "#ff7f00")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_GZMK_effector_memory_like_vs_C.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of GZMK+ effector memory-like Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_GZMK_effector_memory_lik_vs_C.csv" )


#heatmap for comparison across selected DGE results for GZMK+ effector memory-like, vs B

library(dplyr)
DGE_Day14_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day14_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_B/DEGs/combined.csv", row.names = NULL)


combined <- do.call(rbind, list(DGE_Day14_BT1vsB, DGE_Day14_BT2vsB, DGE_Day21_BT1vsB, DGE_Day21_BT2vsB))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Gzmk", "Gzmb", "Il7r", "Cxcr4", "Cxcr3", "Klrg1")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("BT1_vs_B" = "#33a02c", "BT2_vs_B" = "#6a3d9a")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_GZMK_effector_memory_like_vs_B.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of GZMK+ effector memory-like Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_GZMK_effector_memory_like_vs_B.csv" )


#heatmap for comparison across selected DGE results for TPEX-like TEX, vs C 

library(dplyr)
DGE_Day14_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

combined <- do.call(rbind, list(DGE_Day14_BvsC, DGE_Day14_BT1vsC, DGE_Day14_BT2vsC, DGE_Day21_BvsC, DGE_Day21_BT1vsC, DGE_Day21_BT2vsC))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Tcf7", "Ccr7", "Il7r")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("B_vs_C" = "#33a02c", "BT1_vs_C" = "#6a3d9a", "BT2_vs_C" = "#ff7f00")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_TPEX_like_TEX_vs_C.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of TPEX-like TEX Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_TPEX_like_TEX_vs_C.csv" )



#heatmap for comparison across selected DGE results for TPEX-like TEX, vs B 
library(dplyr)
DGE_Day14_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day14_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_B/DEGs/combined.csv", row.names = NULL)


combined <- do.call(rbind, list(DGE_Day14_BT1vsB, DGE_Day14_BT2vsB, DGE_Day21_BT1vsB, DGE_Day21_BT2vsB))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Tcf7", "Ccr7", "Il7r")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("BT1_vs_B" = "#33a02c", "BT2_vs_B" = "#6a3d9a")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_TPEX_like_TEX_vs_B.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of TPEX-like TEX Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_TPEX_like_TEX_vs_B.csv")


#heatmap for comparison across selected DGE results for cycling TEX, vs C 
library(dplyr)
DGE_Day14_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

combined <- do.call(rbind, list(DGE_Day14_BvsC, DGE_Day14_BT1vsC, DGE_Day14_BT2vsC, DGE_Day21_BvsC, DGE_Day21_BT1vsC, DGE_Day21_BT2vsC))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Mki67", "Ccnb2", "Ccna2", "Cdca5", "Top2a", "Cdk1")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("B_vs_C" = "#33a02c", "BT1_vs_C" = "#6a3d9a", "BT2_vs_C" = "#ff7f00")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Cycling_TEX_vs_C.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of Cycling TEX Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Cycling_TEX_vs_C.csv" )



#heatmap for comparison across selected DGE results for cycling TEX, vs B
library(dplyr)
DGE_Day14_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day14_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_B/DEGs/combined.csv", row.names = NULL)


combined <- do.call(rbind, list(DGE_Day14_BT1vsB, DGE_Day14_BT2vsB, DGE_Day21_BT1vsB, DGE_Day21_BT2vsB))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Mki67", "Ccnb2", "Ccna2", "Cdca5", "Top2a", "Cdk1")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("BT1_vs_B" = "#33a02c", "BT2_vs_B" = "#6a3d9a")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Cycling_TEX_vs_B.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of Cycling TEX Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Cycling_TEX_vs_B.csv")


#heatmap for comparison across selected DGE results for IFN-stimulated TEX, vs C
library(dplyr)
DGE_Day14_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

combined <- do.call(rbind, list(DGE_Day14_BvsC, DGE_Day14_BT1vsC, DGE_Day14_BT2vsC, DGE_Day21_BvsC, DGE_Day21_BT1vsC, DGE_Day21_BT2vsC))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Mx1", "Ifi6", "Oas1", "Isg15", "Stat1", "Mx2", "Irf7", "Ifitm1")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("B_vs_C" = "#33a02c", "BT1_vs_C" = "#6a3d9a", "BT2_vs_C" = "#ff7f00")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_IFN_stimulated_TEX_vs_C.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of IFN-stimulated TEX Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_IFN_stimulated_TEX_vs_C.csv" )


#heatmap for comparison across selected DGE results for IFN-stimulated TEX, vs B

library(dplyr)
DGE_Day14_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day14_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_B/DEGs/combined.csv", row.names = NULL)


combined <- do.call(rbind, list(DGE_Day14_BT1vsB, DGE_Day14_BT2vsB, DGE_Day21_BT1vsB, DGE_Day21_BT2vsB))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Mx1", "Ifi6", "Oas1", "Isg15", "Stat1", "Mx2", "Irf7", "Ifitm1")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("BT1_vs_B" = "#33a02c", "BT2_vs_B" = "#6a3d9a")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_IFN_stimulated_TEX_vs_B.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of IFN-stimulated TEX Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_IFN_stimulated_TEX_vs_B.csv")


#heatmap for comparison across selected DGE results for Naive and or TCM-like, vs C
library(dplyr)
DGE_Day14_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

combined <- do.call(rbind, list(DGE_Day14_BvsC, DGE_Day14_BT1vsC, DGE_Day14_BT2vsC, DGE_Day21_BvsC, DGE_Day21_BT1vsC, DGE_Day21_BT2vsC))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Ilr7", "Lef1", "Tcf7", "S1pr1", "Sell")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("B_vs_C" = "#33a02c", "BT1_vs_C" = "#6a3d9a", "BT2_vs_C" = "#ff7f00")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Naive_or_TCM_like_vs_C.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of Naive or TCM-like Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Naive_or_TCM_like_vs_C.csv" )


#heatmap for comparison across selected DGE results for Naive and or TCM-like, vs B

library(dplyr)
DGE_Day14_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day14_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_B/DEGs/combined.csv", row.names = NULL)


combined <- do.call(rbind, list(DGE_Day14_BT1vsB, DGE_Day14_BT2vsB, DGE_Day21_BT1vsB, DGE_Day21_BT2vsB))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Ilr7", "Lef1", "Tcf7", "S1pr1", "Sell")


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("BT1_vs_B" = "#33a02c", "BT2_vs_B" = "#6a3d9a")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Naive_or_TCM_like_vs_B.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC of Naive or TCM-like Genes (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_Naive_or_TCM_like_vs_B.csv")


#heatmap for comparison across selected DGE results for gene list, vs C
library(dplyr)
DGE_Day14_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day14_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BvsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/B_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT1vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_C/DEGs/significant_combined.csv", row.names = NULL)

DGE_Day21_BT2vsC <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_C/DEGs/significant_combined.csv", row.names = NULL)

combined <- do.call(rbind, list(DGE_Day14_BvsC, DGE_Day14_BT1vsC, DGE_Day14_BT2vsC, DGE_Day21_BvsC, DGE_Day21_BT1vsC, DGE_Day21_BT2vsC))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Birc3", "Cflar", "Batf3", "Junb", "Gos2", "Ell2", "Mir155hg", "Tnfrsf8", "Prkcdbp", "Hmsd", "Sh3bp5", "Tnip2", "Traf1", "Ccr7", "Il4l1", "Ddit4", "Arid5a", "Hla-dpa1", "Hla-dpb1", "Hla-dqa1", "Hla-dra", "Hla-drb1", "Cd74", "Enpp2", "Adcy1", "Cxcl10", "Ank3") 


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("B_vs_C" = "#33a02c", "BT1_vs_C" = "#6a3d9a", "BT2_vs_C" = "#ff7f00")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_gene_list_vs_C.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_gene_list_vs_C.csv" )


#heatmap for comparison across selected DGE results for gene list, vs B
library(dplyr)
DGE_Day14_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day14_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day14/BT2_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT1vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT1_vs_B/DEGs/combined.csv", row.names = NULL)

DGE_Day21_BT2vsB <- read.csv("/home/nasim/amgen.abio/amgen_outputs_with_negatives/CustomizedAnalysis_Aggregated_Clusters/DGE_Intra_Timepoints/Day21/BT2_vs_B/DEGs/combined.csv", row.names = NULL)


combined <- do.call(rbind, list(DGE_Day14_BT1vsB, DGE_Day14_BT2vsB, DGE_Day21_BT1vsB, DGE_Day21_BT2vsB))
colnames(combined)[colnames(combined) == "Treatment"] <- "Comparison"


combined <- combined %>%
  mutate(Comparison_ID = paste(Timepoint, Comparison, sep = "_"))



genes_of_interest <- c("Birc3", "Cflar", "Batf3", "Junb", "Gos2", "Ell2", "Mir155hg", "Tnfrsf8", "Prkcdbp", "Hmsd", "Sh3bp5", "Tnip2", "Traf1", "Ccr7", "Il4l1", "Ddit4", "Arid5a", "Hla-dpa1", "Hla-dpb1", "Hla-dqa1", "Hla-dra", "Hla-drb1", "Cd74", "Enpp2", "Adcy1", "Cxcl10", "Ank3") 


filtered_df <- combined %>%
  filter(gene %in% genes_of_interest)


summarized <- filtered_df %>%
  group_by(gene, Comparison_ID, Timepoint, Comparison) %>%  # keep Timepoint info
  slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()


expr_matrix <- summarized %>%
  select(gene, Comparison_ID, avg_log2FC) %>%
  pivot_wider(names_from = Comparison_ID, values_from = avg_log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

col_anno <- summarized %>%
  distinct(Comparison_ID, Timepoint, Comparison) %>%
  column_to_rownames("Comparison_ID")  # rows = columns in matri

anno_colors <- list(
  Timepoint = c("Day14" = "#1f78b4", "Day21" = "#e31a1c"),
  Comparison = c("BT1_vs_B" = "#33a02c", "BT2_vs_B" = "#6a3d9a")
)
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_gene_list_vs_B.pdf'), width = 12, height = 10)
pheatmap(
  as.matrix(expr_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = col_anno,  # ✅ shows Timepoint & Comparison
  annotation_colors = anno_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-5, 5, length.out = 100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "avg_log2FC (Top signal per Comparison)"
)
dev.off()
write.csv(expr_matrix, "/home/nasim/amgen.abio/amgen_outputs_with_negatives/Pheatmap_gene_list_vs_B.csv")


