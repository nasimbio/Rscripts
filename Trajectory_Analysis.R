title: " Trajectory Analysis on T Immune Cells"
author: "Nasim Rahmatpour"
date: "2025-09-18"
output: html_document
---
  
#libraries
set.seed(1234)
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)  

### Task 1: Seurat → Monocle3 cell_data_set Conversion and Trajectory Learning

#give simple name to seurat
seurat_obj <- `3-Data_Annotation_filtered_data`  

#convert seuart to cell data set object
cds <- as.cell_data_set(seurat_obj)
cds

#get the cell metadata, get the gene metadata
colData(cds)
fData(cds)
rownames(fData(cds))[1:20]
fData(cds)$gene_short_name <- rownames(fData(cds))

#get the count
counts(cds)

#using the cluster info of the seurat
#assign the partitions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition

#assign annotation info
list_cluster <- seurat_obj@active.ident
cds@clusters$UMAP$clusters <- list_cluster

#assign umap coordinates
pdf( "/home/rstudio/projects/before_trajectory.pdf", width=15, height =10)
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "New_Annotation",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan', 'white')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names
dev.off()

#learn trajcetory
cds <- learn_graph(cds, use_partition = FALSE)
pdf( "/home/rstudio/projects/after_trajectory.pdf", width=15, height =10)
plot_cells(cds,
           color_cells_by = 'New_Annotation',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)
dev.off()

#order the cells in pseudotime
pdf( "/home/rstudio/projects/after_trajectory_ordered.pdf", width=15, height =10)
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == "C1_Sox5_Stem-like"]))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = TRUE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)
dev.off()

#cells orderd by monocle 3 peseudotime
pdf( "/home/rstudio/projects/pseudotime_boxplot.pdf", width=10, height =10)
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(New_Annotation, monocle3_pseudotime, median), fill = New_Annotation)) +
  geom_boxplot() +
  ylab("Cell types")
dev.off()

#find the genes that their expression is important in trajectory analysis
deg_cells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_cells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

FeaturePlot(seurat_obj, features = c('Tcea1', 'Atp6v1h', 'Rb1cc1'))


# visualizing pseudotime in seurat
pdf( "/home/rstudio/projects/pseudotime_seurat.pdf", width=12, height =10)
seurat_obj$pseudotime <- pseudotime(cds)
Idents(seurat_obj) <- seurat_obj$New_Annotation
FeaturePlot(seurat_obj, features = "pseudotime", label = T)
dev.off()

#save the seurat object that has trajectory info in the metadata and use it later for boxplot
saveRDS(seurat_obj, "/home/rstudio/projects/seurat_obj_trajectory.rds")



### Task 2: Pseudotime Analysis and Visualization
#box plot
pdf( "/home/rstudio/projects/trajectory_boxplot.pdf", width=15, height =10)
df <- seurat_obj_trajectory@meta.data %>%
  transmute(
    pseudotime = pseudotime,
    celltype   = New_Annotation,
    Timepoint  = factor(Timepoint, levels = c("Day14","Day21")),
    Treatment  = factor(Treatment, levels = c("C","B","B+T1","B+T2"))
  ) %>%
  filter(!is.na(pseudotime), !is.na(celltype), !is.na(Timepoint), !is.na(Treatment))
cell_order <- c(
  "C1_Sox5_Stem-like","Proliferating_Mki67","C2_Id3_Sell_Naive",
  "Th17","CD4_Th2","T-regs_Foxp3","CD8_T-cells","C11_NK"
)
df <- df %>% mutate(celltype = factor(celltype, levels = cell_order))
ggplot(df, aes(x = Treatment, y = pseudotime, fill = Treatment)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, color = "grey25") +
  facet_grid(rows = vars(Timepoint), cols = vars(celltype), scales = "free_y")+
  
  #facet_wrap(~ celltype, ncol = 3, scales = "free_y") +   # one panel per cell type
  scale_fill_manual(values = pal, drop = FALSE) +
  labs(x = "Treatment", y = "Monocle3 pseudotime")+
  #title = "Pseudotime across Treatments, faceted by Cell Type") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
dev.off()