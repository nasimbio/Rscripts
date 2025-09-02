# ------------------------------------------------------------------------------
#title: "10x Visium HD Spatial Transcriptomics — Human Lung"
#Author: Nasim Rahmatpour
#Date: "2024-09-04"
# ------------------------------------------------------------------------------


# Project Overview


#Spatial transcriptomic profiling of human lung tissue was performed using the **Visium HD** platform, which offers subcellular resolution. Standard preprocessing, dimensionality reduction, and clustering were conducted to capture spatial gene expression patterns and identify distinct tissue regions.  
#To resolve the cellular composition of each spatial location, we employed the **RCTD** algorithm, leveraging an external single-cell reference atlas for robust deconvolution of Visium HD spots.


---

###Task 1: Preprocessing, Clustering, and Regional Annotation 
  
# packages required for Visium HD
install.packages("hdf5r")
install.packages("arrow")

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(arrow)

localdir <- "/home/nasim/4dmt.abio/NGS087_3/outs/binned_outputs/square_008um/"
NGS087_3 <- Load10X_Spatial(data.dir = localdir, file= "filtered_feature_bc_matrix.h5", assay="Spatial.008um")

dir.create("/home/nasim/4dmt.abio/NGS087_3_outputs")
pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_violin_count.pdf")
vln.plot <- VlnPlot(NGS087_3 , features = "nCount_Spatial.008um", pt.size= 0) + theme(axis.text = element_text(size = 10)) + NoLegend()
count.plot <- SpatialFeaturePlot(NGS087_3 , features = "nCount_Spatial.008um", pt.size= 12) + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot|count.plot
dev.off()

DefaultAssay(NGS087_3) <- "Spatial.008um"
NGS087_3  <- NormalizeData(NGS087_3)

#sketch data for 5000
DefaultAssay(NGS087_3) <- "Spatial.008um"
NGS087_3 <- FindVariableFeatures(NGS087_3)
NGS087_3 <- ScaleData(NGS087_3)
#default number of cells are 5000
# change to 50,0000 cells for test, create a new 'sketch' assay
NGS087_3 <- SketchData(
  object = NGS087_3,
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch analysis to sketched cells
DefaultAssay(NGS087_3) <- "sketch"

# perform clustering workflow
NGS087_3 <- FindVariableFeatures(NGS087_3)
NGS087_3 <- ScaleData(NGS087_3)
NGS087_3 <- RunPCA(NGS087_3 , assay = "sketch", reduction.name = "pca.sketch")
pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Elbow_plot.pdf")
P1 <- ElbowPlot(NGS087_3, reduction = "pca.sketch")
print(P1)
dev.off()
NGS087_3 <- FindNeighbors(NGS087_3 , assay = "sketch", reduction = "pca.sketch", dims = 1:20)
NGS087_3 <- FindClusters(NGS087_3 , cluster.name = "seurat_cluster.sketched", resolution = 1)
NGS087_3 <- RunUMAP(NGS087_3, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:20)

options(future.globals.maxSize = 5 * 1024^3)  # 5 GB



#Project the result of skteching to the rest of the data
NGS087_3 <- ProjectData(
  object = NGS087_3,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:20,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_umap.pdf")
DefaultAssay(NGS087_3) <- "sketch"
Idents(NGS087_3) <- "seurat_cluster.sketched"
p1 <- DimPlot(NGS087_3, reduction = "umap.sketch", label=T,repel=T) + ggtitle("Sketched clustering (5000 cells)") #+ theme(legend.position = "none")

# switch to full dataset
DefaultAssay(NGS087_3) <- "Spatial.008um"
Idents(NGS087_3) <- "seurat_cluster.projected"
p2 <- DimPlot(NGS087_3, reduction = "full.umap.sketch", label=T,repel=T) + ggtitle("Projected clustering (full dataset)") #+ theme(legend.position = "none")

p1 
p2
dev.off()

#make SpatialDim plot of all clusters
pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_SpatialDimPlot_all_clusters.pdf", width = 20,height = 15) 
SpatialDimPlot(NGS087_3, label = T, repel = T, label.size = 4, pt.size.factor = 10) +
  theme(
    legend.text = element_text(size = 15), # Increase legend item text size
    legend.title = element_text(size = 15) # Increase legend title size
  )

#make SpatialDim plot of some clusters
pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_SpatialDimPlot_individual_clusters.pdf") 
Idents(NGS087_3) <- "seurat_cluster.projected"
cells <- CellsByIdentities(NGS087_3, idents = c(0, 1, 2, 3, 4, 5))
p1 <- SpatialDimPlot(NGS087_3,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = T, pt.size.factor = 10 
) + NoLegend()

cells <- CellsByIdentities(NGS087_3, idents = c(6, 7, 8, 9, 10, 11))
p2 <- SpatialDimPlot(NGS087_3,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = T, pt.size.factor = 10 
) + NoLegend()

cells <- CellsByIdentities(NGS087_3, idents = c(12, 13, 14, 15, 16))
p3 <- SpatialDimPlot(NGS087_3,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = T, pt.size.factor = 10
) + NoLegend()


p1
p2
p3

dev.off()

#heatmap to show the top genes of each cluster
# Crete downsampled object to make visualization better
DefaultAssay(NGS087_3) <- "Spatial.008um"
Idents(NGS087_3) <- "seurat_cluster.projected"
object_subset <- subset(NGS087_3, cells = Cells(NGS087_3[["Spatial.008um"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_top10_Heatmap.pdf",  width = 20,height = 15) 
object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top10$gene)

p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top10$gene, size = 10)+
  theme(axis.text = element_text(size = 7)) + 
  theme(axis.text.y = element_text(size = 12))+ # Make y-axis text larger
  NoLegend()
p
dev.off()

write.csv(markers,file = paste0("/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_RNA_markers.csv"),row.names = TRUE)

top10.markers_RNA <- markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top10 <- c()
for(c in unique(top10.markers_RNA$cluster)){
  tmp <- top10.markers_RNA[top10.markers_RNA$cluster==c,]$gene
  top10 <- cbind(top10,tmp)
}
colnames(top10) <- paste0('C_',unique(top10.markers_RNA$cluster))
write.csv(top10,"/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_top10.RNA_markers.csv")

saveRDS(NGS087_3, '/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008.rds')


#feature plot of top genes
pdf(paste0("/home/nasim/4dmt.abio/NGS087_3_outputs/","Spatial.008_top2_markers_featureplot.pdf"),width = 10,height = 10)
DefaultAssay(NGS087_3) <- "Spatial.008um"

P1=FeaturePlot(NGS087_3 , features = c("MT-ND4L", "TMSB4X", "IGKC", "IGHG3", "SCGB3A1", "MUC5AC", "SCGB1A1", "CFDP1", "KRT5"),  ncol = 3)

P2=FeaturePlot(NGS087_3 , features = c("MUC16", "BPIFB1", "PLEKHS1", "ADHFE1", "FGFR4", "CD24", "DYNC1LI2", "IGHA1", "PLXDC2"),  ncol = 3)

P3=FeaturePlot(NGS087_3 , features = c("JCHAIN", "IGKC", "SEC31A", "PTPN1", "IGHD", "TXNDC5", "TNFRSF21", "FAT2", "PTPRS"),  ncol = 3)

P4=FeaturePlot(NGS087_3 , features = c("FOSB", "IGHA1", "IGKC", "DPYSL2", "CLDN18" ),  ncol = 3)

print(P1)
print(P2)
print(P3)
print(P4)
dev.off()

#Spatialfeatureplots
pdf(paste0("/home/nasim/4dmt.abio/NGS087_3_outputs/","Spatial.008_top_markers_spatialfeature1.pdf"),width = 10,height = 10)
p1 <- SpatialFeaturePlot(NGS087_3, features =c("MT-ND4L", "TMSB4X", "IGKC", "IGHG3"),alpha = c(0.1, 1), ncol = 2, pt.size= 12) 
p1
dev.off()


pdf(paste0("/home/nasim/4dmt.abio/NGS087_3_outputs/","Spatial.008_top_markers_spatialfeature2.pdf"),width = 10,height = 10)
p2 <- SpatialFeaturePlot(NGS087_3, features =c("SCGB3A1", "MUC5AC", "SCGB1A1", "CFDP1"),alpha = c(0.1, 1), ncol = 2, pt.size= 12) 
p2
dev.off()

pdf(paste0("/home/nasim/4dmt.abio/NGS087_3_outputs/","Spatial.008_top_markers_spatialfeature3.pdf"),width = 10,height = 10)
p3 <- SpatialFeaturePlot(NGS087_3, features =c( "KRT5", "MUC16", "BPIFB1", "PLEKHS1"),alpha = c(0.1, 1), ncol = 2, pt.size= 12) 
p3
dev.off()

pdf(paste0("/home/nasim/4dmt.abio/NGS087_3_outputs/","Spatial.008_top_markers_spatialfeature4.pdf"),width = 10,height = 10)
p4 <- SpatialFeaturePlot(NGS087_3, features =c( "ADHFE1", "FGFR4", "CD24", "DYNC1LI2"),alpha = c(0.1, 1), ncol = 2, pt.size= 12) 
p4
dev.off()


pdf(paste0("/home/nasim/4dmt.abio/NGS087_3_outputs/","Spatial.008_top_markers_spatialfeature5.pdf"),width = 10,height = 10)
p5 <- SpatialFeaturePlot(NGS087_3, features =c("IGHA1", "PLXDC2","JCHAIN", "IGKC"),alpha = c(0.1, 1), ncol = 2, pt.size= 12) 
p5
dev.off()


pdf(paste0("/home/nasim/4dmt.abio/NGS087_3_outputs/","Spatial.008_top_markers_spatialfeature6.pdf"),width = 10,height = 10)
p6 <- SpatialFeaturePlot(NGS087_3, features =c("SEC31A", "PTPN1", "IGHD", "TXNDC5"),alpha = c(0.1, 1), ncol = 2, pt.size= 12) 
p6
dev.off()


pdf(paste0("/home/nasim/4dmt.abio/NGS087_3_outputs/","Spatial.008_top_markers_spatialfeature7.pdf"),width = 10,height = 10)
p7 <- SpatialFeaturePlot(NGS087_3, features =c("TNFRSF21", "FAT2", "PTPRS", "FOSB"),alpha = c(0.1, 1), ncol = 2, pt.size= 12) 
p7
dev.off()

pdf(paste0("/home/nasim/4dmt.abio/NGS087_3_outputs/","Spatial.008_top_markers_spatialfeature8.pdf"),width = 10,height = 10)
p8 <- SpatialFeaturePlot(NGS087_3, features =c( "IGHA1", "IGKC", "DPYSL2", "CLDN18"),alpha = c(0.1, 1), ncol = 2, pt.size= 12) 
p8
dev.off()


#Here the Bansky analysis starts to determine the anotomical regions
if (!requireNamespace("Banksy", quietly = TRUE)) {
  remotes::install_github("prabhakarlab/Banksy@devel")
}
library(SeuratWrappers)
library(Banksy)

object_Bansky <- RunBanksy(NGS087_3,
                           lambda = 0.8, verbose = TRUE,
                           assay = "Spatial.008um", slot = "data", features = "variable",
                           k_geom = 50
)

DefaultAssay(object_Bansky) <- "BANKSY"
object_Bansky<- RunPCA(object_Bansky, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object_Bansky), npcs = 20)
object_Bansky <- FindNeighbors(object_Bansky, reduction = "pca.banksy", dims = 1:20)
object_Bansky <- FindClusters(object_Bansky, cluster.name = "banksy_cluster", resolution = 0.5)

pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_SpatialDimPlot_tissue_domains_all_clusters.pdf", width = 20,height = 15)
Idents(object_Bansky) <- "banksy_cluster"
p <- SpatialDimPlot(object_Bansky, group.by = "banksy_cluster",  label = T, repel = T,label.size = 4, pt.size.factor = 12) +
  theme(
    legend.text = element_text(size = 15), # Increase legend item text size
    legend.title = element_text(size = 15) # Increase legend title size
  )
p
dev.off()

#make SpatialDim plot of some tissue regions
pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_SpatialDimPlot_tissue_domains_individual_clusters.pdf")
banksy_cells <- CellsByIdentities(object_Bansky, idents = c(0, 1, 2, 3, 4, 5))
p <- SpatialDimPlot(object_Bansky, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = T, pt.size.factor = 12) + NoLegend()


banksy_cells <- CellsByIdentities(object_Bansky, idents = c(6, 7, 8, 9, 10, 11))
p1 <- SpatialDimPlot(object_Bansky, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = T, pt.size.factor = 12) + NoLegend()

banksy_cells <- CellsByIdentities(object_Bansky, idents = c(12, 13, 14, 15))
p2 <- SpatialDimPlot(object_Bansky, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = T, pt.size.factor = 12) + NoLegend()

p
p1
p2

dev.off()

saveRDS(object_Bansky, '/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_Bansky.rds')


###Task 2: Cell Type Deconvolution via RCTD
if (!requireNamespace("spacexr", quietly = TRUE)) {
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
} else {
  # If the package is already installed, you can update it to the latest version
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE, force = TRUE)
}

library(spacexr)


#Read the h5.ad reference object and realted files (gene_names and cell_names), prepare seurat object
#and select the annotaion level of reference
library(Matrix)
library(Seurat)

# Read matrix, need to transverse the matrix, since the matrix coming from h5ad has cell_names as rows and gene_names as columns
expr <- readMM("/home/nasim/4dmt.abio/expr_matrix.mtx")
expr <- t(expr)

#gene_names and cell_names is a dataframe with one column, here get the contents of the only column, it will be saved as vector of characters
gene_names <- read.csv("/home/nasim/4dmt.abio/gene_names.csv", header = TRUE)[[1]]
length(gene_names)  # Should be 3000
rownames(expr) <- gene_names


cell_names <- read.csv("/home/nasim/4dmt.abio/cell_names.csv", header = TRUE)[[1]]
length(cell_names)  # Should be 50000
colnames(expr) <- cell_names

seurat_obj <- CreateSeuratObject(counts = expr)

# Add metadata
meta <- read.csv("/home/nasim/4dmt.abio/obs_metadata.csv", row.names = 1)
seurat_obj <- AddMetaData(seurat_obj, metadata = meta)
print(dim(seurat_obj))

#remove "Unknown Cells" 
seurat_obj <- subset(seurat_obj, subset=ann_level_2!="Unknown")
print(dim(seurat_obj))

#the var meta has gene names in different format, for now add the whole file to the seurat object
var_meta <- read.csv("/home/nasim/4dmt.abio/var_metadata.csv", row.names = 1)

# Make sure rownames match Seurat object features
var_meta <- var_meta[rownames(seurat_obj), , drop = FALSE]
seurat_obj[["RNA"]]@misc$features <- var_meta



#function does all the preprocessig steps for reference to make it ready for RCTD
prepare_RCTD_reference_from_seurat <- function(seurat_obj_subset,
                                               var_metadata_path,
                                               cell_type_col = "ann_level_4",
                                               min_cells = 25,
                                               umi_col = "nCount_RNA") {
  library(Seurat)
  library(Matrix)
  library(spacexr)
  
  # Step 1: Check column exists
  if (!cell_type_col %in% colnames(seurat_obj_subset@meta.data)) {
    stop(paste("Column", cell_type_col, "not found in Seurat metadata"))
  }
  
  # Step 2: Filter valid cell types
  cell_counts <- table(seurat_obj_subset@meta.data[[cell_type_col]])
  valid_types <- names(cell_counts[cell_counts >= min_cells])
  
  # Correct way to subset the Seurat object
  seurat_obj_subset <- subset(seurat_obj_subset, cells = rownames(seurat_obj_subset@meta.data)[seurat_obj_subset@meta.data[[cell_type_col]] %in% valid_types])
  Idents(seurat_obj_subset) <- cell_type_col
  
  # Step 3: Load var_metadata
  var_meta <- read.csv(var_metadata_path, row.names = 1)
  var_meta <- var_meta[rownames(seurat_obj_subset), , drop = FALSE]
  
  # Step 4: Extract and clean counts
  counts <- GetAssayData(seurat_obj_subset, assay = "RNA", slot = "counts", layer = "counts")
  counts@x <- round(counts@x)
  
  # Step 5: Rename genes with gene symbols
  gene_symbols <- var_meta$original_gene_symbols
  names(gene_symbols) <- rownames(var_meta)
  rownames(counts) <- gene_symbols[rownames(counts)]
  counts <- counts[!is.na(rownames(counts)) & !duplicated(rownames(counts)), ]
  
  # Step 6: Prepare cluster (factor of cell types)
  cluster <- as.factor(seurat_obj_subset@meta.data[colnames(counts), cell_type_col])
  names(cluster) <- colnames(counts)  # ✅ Add names!
  
  
  # Step 7: Prepare nUMI
  nUMI <- seurat_obj_subset@meta.data[colnames(counts), umi_col]
  nUMI <- as.integer(round(nUMI))
  names(nUMI) <- colnames(counts)
  
  # Step 8: Filter for valid (non-NA) entries
  valid_cells <- names(nUMI)[!is.na(nUMI)]
  counts <- counts[, valid_cells]
  cluster <- cluster[valid_cells]
  nUMI <- nUMI[valid_cells]
  
  # Final check
  stopifnot(
    sum(is.na(nUMI)) == 0,
    length(nUMI) == ncol(counts),
    length(cluster) == length(nUMI),
    all(names(nUMI) == colnames(counts))
  )
  
  # Step 9: Build Reference
  reference <- Reference(counts, cluster, nUMI)
  return(reference)
}




#Call the function which makes RCTD object for the reference
reference <- prepare_RCTD_reference_from_seurat(
  seurat_obj_subset = seurat_obj_subset,
  var_metadata_path = "/home/nasim/4dmt.abio/var_metadata.csv",
  cell_type_col = "ann_level_4",
  min_cells = 25
)

#preprocessing the query and make query RCTD object and run RCTD deconvolution

# localdir <- "/home/nasim/4dmt.abio/NGS087_3/outs/binned_outputs/square_008um/"
# NGS087_3 <- Load10X_Spatial(data.dir = localdir, file= "filtered_feature_bc_matrix.h5", assay="Spatial.008um")


localdir <- "/home/nasim/4dmt.abio/NGS087_4/outs/binned_outputs/square_008um/"
NGS087_4 <- Load10X_Spatial(data.dir = localdir, file= "filtered_feature_bc_matrix.h5", assay="Spatial.008um")

# DefaultAssay(NGS087_3) <- "Spatial.008um"
# NGS087_3  <- NormalizeData(NGS087_3)


DefaultAssay(NGS087_4) <- "Spatial.008um"
NGS087_4  <- NormalizeData(NGS087_4)


#sketch more cells for better results
# DefaultAssay(NGS087_3) <- "Spatial.008um"
# NGS087_3 <- FindVariableFeatures(NGS087_3)
# NGS087_3 <- ScaleData(NGS087_3)
# #default number of cells are 5000
# # we select 50,0000 cells for test and create a new 'sketch' assay
# NGS087_3 <- SketchData(
#   object = NGS087_3,
#   ncells = 50000,
#   method = "LeverageScore",
#   sketched.assay = "sketch"
# )


DefaultAssay(NGS087_4) <- "Spatial.008um"
NGS087_4 <- FindVariableFeatures(NGS087_4)
NGS087_4 <- ScaleData(NGS087_4)
#default number of cells are 5000
# we select 50,0000 cells for test and create a new 'sketch' assay
NGS087_4 <- SketchData(
  object = NGS087_4,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch analysis to sketched cells
#  DefaultAssay(NGS087_3) <- "sketch"
# # perform clustering workflow
#  NGS087_3 <- FindVariableFeatures(NGS087_3)
#  NGS087_3 <- ScaleData(NGS087_3)
#  NGS087_3 <- RunPCA(NGS087_3 , assay = "sketch", reduction.name = "pca.sketch")
# #pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Elbow_plot.pdf")
# P1 <- ElbowPlot(NGS087_3, reduction = "pca.sketch")
# print(P1)
# #dev.off()
# 
# NGS087_3 <- FindNeighbors(NGS087_3 , assay = "sketch", reduction = "pca.sketch", dims = 1:50)
# NGS087_3 <- FindClusters(NGS087_3 , cluster.name = "seurat_cluster.sketched", resolution = 1)
# NGS087_3 <- RunUMAP(NGS087_3, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)



# switch analysis to sketched cells
DefaultAssay(NGS087_4) <- "sketch"

# perform clustering workflow
NGS087_4 <- FindVariableFeatures(NGS087_4)
NGS087_4 <- ScaleData(NGS087_4)
NGS087_4 <- RunPCA(NGS087_4 , assay = "sketch", reduction.name = "pca.sketch")
#
NGS087_4 <- FindNeighbors(NGS087_4 , assay = "sketch", reduction = "pca.sketch", dims = 1:50)
NGS087_4 <- FindClusters(NGS087_4 , cluster.name = "seurat_cluster.sketched", resolution = 1)
NGS087_4 <- RunUMAP(NGS087_4, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)




#prepare the query for RCTD
# #NGS087_3 <- Spatial.008_Bansky
# DefaultAssay(NGS087_3) <- "sketch"
# counts_hd <- NGS087_3[["sketch"]]$counts
# cells_hd <- colnames(NGS087_3[["sketch"]])
# coords <- GetTissueCoordinates(NGS087_3)[cells_hd, 1:2]
# 
# query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))


DefaultAssay(NGS087_4) <- "sketch"
counts_hd <- NGS087_4[["sketch"]]$counts
cells_hd <- colnames(NGS087_4[["sketch"]])
coords <- GetTissueCoordinates(NGS087_4)[cells_hd, 1:2]

query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))


#run RCTD
myRCTD <- create.RCTD(query, reference, max_cores = 4, UMI_min_sigma = 50)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')


#save
# pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/SpatialRNA_nUMI_histogram1.pdf")
# P1 <- hist(myRCTD@spatialRNA@nUMI, breaks = 50, main = "UMI per Visium HD spot")
# print(P1)
# dev.off()


pdf("/home/nasim/4dmt.abio/NGS087_4_outputs/SpatialRNA_nUMI_histogram1.pdf")
P1 <- hist(myRCTD@spatialRNA@nUMI, breaks = 50, main = "UMI per Visium HD spot")
print(P1)
dev.off()


#NGS087_3 <- AddMetaData(NGS087_3, metadata = myRCTD@results$results_df)
NGS087_4 <- AddMetaData(NGS087_4, metadata = myRCTD@results$results_df)

# Increase the max size limit for globals
options(future.globals.maxSize = 1e9)  # 1 GB limit


#Project to the rest
# DefaultAssay(NGS087_3) <- "sketch"
# NGS087_3$first_type <- as.character(NGS087_3$first_type)
# NGS087_3$first_type[is.na(NGS087_3$first_type)] <- "Unknown"
# object <- ProjectData(
#   object = NGS087_3,
#   assay = "Spatial.008um",
#   full.reduction = "full.pca.sketch",
#   sketched.assay = "sketch",
#   sketched.reduction = "pca.sketch",
#   umap.model = "umap.sketch",
#   dims = 1:50,
#   refdata = list(full_first_type = "first_type")
# )

#Project to the rest
DefaultAssay(NGS087_4) <- "sketch"
NGS087_4$first_type <- as.character(NGS087_4$first_type)
NGS087_4$first_type[is.na(NGS087_4$first_type)] <- "Unknown"
object <- ProjectData(
  object = NGS087_4,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)


p1 <- SpatialDimPlot(object, group.by = "first_type", pt.size.factor = 10)
p2 <- SpatialDimPlot(object, group.by = "second_type", pt.size.factor = 10)
p3 <- SpatialDimPlot(object, group.by = "full_first_type", pt.size.factor = 10)
p1 | p2
p3

SpatialDimPlot(object, group.by = "full_first_type", label = TRUE)
DimPlot(object, group.by = "full_first_type", reduction = "full.umap.sketch", label = TRUE)
DimPlot(object, group.by = "full_first_type", reduction = "full.pca.sketch", label = TRUE)


#save 
pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Deconvolution_ann_level_4_first_type.pdf")
DefaultAssay(object) <- "Spatial.008um"
#Idents(object) <- "full_first_type"
Idents(object) <- "first_type"

#select the cells
cells <- CellsByIdentities(object)
#excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
p1 <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = F, pt.size.factor = 20)
p1
dev.off()

pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Deconvolution_ann_level_4_second_type.pdf")
DefaultAssay(object) <- "Spatial.008um"
#Idents(object) <- "full_first_type"
Idents(object) <- "second_type"

#selct the cells
cells <- CellsByIdentities(object)
#excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
p2 <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = F, pt.size.factor = 20)
p2
dev.off()

pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Deconvolution_ann_level_4_full_first_type.pdf")
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "full_first_type"

cells <- CellsByIdentities(object)
#excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
p3 <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = F, pt.size.factor = 20)
p3
dev.off()





pdf("/home/nasim/4dmt.abio/NGS087_4_outputs/Deconvolution_ann_level_4_first_type.pdf")
DefaultAssay(object) <- "Spatial.008um"
#Idents(object) <- "full_first_type"
Idents(object) <- "first_type"

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(object)
#excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
p1 <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = F, pt.size.factor = 20)
p1
dev.off()
# 
# 
# 
pdf("/home/nasim/4dmt.abio/NGS087_4_outputs/Deconvolution_ann_level_4_second_type.pdf")
DefaultAssay(object) <- "Spatial.008um"
#Idents(object) <- "full_first_type"
Idents(object) <- "second_type"

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(object)
#excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
p2 <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = F, pt.size.factor = 20)
p2
dev.off()
# 
# 
pdf("/home/nasim/4dmt.abio/NGS087_4_outputs/Deconvolution_ann_level_4_full_first_type.pdf")
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "full_first_type"


cells <- CellsByIdentities(object)
#excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
p3 <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FF0000", "grey50"), facet.highlight = T, combine = F, pt.size.factor = 20)
p3
dev.off()
# 

#cell numbers
write.csv(table(object@meta.data$first_type), '/home/nasim/4dmt.abio/NGS087_3_outputs/query_celltype_ann_level_4_first_type.csv', row.names = FALSE )

write.csv(table(object@meta.data$second_type), '/home/nasim/4dmt.abio/NGS087_3_outputs/query_celltype_ann_level_4_second_type.csv', row.names = FALSE )

write.csv(table(object@meta.data$full_first_type), '/home/nasim/4dmt.abio/NGS087_3_outputs/query_celltype_ann_level_4_full_first_type.csv', row.names = FALSE )

write.csv(table(seurat_obj_subset@meta.data$ann_level_4), '/home/nasim/4dmt.abio/NGS087_3_outputs/reference_ann_level_4.csv', row.names = FALSE )



write.csv(table(object@meta.data$first_type), '/home/nasim/4dmt.abio/NGS087_4_outputs/query_celltype_ann_level_4_first_type.csv', row.names = FALSE )

write.csv(table(object@meta.data$second_type), '/home/nasim/4dmt.abio/NGS087_4_outputs/query_celltype_ann_level_4_second_type.csv', row.names = FALSE )

write.csv(table(object@meta.data$full_first_type), '/home/nasim/4dmt.abio/NGS087_4_outputs/query_celltype_ann_level_4_full_first_type.csv', row.names = FALSE )

write.csv(table(seurat_obj_subset@meta.data$ann_level_4), '/home/nasim/4dmt.abio/NGS087_4_outputs/reference_ann_level_4.csv', row.names = FALSE )

#save
saveRDS(object, '/home/nasim/4dmt.abio/NGS087_3_outputs/Spatial.008_deconvolution.rds')
saveRDS(myRCTD, '/home/nasim/4dmt.abio/NGS087_3_outputs/RCTD_object.rds')
saveRDS(object, '/home/nasim/4dmt.abio/NGS087_4_outputs/Spatial.008_deconvolution.rds')
saveRDS(myRCTD, '/home/nasim/4dmt.abio/NGS087_4_outputs/RCTD_object.rds')

#extract weights and add to seurat object
# Extract the weights matrix (genes x spots)
weight_df <- as.data.frame(myRCTD@results$weights)

# Check structure
head(rownames(weight_df))  # should match colnames(object)
head(colnames(weight_df))  # should be cell types

# Make sure the spot barcodes align
shared_barcodes <- intersect(rownames(weight_df), colnames(object))

library(dplyr)
# Convert Seurat metadata to a data.frame with barcode as a column
meta_df <- object@meta.data %>% tibble::rownames_to_column("barcode")

# Convert weight_df to a data.frame and add barcodes as a column
weight_df2 <- weight_df %>% tibble::rownames_to_column("barcode")

# Join by barcode — retains all Seurat barcodes, fills missing weights with NA
meta_joined <- left_join(meta_df, weight_df2, by = "barcode")

rownames(meta_joined) <- meta_joined$barcode
#remove barcode column
object@meta.data <- meta_joined %>% select(-barcode)



#weight plots
pdf("/home/nasim/4dmt.abio/NGS087_3_outputs/Deconvolution_celltype_weights.pdf")
pdf("/home/nasim/4dmt.abio/NGS087_4_outputs/Deconvolution_celltype_weights.pdf")
object$Multiciliated[is.na(object$Multiciliated)] <- 0
p1 <- SpatialFeaturePlot(object, features = "Multiciliated" , pt.size= 20) + theme(legend.position = "right")

object$`SMG duct`[is.na(object$`SMG duct`)] <- 0 
p2 <- SpatialFeaturePlot(object, features = "SMG duct" , pt.size= 20) + theme(legend.position = "right")


object$`Basal resting`[is.na(object$`Basal resting`)] <- 0 
p3 <- SpatialFeaturePlot(object, features = "Basal resting" , pt.size= 20) + theme(legend.position = "right")

object$Goblet[is.na(object$Goblet)] <- 0 
p4 <- SpatialFeaturePlot(object, features = "Goblet" , pt.size= 20) + theme(legend.position = "right")


object$Suprabasal[is.na(object$Suprabasal)] <- 0 
p5 <- SpatialFeaturePlot(object, features = "Suprabasal" , pt.size= 20) + theme(legend.position = "right")

object$`Transitional Club-AT2`[is.na(object$`Transitional Club-AT2`)] <- 0 
p6 <- SpatialFeaturePlot(object, features = "Transitional Club-AT2" , pt.size= 20) + theme(legend.position = "right")


object$Club[is.na(object$Club)] <- 0 
p7 <- SpatialFeaturePlot(object, features = "Club" , pt.size= 20) + theme(legend.position = "right")

object$`SMG serous`[is.na(object$`SMG serous`)] <- 0 
p8 <- SpatialFeaturePlot(object, features = "SMG serous" , pt.size= 20) + theme(legend.position = "right")

object$`Hillock-like`[is.na(object$`Hillock-like`)] <- 0 
p9 <- SpatialFeaturePlot(object, features = "Hillock-like" , pt.size= 20) + theme(legend.position = "right")

object$Ionocyte[is.na(object$Ionocyte)] <- 0 
p10 <- SpatialFeaturePlot(object, features = "Ionocyte" , pt.size= 20) + theme(legend.position = "right")

object$`SMG mucous`[is.na(object$`SMG mucous`)] <- 0 
p11 <- SpatialFeaturePlot(object, features = "SMG mucous" , pt.size= 20) + theme(legend.position = "right")

object$Deuterosomal[is.na(object$Deuterosomal)] <- 0 
p12 <- SpatialFeaturePlot(object, features = "Deuterosomal" , pt.size= 20) + theme(legend.position = "right")



print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
print(p10)
print(p11)
print(p12)
dev.off()