# ------------------------------------------------------------------------------
#title: "Gene and Protein Expression Profiling from FFPE Tissues — Visium CytAssist"
#Author: Nasim Rahmatpour
#Date: "2025-03-31"
# ---

# Project Overview

#Spatial transcriptomics data generated with **10x Visium CytAssist** platform enables simultaneous profiling of **mRNA** and **protein** markers in FFPE tissue sections. This project aims to:
  
# - Perform separate analysis of spatial RNA and protein expression
# - Visualize clustering and spatial patterns of cellular compartments
# - Infer cell type composition in each spot by leveraging a matched scRNA-seq reference through **SPOTlight deconvolution**
  
---
  
###Task 1: Spatial RNA and Protein Analysis  

library(ggplot2)
library(Matrix)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
library(Signac)

#make seurat
counts <- Read10X_h5(filename=paste0("./filtered_feature_bc_matrix.h5"))
Spatial_GEX <- CreateSeuratObject(counts = counts$`Gene Expression`, project = "CytAssist, RNA_ADT") 
Spatial_GEX[['Protein']] = CreateAssayObject(counts = counts$`Antibody Capture`) 

#read image
# to make this work I removed one of the csv files in spatial folder "tissue_positions.csv"
image= Read10X_Image(
  image.dir = "./spatial/",
  #image.name = "tissue_lowres_image.png",
  filter.matrix = TRUE,
)
# Correct the image data to match the Seurat object
image <- image[Cells(Spatial_GEX)]
DefaultAssay(image) <- DefaultAssay(Spatial_GEX)
# Add the image to the object
Spatial_GEX[["slice"]] <- image

Spatial_GEX[["percent.mt"]] <- PercentageFeatureSet(object=Spatial_GEX, pattern = "^MT-")
Spatial_GEX[["percent.hb"]] <- PercentageFeatureSet(object=Spatial_GEX, pattern = "^HB")
Spatial_GEX[["percent.ribo"]] <- PercentageFeatureSet(object=Spatial_GEX, pattern = "^RP[SL]")

#these need to be changed before making SpatialFeaturePlot
Spatial_GEX@images[["slice"]]@coordinates[["tissue"]] <- as.integer(Spatial_GEX@images[["slice"]]@coordinates[["tissue"]])
Spatial_GEX@images[["slice"]]@coordinates[["row"]] <- as.integer(Spatial_GEX@images[["slice"]]@coordinates[["row"]])
Spatial_GEX@images[["slice"]]@coordinates[["col"]] <- as.integer(Spatial_GEX@images[["slice"]]@coordinates[["col"]])
Spatial_GEX@images[["slice"]]@coordinates[["imagerow"]] <- as.integer(Spatial_GEX@images[["slice"]]@coordinates[["imagerow"]])
Spatial_GEX@images[["slice"]]@coordinates[["imagecol"]] <- as.integer(Spatial_GEX@images[["slice"]]@coordinates[["imagecol"]])


#figures
dir.create("./GEX_PEX_out/")
pdf(paste0("./GEX_PEX_out/","features.pdf"), width=20, height=10)
p1 <- VlnPlot(Spatial_GEX, features = "nCount_RNA", pt.size = 0.1) + NoLegend() + theme(axis.text = element_text(size = 10)) + theme(plot.title = element_text(size = 20))
p2 <- VlnPlot(Spatial_GEX, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() + theme(axis.text = element_text(size = 10)) + theme(plot.title = element_text(size = 20))

p3 <- VlnPlot(Spatial_GEX, features = "nCount_Protein", pt.size = 0.1) + NoLegend() + theme(axis.text = element_text(size = 10)) + theme(plot.title = element_text(size = 20))
p4 <- VlnPlot(Spatial_GEX, features = "nFeature_Protein", pt.size = 0.1) + NoLegend() + theme(axis.text = element_text(size = 10)) + theme(plot.title = element_text(size = 20))



p5 <- SpatialFeaturePlot(Spatial_GEX, features = "nCount_RNA", pt.size = 5) + theme(legend.position = "right")+ theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size=20)) 
p6 <- SpatialFeaturePlot(Spatial_GEX, features = "nFeature_RNA", pt.size = 5) + theme(legend.position = "right") + theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size=20))



p7 <- SpatialFeaturePlot(Spatial_GEX, features = "nCount_Protein", pt.size = 5) + theme(legend.position = "right")+ theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size=20)) 
p8 <- SpatialFeaturePlot(Spatial_GEX, features = "nFeature_Protein", pt.size = 5) + theme(legend.position = "right") + theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size=20))


p9 <- SpatialFeaturePlot(Spatial_GEX, features = "percent.mt", pt.size = 5) 
p10 <- SpatialFeaturePlot(Spatial_GEX, features = "percent.hb", pt.size = 5) 
p11 <- SpatialFeaturePlot(Spatial_GEX, features = "percent.ribo", pt.size = 5) 
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
#wrap_plots(plot1, plot2)
#dev.off()

meta <- Spatial_GEX@meta.data
median.value <- apply(meta[,c('nCount_RNA','nFeature_RNA','nCount_Protein','nFeature_Protein','percent.mt','percent.hb',  'percent.ribo')],2,median)
write.csv(median.value, paste0("./GEX_PEX_out//median.csv"))

Spatial_GEX <- SCTransform(Spatial_GEX, assay = "RNA", verbose = FALSE)
Spatial_GEX <- RunPCA(Spatial_GEX, assay = "SCT", verbose = FALSE)
Spatial_GEX <- FindNeighbors(Spatial_GEX, reduction = "pca", dims = 1:30)
Spatial_GEX <- FindClusters(Spatial_GEX, verbose = FALSE)
Spatial_GEX<- RunUMAP(Spatial_GEX, reduction = "pca", dims = 1:30)

pdf(paste0("./GEX_PEX_out/","umap.pdf"), width=20, height=10)
p1 <- DimPlot(Spatial_GEX, reduction = "umap", label = TRUE, label.size = 5)
p2 <- SpatialDimPlot(Spatial_GEX, label = TRUE, label.size = 5, pt.size.factor=5)
#p3 <- SpatialDimPlot(Spatial_GEX, cells.highlight = CellsByIdentities(object = Spatial_GEX, idents = c(0, 1, 2, 3,
#4, 5, 6, 7, 8, 9, 10, 11, 12)), facet.highlight = TRUE, ncol = 1)
p3 <- SpatialDimPlot(Spatial_GEX, cells.highlight = CellsByIdentities(object = Spatial_GEX, idents = c(0,1)), facet.highlight = TRUE, ncol = 2, pt.size.factor=5)
p4 <- SpatialDimPlot(Spatial_GEX, cells.highlight = CellsByIdentities(object = Spatial_GEX, idents = c(2,3)), facet.highlight = TRUE, ncol = 2, pt.size.factor=5)
p5 <- SpatialDimPlot(Spatial_GEX, cells.highlight = CellsByIdentities(object = Spatial_GEX, idents = c(4,5)), facet.highlight = TRUE, ncol = 2, pt.size.factor=5)
p6 <- SpatialDimPlot(Spatial_GEX, cells.highlight = CellsByIdentities(object = Spatial_GEX, idents = c(6,7)), facet.highlight = TRUE, ncol = 2, pt.size.factor=5)
p7 <- SpatialDimPlot(Spatial_GEX, cells.highlight = CellsByIdentities(object = Spatial_GEX, idents = c(8,9)), facet.highlight = TRUE, ncol = 2, pt.size.factor=5)
p8 <- SpatialDimPlot(Spatial_GEX, cells.highlight = CellsByIdentities(object = Spatial_GEX, idents = c(10,11)), facet.highlight = TRUE, ncol = 2, pt.size.factor=5)
p9 <- SpatialDimPlot(Spatial_GEX, cells.highlight = CellsByIdentities(object = Spatial_GEX, idents = c(12)), facet.highlight = TRUE, ncol = 1, pt.size.factor=5)
print(p1 + p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
dev.off()

#finding the cluster markers
markers <- FindAllMarkers(Spatial_GEX, only.pos = TRUE)
write.csv(markers, file =paste0("./GEX_PEX_out/", "cluster_markers.csv"))

#save the top20 cluster markers
top20_markers <- markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top20<- c()
for (c in unique(top20_markers$cluster)){
  tmp <- top20_markers[top20_markers$cluster==c,]$gene
  top20 <- cbind(top20, tmp)
}
colnames(top20) <- paste0("C_", unique(top20_markers$cluster))
write.csv(top20, paste0("./GEX_PEX_out/", "top20_cluster_markers.csv"))

#finding the variable features, spatial heterogenity
Spatial_GEX<- FindSpatiallyVariableFeatures(Spatial_GEX, assay = "SCT", features = VariableFeatures(Spatial_GEX)[1:1000], selection.method = "markvariogram")

pdf(paste0("./GEX_PEX_out/","Spatial_heterogenity_genes.pdf"))
SpatialFeaturePlot(Spatial_GEX , features = top.features, pt.size = 5, ncol = 3, alpha = c(0.1, 1))
dev.off()

saveRDS(Spatial_GEX, "./GEX_PEX_out/Sparial_GEX.rds" )

#expected gene markers
pdf(paste0("./GEX_PEX_out/","expected_genemarkers_RNAs.pdf"))
P1=SpatialFeaturePlot(Spatial_GEX , features = c("CD4","CD3D", "FOXP3"), pt.size = 5, ncol = 3)
P2=SpatialFeaturePlot(Spatial_GEX , features = c("CD8A", "CD8B", "NCAM1"), pt.size = 5, ncol = 2)
P3= SpatialFeaturePlot(Spatial_GEX , features = c("CD19","CD79A", "MS4A1","BLK", "MZBI", "JCHAIN"), pt.size = 5, ncol = 3)
P4= SpatialFeaturePlot(Spatial_GEX , features = c("S100A9","S100A8","CSF3R","CSF1R","LILRA4", "TPSAB1", "KIT"), pt.size = 5, ncol = 3)
P5= SpatialFeaturePlot(Spatial_GEX , features = c("CD14","CD68", "MS4A6A"), pt.size = 5, ncol = 3)
P6= SpatialFeaturePlot(Spatial_GEX , features = c("COL11A1","COL1A2", "EPCAM", "PECAM1"), pt.size = 5, ncol = 3)
print(P1)
print(P2)
print(P3)
print(P4)
print(P5)
print(P6)
dev.off()

#Protein data (ADT data) 
pdf(paste0("./GEX_PEX_out/","expected_genemarkers_Proteins_Raw.pdf"))
DefaultAssay(Spatial_GEX) <- "Protein"
P1=SpatialFeaturePlot(Spatial_GEX , features = c("CD163.1", "CR2.1", "PCNA.1"), pt.size = 5, ncol = 3)
P2=SpatialFeaturePlot(Spatial_GEX , features = c("VIM.1", "KRT5.1", "CD68.1"), pt.size = 5, ncol = 3)
P3=SpatialFeaturePlot(Spatial_GEX , features = c("CEACAM8.1", "PTPRC.1", "HLA-DRA"), pt.size = 5, ncol = 3)
P4=SpatialFeaturePlot(Spatial_GEX , features = c("PAX5.1", "SDC1.1", "PTPRC.2"), pt.size = 5, ncol = 3)
P5=SpatialFeaturePlot(Spatial_GEX , features = c("CD8A.1", "BCL2.1"), pt.size = 5, ncol = 3)
P6=SpatialFeaturePlot(Spatial_GEX , features = c("mouse-IgG2a", "mouse-IgG1k", "mouse-IgG2bk", "rat-IgG2a"), pt.size = 5, ncol = 3)
P7=SpatialFeaturePlot(Spatial_GEX , features = c("CD19.1", "PDCD1.1", "ACTA2.1"), pt.size = 5, ncol = 3)
P8=SpatialFeaturePlot(Spatial_GEX , features = c("FCGR3A.1", "ITGAX.1", "CXCR5.1"), pt.size = 5, ncol = 3)
P9=SpatialFeaturePlot(Spatial_GEX , features = c("EPCAM.1", "MS4A1.1", "CD3E.1"), pt.size = 5, ncol = 3)
P10=SpatialFeaturePlot(Spatial_GEX , features = c("CD14.1", "CD40.1", "PECAM1.1", "CD4.1"), pt.size = 5, ncol = 3)
print(P1)
print(P2)
print(P3)
print(P4)
print(P5)
print(P6)
print(P7)
print(P8)
print(P9)
print(P10)
dev.off()

Spatial_GEX <- NormalizeData(Spatial_GEX, normalization.method = "CLR", margin = 2)

pdf(paste0("./GEX_PEX_out/","expected_genemarkers_Proteins_Normalized.pdf"))
DefaultAssay(Spatial_GEX) <- "Protein"
P1=SpatialFeaturePlot(Spatial_GEX , features = c("CD163.1", "CR2.1", "PCNA.1"), pt.size = 5, ncol = 3)
P2=SpatialFeaturePlot(Spatial_GEX , features = c("VIM.1", "KRT5.1", "CD68.1"), pt.size = 5, ncol = 3)
P3=SpatialFeaturePlot(Spatial_GEX , features = c("CEACAM8.1", "PTPRC.1", "HLA-DRA"), pt.size = 5, ncol = 3)
P4=SpatialFeaturePlot(Spatial_GEX , features = c("PAX5.1", "SDC1.1", "PTPRC.2"), pt.size = 5, ncol = 3)
P5=SpatialFeaturePlot(Spatial_GEX , features = c("CD8A.1", "BCL2.1"), pt.size = 5, ncol = 3)
P6=SpatialFeaturePlot(Spatial_GEX , features = c("mouse-IgG2a", "mouse-IgG1k", "mouse-IgG2bk", "rat-IgG2a"), pt.size = 5, ncol = 3)
P7=SpatialFeaturePlot(Spatial_GEX , features = c("CD19.1", "PDCD1.1", "ACTA2.1"), pt.size = 5, ncol = 3)
P8=SpatialFeaturePlot(Spatial_GEX , features = c("FCGR3A.1", "ITGAX.1", "CXCR5.1"), pt.size = 5, ncol = 3)
P9=SpatialFeaturePlot(Spatial_GEX , features = c("EPCAM.1", "MS4A1.1", "CD3E.1"), pt.size = 5, ncol = 3)
P10=SpatialFeaturePlot(Spatial_GEX , features = c("CD14.1", "CD40.1", "PECAM1.1", "CD4.1"), pt.size = 5, ncol = 3)
print(P1)
print(P2)
print(P3)
print(P4)
print(P5)
print(P6)
print(P7)
print(P8)
print(P9)
print(P10)
dev.off()


#Task 2: Cell Type Deconvolution with SPOTlight

library(Seurat)
library(ggplot2)
library(dplyr)
library(purrr)
library(SPOTlight)
library(NMF)
library(nnls)
library(scrattch.io)
library(ggplot2)
library(scater)
library(scran)
library(SingleCellExperiment)
library(readr)
library(tidyverse)
library(Matrix)
library(irlba)
library(bigstatsr)
library(RSpectra)

#set the common parameters
clust_vr <- "final_celltype"
hvg <- 3000
FC <- 0.5
pct1 <- 0.9
date <- Sys.Date()
sathe_meta <- read_csv("./seurat_analysis/cell_labels.csv")
sathe_cells <- split(sathe_meta, sathe_meta$orig.ident)
sathe_matrices <- lapply(names(sathe_cells), function(x){
  print(x)
  mat <- Read10X(file.path('./seurat_analysis/sathe_matrix', x))
  cells <- substr(sathe_cells[[x]]$cell_barcode, 1, 16)
  colnames(mat) <- substr(colnames(mat), 1, 16)
  mat <- mat[, cells]
  colnames(mat) <- sathe_cells[[x]]$cell_barcode
  mat
})
sathe_matrix <- Reduce(cbind, sathe_matrices)
sathe_metadata <- bind_rows(sathe_cells) %>% column_to_rownames("cell_barcode")

sathedata <- CreateSeuratObject(sathe_matrix, meta.data = sathe_metadata)
#sathedata <- SCTransform(sathedata, assay = "RNA", verbose = FALSE)
sathedata <- NormalizeData(sathedata, verbose = FALSE)
sathedata <- FindVariableFeatures(sathedata, selection.method = "vst", nfeatures = 3000)
sathedata <- ScaleData(sathedata, verbose = FALSE)

sathedata <- RunPCA(sathedata, npcs = 30, verbose = FALSE)
sathedata <- RunUMAP(sathedata, reduction = "pca", dims = 1:30, verbose = TRUE)

#make the seurat 
sathe_meta <- read_csv("./seurat_analysis/cell_labels.csv")
sathe_cells <- split(sathe_meta, sathe_meta$orig.ident)
ids <- names(sathe_cells)

for (file in ids){
  prefix <- paste0("Data_",file)
  object <- Read10X(file.path('./seurat_analysis/sathe_matrix', file))
  objectdata <- CreateSeuratObject(counts = object, project = file) 
  assign(prefix, objectdata)
}

other <- ls(pattern="Data_")[-1]
m <- function(x){eval(as.name(x))}
data <- merge(x = get(ls(pattern="Data_")[1]),
              y = c(lapply(X = other, function(x) m(x))))


#get the gene markers
Seurat::Idents(object = sathedata) <- sathedata@meta.data[, clust_vr]
cluster_markers_all <- Seurat::FindAllMarkers(object = sathedata, 
                                              verbose = TRUE, 
                                              only.pos = TRUE,
                                              slot = "data", 
                                              min.pct = 0.9,
                                              max.cells.per.ident = 100)

#look at some quality of the gene markers
table(cluster_markers_all$cluster)
min(cluster_markers_all$avg_log2FC)
min(cluster_markers_all$pct.1)
min(cluster_markers_all$pct.2)

#save the single cell seurat file
saveRDS(object = sathedata, file = "./GEX_PEX_out/sathedata.RDS")
sathedata <- readRDS(file = "./GEX_PEX_out/sathedata.RDS")

#save the gene markers file
saveRDS(object = cluster_markers_all, file = "./GEX_PEX_out/markers_ref_majorcells.RDS")
cluster_markers_all <- readRDS(file = "./GEX_PEX_out/markers_ref_majorcells.RDS")

#filter gene markers
cluster_markers_filt <- cluster_markers_all %>% 
  filter(avg_log2FC > FC & pct.1 > pct1)
cluster_markers_filt %>% 
  dplyr::count(cluster) %>% 
  data.frame()

table(unique(cluster_markers_all$cluster) %in% unique(cluster_markers_filt$cluster))

#convert the single cell seurat object to single cell object
sce <- as.SingleCellExperiment(sathedata)
colLabels(sce) <- colData(sce)$final_celltype

#downsample the reference, 
#split is like making a dictionary in list format; assigning values to keys. sample does the reverse
# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$final_celltype)
# downsample to at most 100 per identity & subset
n_cells <- 100
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]
dim(sce)

#deconvolution
res <- SPOTlight(
  x = sce,
  y = Spatial_GEX@assays$SCT@counts,
  groups = sce$final_celltype,
  mgs = cluster_markers_filt,
  hvg = sathedata@assays$RNA@var.features,
  weight_id = "avg_log2FC",
  group_id = "cluster",
  gene_id = "gene")

#extract deconvolution matrix
head(mat <- res$mat)[, seq_len(3)]
colnames(mat) <- unique(sce$final_celltype)
# Extract NMF model fit
mod <- res$NMF

plotTopicProfiles(
  x = mod,
  y =(sce$final_celltype),
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

#how the individual topic profiles of each cell within each cell-type behave
plotTopicProfiles(
  x = mod,
  y = sce$final_celltype,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 9)

#see which genes are learned for each topic, higher number shows each gene is more relevant to each topic.
library(NMF)
sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)

library(ggcorrplot)
plotCorrelationMatrix(mat)

#colocalization
plotInteractions(mat, which = "heatmap", metric = "prop")
plotInteractions(mat, which = "heatmap", metric = "jaccard")
plotInteractions(mat, which = "network")

str(mat)
typeof(mat)
class(mat)

pdf("./GEX_PEX_out/Deconvolution.pdf")
ct <- colnames(mat)
mat[mat < 0.1] <- 0
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
  "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

plotSpatialScatterpie(
  x= Spatial_GEX,
  #y = mat.new[1,1],
  y = mat,
  cell_types = colnames(mat),
  #cell_types = as.vector(colnames(mat)),
  #cell_types=y[,c(1,2)],
  img = FALSE,
  #Coord_flip(),
  #slice="/home/nasim/Visium_737/spatial/tissue_lowres_image.png",
  #slice=slice1,
  scatterpie_alpha = 1,
  pie_scale = 0.3) + 
  scale_y_reverse()+
  
  scale_fill_manual(
    values = pal,
    breaks = names(pal))


dev.off()
#png("home/nasim/nasim.png")
#im = readImage("/home/nasim/test.png")
#flipImage(im, mode = 'vertical')
#dev.off()

Spatial_GEX@meta.data <- cbind(Spatial_GEX@meta.data, mat)

#celltype scores
pdf("./GEX_PEX_out/celltypes_scores.pdf")
SpatialFeaturePlot(Spatial_GEX, features = "B", pt.size = 5, image.alpha =  1)
SpatialFeaturePlot(Spatial_GEX, features = "CD4", pt.size = 5, image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "CD8", pt.size = 5,image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "DC", pt.size = 5, image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "endothelial", pt.size = 5, image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "epithelial", pt.size = 5, image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "fibroblasts", pt.size = 5, image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "macrophage", pt.size = 5, image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "mast", pt.size = 5, image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "NK", pt.size = 5, image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "pericytes",pt.size = 5,  image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "plasma", pt.size = 5, image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = "Treg", pt.size = 5, image.alpha = 1)
SpatialFeaturePlot(Spatial_GEX, features = c("B","CD4","CD8", "DC", "endothelial",  "epithelial", "macrophage", "mast", "NK", "pericytes", "plasma", "Treg") , pt.size = 5, image.alpha = 1, ncol=4) + theme(legend.text = element_text(size = 10)) + theme(legend.title = element_text(size=10)) 
dev.off()
mat.new <- as.data.frame(mat)
