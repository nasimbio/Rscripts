# ------------------------------------------------------------------------------
#title: "10x Visium Spatial Transcriptomics — Breast Cancer"
#Author: Nasim Rahmatpour
#Date: "2/25/2022"
# ------------------------------------------------------------------------------


# Project Overview

#Starting from a slice of human breast cancer tissue, spatial transcriptomics was performed using the 10x Genomics Visium platform. The goal was to uncover spatial gene expression patterns and identify the regional distribution and co-localization of cell types across the tissue section.
#A matched single-cell RNA-seq reference from a breast cancer atlas (Wu et al., 2021) was used to deconvolute spatial gene expression and estimate cell type proportions for each spatial spot.

---

  
###Task 1: Preprocessing and Spatial Clustering    
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

ids <- c("423", "625", "599", "737")

#create directories
setwd("/home/nasim")
for (i in ids){
  dir.create(paste0("~/Visium_", i))
}

#path to the images
image_paths <- c("~/Visium_423/spatial",
                 "~/Visium_625/spatial",
                 "~/Visium_599/spatial",
                 "~/Visium_737/spatial")

#read the images through loop and save all in one list
all_images <- list()
for (i in 1:length(ids)){
  all_images[[i]] <- Read10X_Image(
    #or use
    #image.dir=image_paths[i],
    image.dir= paste0("~/Visium_", ids[[i]],"/spatial"),
    image.name = "tissue_lowres_image.png",
    filter.matrix = TRUE,
  )
}

#read the images in the loop and save each separately by assign
for (i in ids){
  image = Read10X_Image(
    #image_path[i] does not work here because i in this loop is not a number/index
    #image.dir=image_paths[i],
    image.dir= paste0("~/Visium_",i,"/spatial"),
    image.name = "tissue_lowres_image.png",
    filter.matrix = TRUE,
  )
  assign(i,image)
}

#read the images in function
read_images <- function(x){
  x <- Read10X_Image(
    image.dir= paste0(x,"/spatial"),
    image.name = "tissue_lowres_image.png",
    filter.matrix = TRUE,
  )
}

#call the function
image_A <- read_images("Visium_423")
image_B <- read_images("Visium_625")
image_C <- read_images("Visium_599")
image_D <- read_images("Visium_737")

#path to the matrix files
matrix_path <- c("~/Visium_423",
                 "~/Visium_625",
                 "~/Visium_599",
                 "~/Visium_737")

#save all seurats in a list, add metadata to each at the same time
all_seurats <- list()
for (i in 1:length(ids)){
  all_seurats[[i]] <- Load10X_Spatial(
    #or use
    #data.dir = paste0("~/Visium_", ids[[i]]),
    data.dir = matrix_path[i],
    filename = "filtered_feature_bc_matrix.h5",
    assay= "Spatial",
    #slice= "D1",
    filter.matrix = TRUE,
    to.upper = FALSE,
    image=NULL)
  all_seurats[[i]][["percent.hb"]] <- PercentageFeatureSet(object=all_seurats[[i]], pattern="^HB")
  all_seurats[[i]][["percent.mt"]] <- PercentageFeatureSet(object=all_seurats[[i]], pattern="^MT-")
  
}

#save each seurat separately, add the metadata at the same time
for (i in ids){
  seurat = Load10X_Spatial(
    #the following does not work because i is not a number 
    #data.dir =matrix_path[i],
    data.dir = paste0("~/Visium_", i),
    filename = "filtered_feature_bc_matrix.h5",
    assay= "Spatial",
    #slice= "A1",
    filter.matrix = TRUE,
    to.upper = FALSE,
    image=NULL)
  seurat[["percent.hb"]] <- PercentageFeatureSet(object=seurat,pattern = "^HB")
  seurat[["percent.mt"]] <- PercentageFeatureSet(object=seurat, pattern = "^MT-")
  assign(i, seurat)
  
}

#read seurats in a function
read_seurats <- function(x){
  x <- Load10X_Spatial(
    data.dir = paste0("~/", x),
    filename = "filtered_feature_bc_matrix.h5",
    assay= "Spatial",
    #slice= "A1",
    filter.matrix = TRUE,
    to.upper = FALSE,
    image=NULL)
  #if you don't add the return(x), it will make a dataframe with one column variable
  x[["percent.hb"]] <- PercentageFeatureSet(object=x, pattern = "^HB")
  x[["percent.mt"]] <- PercentageFeatureSet(object=x, pattern="^MT-")
  x[["percent.mt"]] <- PercentageFeatureSet(object=x, pattern="^RP[SL]")
  return(x)
  
}


#call the function
Visium_423 <- read_seurats("Visium_423")
Visium_625 <- read_seurats("Visium_625")
Visium_599 <- read_seurats("Visium_599")
Visium_737 <- read_seurats("Visium_737")


#function to add columns
add_meta <- function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(object=x, pattern = "^MT-")
  x[["percent.hb"]] <- PercentageFeatureSet(object=x, pattern = "^HB")
  x[["percent.ribo"]] <- PercentageFeatureSet(object=x, pattern = "^RP[SL]")
  return(x)
}

#call the function
Visium_423 <- add_meta(Visium_423)
Visium_625 <- add_meta(Visium_625)
Visium_599 <- add_meta(Visium_599)
Visium_737 <- add_meta(Visium_737)

#also can apply the defined function on you list elements
all_seurats <- lapply(all_seurats, add_meta)

#or
all_seurats <- lapply(all_seurats, FUN = function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(object=x, pattern = "^MT-")
  x[["percent.hb"]] <- PercentageFeatureSet(object=x, pattern = "^HB")
  x[["percent.ribo"]] <- PercentageFeatureSet(object=x, pattern = "^RP[SL]")
  return(x)
})


#figures
pdf(paste0("~/Visium_737/","p1_p2_p3_4_p5_p6_p7.pdf"), width=20, height=10)
p1 <- VlnPlot(Visium_737, features = "nCount_Spatial", pt.size = 0.1) + NoLegend() + theme(axis.text = element_text(size = 20)) + theme(plot.title = element_text(size = 20))
p2 <- VlnPlot(Visium_737, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend() + theme(axis.text = element_text(size = 20)) + theme(plot.title = element_text(size = 20))
p3 <- SpatialFeaturePlot(Visium_737, features = "nCount_Spatial", pt.size = 2) + theme(legend.position = "right")+ theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size=20)) 
p4 <- SpatialFeaturePlot(Visium_737, features = "nFeature_Spatial", pt.size = 2) + theme(legend.position = "right") + theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size=20))
p5 <- SpatialFeaturePlot(Visium_737, features = "percent.mt", pt.size = 2) 
p6 <- SpatialFeaturePlot(Visium_737, features = "percent.hb", pt.size = 2) 
p7 <- SpatialFeaturePlot(Visium_737, features = "percent.ribo", pt.size = 2) 
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
#wrap_plots(plot1, plot2)
dev.off()



#get the median in function
get_median <- function(x){
  metax <- x@meta.data
  print(dim(metax))
  median.value <- apply(metax[,c('nCount_Spatial','nFeature_Spatial','percent.mt','percent.hb','percent.ribo')], 2, median)
  return(median.value)
  #the following inside the function does not work, do it when calling the function
  #write.csv(median.value, paste0("~/", x, "median_f.csv"))
}

#call the function
write.csv(get_median(`Visium_423`), paste0("~/Visium_423/median_f.csv"))
write.csv(get_median(`Visium_625`), paste0("~/Visium_625/median_f.csv"))
write.csv(get_median(`Visium_599`), paste0("~/Visium_599/median_f.csv"))
write.csv(get_median(`Visium_737`), paste0("~/Visium_737/median_f.csv"))

#save the seurat objects
saveRDS(Visium_423, "~/Visium_423/Visium_423_v1.rds")
saveRDS(Visium_625, "~/Visium_625/Visium_625_v1.rds")
saveRDS(Visium_599, "~/Visium_599/Visium_599_v1.rds")
saveRDS(Visium_737, "~/Visium_737/Visium_737_v1.rds")

#function for standard analysis
standard_analysis <- function(x){
  x <- SCTransform(x, assay = "Spatial", verbose = FALSE)
  x <- RunPCA(x, assay = "SCT", verbose = FALSE)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:30)
  x <- FindClusters(x, verbose = FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:30)
  return(x)
}

#call the function
Visium_423_v2 <- standard_analysis(`Visium_423`)
Visium_625_v2 <- standard_analysis(`Visium_625`)
Visium_599_v2 <- standard_analysis(`Visium_599`)
Visium_737_v2 <- standard_analysis(`Visium_737`)

#you can also lapply the defined function on the list elements
all_seurats <- lapply(all_seurats, standard_analysis)

#or
all_seurats <- lapply(all_seurats, FUN = function(x){
  x <- SCTransform(x, assay = "Spatial", verbose = FALSE)
  x <- RunPCA(x, assay = "SCT", verbose = FALSE)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:30)
  x <- FindClusters(x, verbose = FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:30)
})  

#make the umap}
pdf(paste0("~/Visium_423/","umap.pdf"), width=20, height=10)
p1 <- DimPlot(Visium_423_v2, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(Visium_423_v2, label = TRUE, label.size = 3)
p3 <- SpatialDimPlot(Visium_423_v2, cells.highlight = CellsByIdentities(object = Visium_423_v2, idents = c(0, 1, 2, 3,
                                                                                                           4, 5, 6, 7, 8, 9)), facet.highlight = TRUE, ncol = 3)
print(p1 + p2)
print(p3)
dev.off()

pdf(paste0("~/Visium_737/","umap.pdf"), width=20, height=10)
p1 <- DimPlot(Visium_737_v2, reduction = "umap", label = TRUE, label.size = 5)
p2 <- SpatialDimPlot(Visium_737_v2, label = TRUE, label.size = 5)
p3 <- SpatialDimPlot(Visium_737_v2, cells.highlight = CellsByIdentities(object = Visium_737_v2, idents = c(0, 1, 2, 3,
                                                                                                           4, 5, 6, 7, 8, 9, 10, 11, 12, 13)), facet.highlight = TRUE, ncol = 3)
print(p1 + p2)
print(p3)
dev.off()

SpatialDimPlot(Visium_423_v2, interactive = TRUE)
SpatialDimPlot(Visium_625_v2, interactive = TRUE)
SpatialDimPlot(Visium_599_v2, interactive = TRUE)
SpatialDimPlot(Visium_737_v2, interactive = TRUE)

#find specific genes
SpatialFeaturePlot(Visium_423_v2, features = "CD8A", interactive = TRUE)
SpatialFeaturePlot(Visium_625_v2, features = "CD8A", interactive = TRUE)
SpatialFeaturePlot(Visium_599_v2, features = "CD8A", interactive = TRUE)
SpatialFeaturePlot(Visium_737_v2, features = "CD8A", interactive = TRUE)

#this is usefull, showing the umap and the image side by side and you can highlight a specific cluster
LinkedDimPlot(Visium_423_v2)
LinkedDimPlot(Visium_625_v2)
LinkedDimPlot(Visium_599_v2)
LinkedDimPlot(Visium_737_v2)

#save 
saveRDS(Visium_423_v2, "~/Visium_423/Visium_423_v2.rds")
saveRDS(Visium_625_v2, "~/Visium_625/Visium_625_v2.rds")
saveRDS(Visium_599_v2, "~/Visium_599/Visium_599_v2.rds")
saveRDS(Visium_737_v2, "~/Visium_737/Visium_737_v2.rds")


#finding the variable features, spatial heterogenity
Visium_737_v2 <- FindSpatiallyVariableFeatures(Visium_737_v2, assay = "SCT", features = VariableFeatures(Visium_737_v2)[1:1000], selection.method = "markvariogram")

#Save them
variable.features_737 <- SpatiallyVariableFeatures(Visium_737_v2,selection.method = "markvariogram")
write.csv(variable.features_737, paste0("~/Visium_737/", "Spatial_markers_737.csv"),row.names = F,col.names = F)
top.features <- head(SpatiallyVariableFeatures(Visium_737_v2, selection.method = "markvariogram"), 10)
top10 <- as.data.frame(top.features[1:10])
write.csv(top10, paste0("~/Visium_737/", "top10_spatial_markers_737.csv"),row.names = F,col.names = F)

#visiualize the genes
SpatialFeaturePlot(Visium_737_v2 , features = top.features, ncol = 5, alpha = c(0.1, 1))


###Task 2: Spatial Deconvolution with SPOTlight

library(Seurat)
library(ggplot2)
library(dplyr)
library(purrr)
library(SPOTlight)
library(NMF)
library(nnls)
library(scrattch.io)
library(scater)
library(scran)
library(SingleCellExperiment)

#set the common parameters
clust_vr <- "celltype_major"
hvg <- 3000
FC <- 1
pct1 <- 0.9
date <- Sys.Date()

#set the data directory
setwd("/home/nasim")
dir.create(paste0("~/BRCA_scRNASeq"))
data_dir <- ("/home/nasim/BRCA_scRNASeq")

#read the 10X data, gen.column=1 is needed since the gene matrix has only one column
ref_data <- Read10X(data.dir = data_dir, gene.column=1)

#read the metadta
metadata <- read.csv(file ="metadata.csv", row.names = 1)

#filter the metadata
metadata1 <- subset(metadata, celltype_subset!="T_cells_c11_MKI67" & celltype_subset!="Cycling_Myeloid" & celltype_subset!="Cycling PVL")

#make seurat object
breast_ref <- Seurat::CreateSeuratObject(counts = ref_data[, rownames(metadata1)], 
                                         project = "visium_breast", 
                                         assay = "RNA", 
                                         meta.data = metadata1)
#normalize
breast_ref <- NormalizeData(breast_ref)

#get highly variable genes
breast_ref  <- FindVariableFeatures(breast_ref , selection.method = "vst", nfeatures = 3000)

#scale it
breast_ref <- ScaleData(breast_ref)

#get the gene markers
Seurat::Idents(object = breast_ref) <- breast_ref@meta.data[, clust_vr]
cluster_markers_all <- Seurat::FindAllMarkers(object = breast_ref, 
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
saveRDS(object = breast_ref, file = "/home/nasim/Visium_737/breast_ref_majorcells.RDS")
breast_ref <- readRDS(file = "/home/nasim/Visium_737/breast_ref_majorcells.RDS")

#save the gene markers file
saveRDS(object = cluster_markers_all, file = "/home/nasim/Visium_737/markers_breast_ref_majorcells.RDS")
cluster_markers_all <- readRDS(file = "/home/nasim/Visium_737/markers_breast_ref_majorcells.RDS")

#filter gene markers
cluster_markers_filt <- cluster_markers_all %>% 
  filter(avg_log2FC > FC & pct.1 > pct1)
cluster_markers_filt %>% 
  count(cluster) %>% 
  data.frame()

table(unique(cluster_markers_all$cluster) %in% unique(cluster_markers_filt$cluster))

#convert the single cell seurat object to single cell object
sce <- as.SingleCellExperiment(breast_ref)
colLabels(sce) <- colData(sce)$celltype_major

#downsample the reference, 
#split is like making a dictionary in list format; assigning values to keys. sample does the reverse

# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$celltype_major)
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
  y = Visium_737@assays$Spatial@counts,
  groups = sce$celltype_major,
  mgs = cluster_markers_filt,
  hvg = breast_ref@assays$RNA@var.features,
  weight_id = "avg_log2FC",
  group_id = "cluster",
  gene_id = "gene")

#extract deconvolution matrix
head(mat <- res$mat)[, seq_len(3)]
colnames(mat) <- unique(sce$celltype_major)
# Extract NMF model fit
mod <- res$NMF

plotTopicProfiles(
  x = mod,
  y =(sce$celltype_major),
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

#how the individual topic profiles of each cell within each cell-type behave
plotTopicProfiles(
  x = mod,
  y = sce$celltype_major,
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

#later make imagecol to negative to get right direction of the image
#pdf("/home/nasim/Visium_737/Deconvolution.pdf")
ct <- colnames(mat)
mat[mat < 0.1] <- 0
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
  "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

plotSpatialScatterpie(
  x= Visium_737_v2,
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
  pie_scale = 0.4) + 
  scale_y_reverse()+
  
  scale_fill_manual(
    values = pal,
    breaks = names(pal))


#dev.off()

Visium_737_v2@meta.data <- cbind(Visium_737_v2@meta.data, mat)

#celltype scores
pdf("/home/nasim/Visium_737/celltypes_scores.pdf")
SpatialFeaturePlot(Visium_737_v2, features = "B-cells", alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = "B-cells", image.alpha = 1)
SpatialFeaturePlot(Visium_737_v2, features = "CAFs", image.alpha = 1)
SpatialFeaturePlot(Visium_737_v2, features = "Cancer Epithelial", image.alpha = 1)
SpatialFeaturePlot(Visium_737_v2, features = "Endothelial", image.alpha = 1)
SpatialFeaturePlot(Visium_737_v2, features = "Myeloid", image.alpha = 1)
SpatialFeaturePlot(Visium_737_v2, features = "Normal Epithelial", image.alpha = 1)
SpatialFeaturePlot(Visium_737_v2, features = "Plasmablasts", image.alpha = 1)
SpatialFeaturePlot(Visium_737_v2, features = "PVL", image.alpha = 1)
SpatialFeaturePlot(Visium_737_v2, features = "T-cells", image.alpha = 1)
SpatialFeaturePlot(Visium_737_v2, features = c("B-cells","CAFs","Cancer Epithelial", "Endothelial", "Myeloid",  "Normal Epithelial", "Plasmablasts", "PVL", "T-cells") , image.alpha = 1, ncol=4) + theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size=20)) 
dev.off()

#gene markers expressions
#p <- SpatialFeaturePlot(Visium_737_v2, features = 'BANK1', alpha = c(0.1, 1))
#p <- SpatialFeaturePlot(Visium_737_v2, features = 'BANK1', image.alpha = 1)
#p$layers[[1]]$data$imagerow <- -p$layers[[1]]$data$imagerow 
#p
pdf("/home/nasim/Visium_737/genemarkers.pdf")
SpatialFeaturePlot(Visium_737_v2, features = 'BANK1', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'CD3D', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'MS4A6A', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'CRIP2', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'COL11A1', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'SPARCL1', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'FABP4', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'IGFBP7', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'EPCAM', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'MKI67', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'CD68', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'MS4A1', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'JCHAIN', alpha = c(0.1, 1))
SpatialFeaturePlot(Visium_737_v2, features = 'PDGFRB', alpha = c(0.1, 1)) + theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size=20)) 
dev.off()

sInfo <- sessionInfo()

#make mat which is matrix array to data frame
mat.new <- as.data.frame(mat)

#B
B <- subset(mat.new, mat.new$`B-cells` > 0.6 & mat.new$Plasmablasts > 0.1)
B <- write.csv(B, file="/home/nasim/Visium_737/B.csv")

#PVL
PVL <- subset(mat.new, mat.new$PVL > 0.3 & mat.new$`Normal Epithelial`> 0.10 & mat.new$CAFs > 0.153)
PVL <- write.csv(PVL, file= "/home/nasim/Visium_737/PVL.csv")

#CAFs
CAFs <- subset(mat.new, mat.new$CAFs > 0.35 & mat.new$`B-cells`> 0.16 &mat.new$Myeloid==0 & mat.new$Plasmablasts ==0 & mat.new$`T-cells`> 0.13)
write.csv(CAFs, file="/home/nasim/Visium_737/CAFs.csv")

#cancerepithelial
cancerepithelial <- subset(mat.new, mat.new$`Cancer Epithelial` != 0 & mat.new$`B-cells` !=0 & mat.new$`Normal Epithelial`!= 0 & mat.new$CAFs==0 & mat.new$Myeloid==0 & mat.new$Plasmablasts==0 & mat.new$PVL==0 & mat.new$`T-cells`==0)
write.csv(cancerepithelial, file ="/home/nasim/Visium_737/cancerepithelial.csv")

#T-cells
T <- subset(mat.new, mat.new$`Cancer Epithelial` == 0 & mat.new$`B-cells` ==0 & mat.new$`Normal Epithelial`== 0 & mat.new$Endothelial==0 & mat.new$CAFs!=0 & mat.new$Myeloid==0 & mat.new$Plasmablasts==0 & mat.new$PVL==0 & mat.new$`T-cells` > 0.6)
write.csv(T, file ="/home/nasim/Visium_737/T.csv")

#the colors
#000000" = black"
#004949" = dark green
#009292"= blue green
"#ff6db6" = "pink"
"#ffb6db"= "ligt pink"
"#490092"= "dark purple"
"#006ddb"= "blue"
"#b66dff"= "light purple"
"#6db6ff"= "light blue"
"#b6dbff"= "pale blue"
"#920000"= "red"
"#924900"= "brown"
"#db6d00"= "orange"
"#24ff24"= "light green"
"#ffff6d"=  "yellow"

#piechart B, PVL, CAFs, cancer epithelial, T
pdf("/home/nasim/Visium_737/piecharts.pdf")
ggplot(B, aes(x = "", y=values, fill=celltypes)) +
  geom_col() +
  # to change the orientation of the pies use the direction in coor_polar
  coord_polar(theta = "y", direction = -1)+
  geom_text(aes(label = values),
            position = position_stack(vjust = 0.6), color="white") +
  scale_fill_manual(values=c("#000000", "#920000"))+ theme_void()

ggplot(PVL, aes(x = "", y=values, fill= celltypes)) +
  geom_col() +
  # to change the orientation of the pies use the direction in coor_polar
  coord_polar(theta = "y", direction = -1)+
  geom_text(aes(label = values),
            position = position_stack(vjust = 0.5), color="white") +
  scale_fill_manual(values=c("#004949","#6db6ff","#924900")) + theme_void()

ggplot(CAFs, aes(x = "", y=values, fill= celltypes)) +
  geom_col() +
  # to change the orientation of the pies use the direction in coor_polar
  coord_polar(theta = "y", direction = -1)+
  geom_text(aes(label = values),
            position = position_stack(vjust = 0.5), color="red") +
  scale_fill_manual(values=c("#000000","#004949","#ffff6d")) + theme_void()

ggplot(cancerepithelial, aes(x = "", y=values, fill= celltypes)) +
  geom_col() +
  # to change the orientation of the pies use the direction in coor_polar
  coord_polar(theta = "y", direction = -1)+
  geom_text(aes(label = values),
            position = position_stack(vjust = 0.5), color="white") +
  scale_fill_manual(values=c("#000000","#ff6db6","#6db6ff")) + theme_void()

ggplot(T, aes(x = "", y=values, fill= celltypes)) +
  geom_col() +
  # to change the orientation of the pies use the direction in coor_polar
  coord_polar(theta = "y", direction = -1)+
  geom_text(aes(label = values),
            position = position_stack(vjust = 0.5), color="red") +
  scale_fill_manual(values=c("#004949","#ffff6d")) + theme_void()

dev.off()

#save the new file
saveRDS(Visium_737_v2, "~/Visium_737/Visium_737_v3.rds")