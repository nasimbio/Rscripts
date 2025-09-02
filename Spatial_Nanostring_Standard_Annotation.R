# ------------------------------------------------------------------------------
#title: "Nanostring CosMx Spatial Transcriptomics — Tumor vs. Control Analysis"
#Author: Nasim Rahmatpour
#Date: "2025-03-31"
# ---

# Project Overview
#This project analyzes CosMx spatial transcriptomics data from two slides, each containing two distinct tissue regions, totaling **three tumor** and **one control** sample. The analysis integrates spatial gene expression with cell type annotation via reference single-cell RNA-seq data, aiming to explore spatial clustering and cell type composition across conditions.

#Spatial transcriptomic profiles were generated using **Nanostring CosMx SMI** technology. Two slides were analyzed, each containing multiple tissue sections. A reference scRNA-seq dataset with annotated cell types was used to transfer labels to spatial cells using **Seurat’s `MapQuery()`** method. The focus was to evaluate how cell types are distributed in tumor versus control and identify key expression patterns across spatial regions.


---

###Task 1: Preprocessing and Spatial Clustering 

#install.packages("data.table")
#install.packages("progressr")
#install.packages("future")
library(Seurat)
library(data.table)
library(progressr)
library(future)
library(sp)
library(dplyr)

#making seurat for KOKOM1
data.dir <- "/home/nasim/augusta.abio/Liver_KOKOM1/KOKOM1"
fov_df <- fread(file.path(data.dir,"KOKOM1_fov_positions_file.csv.gz"))
head(fov_df)
# Count unique FOVs
length(unique(fov_df$FOV)) 

meta <- fread(file.path(data.dir, "KOKOM1_metadata_file.csv.gz"))
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell_id

nano.obj_KOKOM1 <- LoadNanostring(data.dir = "/home/nasim/augusta.abio/Liver_KOKOM1/KOKOM1", fov = "KOKOM1")
nano.obj_KOKOM1 <- AddMetaData(object=nano.obj_KOKOM1, metadata = meta)


#add ploygons
polygons <- fread(file.path(data.dir, "KOKOM1-polygons.csv.gz"))
head(colnames(nano.obj_KOKOM1))
head(nano.obj_KOKOM1@meta.data$cell_ID)
head(nano.obj_KOKOM1@meta.data$cell_id)



# adding spatial location directly from metada info
coords <- nano.obj_KOKOM1@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
colnames(coords) <- c("SPATIAL_1", "SPATIAL_2")
coords<- coords[colnames(nano.obj_KOKOM1), , drop = FALSE]

# Attach spatial coordinates to Seurat object
nano.obj_KOKOM1[["spatial"]] <- CreateDimReducObject(
  embeddings = as.matrix(coords),
  key = "SPATIAL_",
  assay = DefaultAssay(nano.obj_KOKOM1)
)

DimPlot(nano.obj_KOKOM1, reduction = "spatial")



molecules <- fread(file.path(data.dir,"KOKOM1_tx_file.csv.gz"))
nano.obj_KOKOM1@misc$molecules <- molecules
nano.obj_KOKOM1@misc$fov_info <- fov_df

#add Group based on FOV
nano.obj_KOKOM1@meta.data$Group <- case_when(
  nano.obj_KOKOM1@meta.data$fov >= 1 & nano.obj_KOKOM1@meta.data$fov <= 123 ~ "Group3",
  nano.obj_KOKOM1@meta.data$fov >= 124 & nano.obj_KOKOM1@meta.data$fov <= 213 ~ "Group4",
  TRUE ~ "Unknown"
)

#add Group info
nano.obj_KOKOM1@meta.data$Group_info <- ""
nano.obj_KOKOM1@meta.data$Group_info[nano.obj_KOKOM1@meta.data$Group=="Group3"] <- "KO tumor in WT mice"
nano.obj_KOKOM1@meta.data$Group_info[nano.obj_KOKOM1@meta.data$Group=="Group4"] <- "KO tumor in KO mice"

nano.obj_KOKOM1[["percent.mt"]] <- PercentageFeatureSet(nano.obj_KOKOM1, pattern = "^mt-")

#making seurat for WTKM1

data.dir1 <- "/home/nasim/augusta.abio/Liver_WTKOM1/WTKOM1"
fov_df <- fread(file.path(data.dir1,"WTKOM1_fov_positions_file.csv.gz"))
head(fov_df)
length(unique(fov_df$FOV)) 

meta <- fread(file.path(data.dir1, "WTKOM1_metadata_file.csv.gz"))
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell_id

nano.obj_WTKOM1 <- LoadNanostring(data.dir = "/home/nasim/augusta.abio/Liver_WTKOM1/WTKOM1", fov = "WTKOM1")
nano.obj_WTKOM1 <- AddMetaData(object=nano.obj_WTKOM1, metadata = meta)



#add ploygons
polygons <- fread(file.path(data.dir1, "WTKOM1-polygons.csv.gz"))


head(colnames(nano.obj_WTKOM1))
head(nano.obj_WTKOM1@meta.data$cell_ID)
head(nano.obj_WTKOM1@meta.data$cell_id)




# this is for adding sptial location directly from metada info
coords <- nano.obj_WTKOM1@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
colnames(coords) <- c("SPATIAL_1", "SPATIAL_2")
coords<- coords[colnames(nano.obj_WTKOM1), , drop = FALSE]

# Attach spatial coordinates to Seurat object
nano.obj_WTKOM1[["spatial"]] <- CreateDimReducObject(
  embeddings = as.matrix(coords),
  key = "SPATIAL_",
  assay = DefaultAssay(nano.obj_WTKOM1)
)

DimPlot(nano.obj_WTKOM1, reduction = "spatial")



molecules <- fread(file.path(data.dir1,"WTKOM1_tx_file.csv.gz"))
nano.obj_WTKOM1@misc$molecules <- molecules
nano.obj_WTKOM1@misc$fov_info <- fov_df

#add Group based on FOV
nano.obj_WTKOM1@meta.data$Group <- case_when(
  nano.obj_WTKOM1@meta.data$fov >= 1 & nano.obj_WTKOM1@meta.data$fov <= 156 ~ "Group1",
  nano.obj_WTKOM1@meta.data$fov >= 157 & nano.obj_WTKOM1@meta.data$fov <= 331 ~ "Group2",
  TRUE ~ "Unknown"
)

#add Group info
nano.obj_WTKOM1@meta.data$Group_info <- ""
nano.obj_WTKOM1@meta.data$Group_info[nano.obj_WTKOM1@meta.data$Group=="Group1"] <- "WT tumor in WT mice"
nano.obj_WTKOM1@meta.data$Group_info[nano.obj_WTKOM1@meta.data$Group=="Group2"] <- "WT tumor in KO mice"

nano.obj_WTKOM1[["percent.mt"]] <- PercentageFeatureSet(nano.obj_WTKOM1, pattern = "^mt-")

out_dir <- "/home/nasim/augusta.abio/outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(future)
plan(sequential)
options(future.globals.maxSize = 100 * 1024^3)  # 100GB
nano.merged <- merge(nano.obj_WTKOM1, nano.obj_KOKOM1)

nano.merged <- JoinLayers(nano.merged)

#for merged data
# Build spatial coordinates matrix
coords <- nano.merged@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
colnames(coords) <- c("SPATIAL_1", "SPATIAL_2")
coords <- coords[colnames(nano.merged), , drop = FALSE]

# Re-create spatial reduction
nano.merged[["spatial"]] <- CreateDimReducObject(
  embeddings = as.matrix(coords),
  key = "SPATIAL_",
  assay = DefaultAssay(nano.merged)
)

#saveRDS(nano.merged, file.path(out_dir, "nano.merged.rds"))
saveRDS(nano.obj_WTKOM1, file.path(out_dir, "WTKOM1.rds"))
saveRDS(nano.obj_KOKOM1, file.path(out_dir, "KOKOM1.rds"))

#Before filter Stat
pdf(file.path(out_dir, "Beforefilter_Violin.pdf"), width = 20,height = 10)
tf <- function(x){tapply(x, Idents(nano.merged), median)}
Idents(nano.merged) <- nano.merged$Group
median.value <- apply(nano.merged@meta.data[,c("nCount_RNA", "nFeature_RNA", "percent.mt")],2,tf)
write.csv(median.value, file.path(out_dir,'Beforefilter.Median.csv'))
saveRDS(nano.merged, file.path(out_dir, '0-Data_Merged_RawObject.rds'))
VlnPlot(nano.merged, features = c('nCount_RNA','nFeature_RNA','nCount_HTO','percent.mt'),group.by='orig.ident', pt.size = 0,ncol = 3)
write.csv(table(nano.merged@meta.data$Group),file.path(out_dir,'Beforefilter.cellnumber_per_Group.csv'))
dev.off()

pdf(file.path(out_dir, "Beforefilter_Spatial.pdf"), width = 15,height = 10)
library(ggplot2)
DimPlot(
  nano.merged,
  reduction = "spatial",
  group.by = "Run_Tissue_name",  # or group by anything you like
  split.by = "Group",             # THIS is the key
  pt.size = 1.2
) +
  ggtitle("Spatial Clusters Split by Experimental Group")

dev.off()


#filter
nano.obj.filtered <- subset(
  nano.merged,
  subset = nCount_RNA > 150 & nFeature_RNA > 50
)

# Extract filtered cell names
filtered_cells <- colnames(nano.obj.filtered)

# Get the original coordinates
spatial_coords <- nano.merged@reductions$spatial@cell.embeddings

# Subset to filtered cells only
spatial_coords_filtered <- spatial_coords[filtered_cells, , drop = FALSE]

nano.obj.filtered[["spatial"]] <- CreateDimReducObject(
  embeddings = spatial_coords_filtered,
  key = "SPATIAL_",
  assay = DefaultAssay(nano.obj.filtered)
)




#After filter Stat
pdf(file.path(out_dir, "Afterfilter_Violin.pdf"), width = 20,height = 10)
tf1 <- function(x){tapply(x, Idents(nano.obj.filtered), median)}
Idents(nano.obj.filtered) <- nano.obj.filtered$Group
median.value <- apply(nano.obj.filtered@meta.data[,c("nCount_RNA", "nFeature_RNA", "percent.mt")],2,tf1)
write.csv(median.value, file.path(out_dir,'Afterfilter.Median.csv'))
saveRDS(nano.obj.filtered, file.path(out_dir, '1-Data_Merged_QC_Filtered.rds'))
VlnPlot(nano.obj.filtered, features = c('nCount_RNA','nFeature_RNA','nCount_HTO','percent.mt'),group.by='orig.ident', pt.size = 0,ncol = 3)
write.csv(table(nano.obj.filtered@meta.data$Group),file.path(out_dir,'AFterfilter.cellnumber_per_Group.csv'))
dev.off()

pdf(file.path(out_dir, "Afterfilter_Spatial.pdf"), width = 15,height = 10)
library(ggplot2)
DimPlot(
  nano.obj.filtered,
  reduction = "spatial",
  group.by = "Run_Tissue_name",  # or group by anything you like
  split.by = "Group",             # THIS is the key
  pt.size = 1.2
) +
  ggtitle("Spatial Clusters Split by Experimental Group")

dev.off()

#standard
`1-Data_Merged_QC_Filtered` <- NormalizeData(`1-Data_Merged_QC_Filtered`, normalization.method = "LogNormalize", scale.factor = 10000)
VariableFeatures(`1-Data_Merged_QC_Filtered`) <- rownames(`1-Data_Merged_QC_Filtered`)
`1-Data_Merged_QC_Filtered` <- ScaleData(`1-Data_Merged_QC_Filtered`)

`1-Data_Merged_QC_Filtered` <- RunPCA(`1-Data_Merged_QC_Filtered`, features = VariableFeatures(`1-Data_Merged_QC_Filtered`))

elbow_path <- file.path(out_dir, "Elbow_plot.pdf")
pdf(elbow_path, width = 15, height = 10)
p1 <- ElbowPlot(`1-Data_Merged_QC_Filtered`,ndims = 50)
print(p1)
dev.off()

heatmap_path <- file.path(out_dir, "PC_heatmap_plot.pdf")
pdf(heatmap_path, width = 15, height = 20)
p2 <-DimHeatmap(`1-Data_Merged_QC_Filtered`, dims = 1:30, cells = 500, balanced = TRUE)
print(p2)
dev.off()

`1-Data_Merged_QC_Filtered` <- FindNeighbors(`1-Data_Merged_QC_Filtered`, dims = 1:30)
`1-Data_Merged_QC_Filtered` <- FindClusters(`1-Data_Merged_QC_Filtered`, resolution = 0.5)
`1-Data_Merged_QC_Filtered` <- RunUMAP(`1-Data_Merged_QC_Filtered`, dims = 1:30)
#`1-Data_Merged_QC_Filtered`<- RunTSNE(`1-Data_Merged_QC_Filtered`, dims = 1:30)

plots_path <- file.path(out_dir, "umap_spatial_plots.pdf")
pdf(plots_path, width = 12, height = 10)
DimPlot(`1-Data_Merged_QC_Filtered`,group.by = 'seurat_clusters',reduction = 'umap',label=T,repel=T)
DimPlot(`1-Data_Merged_QC_Filtered`, group.by = 'seurat_clusters', split.by= 'Group',reduction = 'umap',label=T,repel=T)


DimPlot(`1-Data_Merged_QC_Filtered`, group.by = 'seurat_clusters', split.by= 'Run_Tissue_name',reduction = 'umap',label=T,repel=T)

DimPlot(`1-Data_Merged_QC_Filtered`, reduction = "spatial", group.by = "seurat_clusters", split.by = "Group", pt.size = 2) +
  ggtitle("Spatial Distribution of Clusters")

DimPlot(`1-Data_Merged_QC_Filtered`, reduction = "spatial", group.by = "seurat_clusters",  split.by = "Run_Tissue_name", pt.size = 1.2
) +
  ggtitle("Spatial Clusters Split by Experimental Group")


dev.off()

#Find gene markers for each cluster
markers <- FindAllMarkers(`1-Data_Merged_QC_Filtered` , only.pos = TRUE)

write.csv(markers, file.path(out_dir,"RNA_markers.csv"), row.names= TRUE)

#top20
top20.markers_RNA <- markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top20 <- c()
for(c in unique(top20.markers_RNA$cluster)){
  tmp <- top20.markers_RNA[top20.markers_RNA$cluster==c,]$gene
  top20 <- cbind(top20,tmp)
}
colnames(top20) <- paste0('C_',unique(top20.markers_RNA$cluster))
write.csv(top20,file.path(out_dir,"top20_RNA_markers.csv"))


#save
saveRDS(`1-Data_Merged_QC_Filtered`, file.path(out_dir, "2-Data_Merged_DimR.rds"))

#canonical markers feature plot
pdf(file.path(out_dir,"canonical_markers_featureplot.pdf"),width = 10,height = 10)

#CD8 Tcells
# P1=FeaturePlot(`2-Data_Merged_DimR` , features = c("Gpr183", "Cx3cr1", "Rgs1", "Gzmk", "Cxcr5", "Cd6", "Klrb1", "Hspa1a", "Mki67"),  ncol = 3)

#CD4 Tcell
# P2=FeaturePlot(`2-Data_Merged_DimR` , features = c("Ccr7", "Anxa1", "Rgs1", "Tcf7", "Cxcr5", "Gzmk", "Cxcr6", "Cxcl13", "Foxp3", "Il10", "Ctla4", "Hspa1a", "Mki67", "Fcgr3a", "Cd27"),  ncol = 3)

#B cells
# P3=FeaturePlot(`2-Data_Merged_DimR`, features = c("Cd79b", "Cd56", "Gpr183", "Hspa1a", "Mki67", "Igha1", "Ighg1"),  ncol = 3)


#Myeloid cells
# P4=FeaturePlot(`2-Data_Merged_DimR` , features = c("Lamp3", "Mki67", "Clec9a", "C1qc", "Clec10a", "Timp1", "Lilra4", "Cd14", "Fcgr3a", "Nlrp3", "Il1b", "Cxcl12", "C1ac", "Spp1", "Cpa3"),  ncol = 3)

#ILCs
# P5=FeaturePlot(`2-Data_Merged_DimR` , features = c("Kit", "Fcgr3a", "Sell", "Dusp4", "Icam1", "Mki67", "Il7r"),  ncol = 3)


#CD4, NaiveT
P1=FeaturePlot(`2-Data_Merged_DimR` , features = c("Tcf7", "Lef1", "Ccr7", "Sell", "Mal", "Il7r", "Trac", "Ltb","Cd11a", "Cd14", "Cd19", "Cd25", "Cd27", "Cd28", "Cd44", "Cd45ra", "Cd45ro", "Cd57", "Cd621", "Cd122", "Cd127", "Tfn-y", "Il-2", "Tnf"),  ncol = 3)
#CD4
P2=FeaturePlot(`2-Data_Merged_DimR` , features = c("Cd4","Cd3","Cd3e","Cd3d", "Cd2", "Ccr1", "Ccr5", "Cxcr3", "Stat4", "T-bet", "Ccr3", "Ccr4", "Ccr8"," Cxcr4", "Stat5", "Stat6", "Gata-3", "Ccr6", "roryT", "Ccr10", "Fgf", "Ahr9", "Btla", "Bcl-6", "Stat","Il2ra", "Foxp3","Blk"),  ncol = 3)


#CD8
P3=FeaturePlot(`2-Data_Merged_DimR`, features = c("Prf1", "Klrg1", "Gzmb", "Cd8b1","Cd8a", "Ca8b", "Cx3cr1", "Gzmh", "Tbx21", "Eomes","Zfp36l2", "Gzmk", "Gzma", "Ccl5", "Infg", "Runx3",  "Cd69", "Itgae"),  ncol = 3)

#other T markers
P4=FeaturePlot(`2-Data_Merged_DimR`, features = c("Ctla4", "Icos", "Tigit", "Pdcd1"),  ncol = 3)

#NK
P5=FeaturePlot(`2-Data_Merged_DimR`, features = c("Klrd1", "Tyrobp", "Nkg7", "Ncam1", "Klra3", "Klra7", "Klrg1", "Kir2dl1", "Kir2dl3", "Kir2dl4", "Kir3dl1", "Kir3dl2", "Kir3dl3", "Kir2ds4"),  ncol = 3)

P6=FeaturePlot(`2-Data_Merged_DimR`, features = c("Ifit1", "Ifit2", "Ifit3"),  ncol = 1)


#Monocytes
P7=FeaturePlot(`2-Data_Merged_DimR`, features = c("Ctsa", "Fcn1", "Neat1", "Lyz", "Psap", "Aif1", "Mnda", "Serpina1", "Tyrobp"),  ncol = 3)

#Monocytes
P8=FeaturePlot(`2-Data_Merged_DimR` , features = c("Cd14", "Itgam","S100a8", "S100a9", "Fcn1", "Vcan", "Il1b"), ncol=3)

#Monocytes
P9=FeaturePlot(`2-Data_Merged_DimR`, features = c("Fcgr3a", "Cd68", "Vmo1", "Rhoc", "Ms4a7", "Ifitm3", "Aif1", "Lst1"),  ncol = 3)
#Neutrophils
P10=FeaturePlot(`2-Data_Merged_DimR`, features = c("Ceacam8", "Arg1", "Flt1", "Icam1", "Il17ra", "Cxcr4"),  ncol = 3)

#B cells
P11=FeaturePlot(`2-Data_Merged_DimR`, features = c("Cd79a", "Ms4a1", "Bank1", "Cd19", "Cd79b", "Igkc"),  ncol = 3)

#Plasma:
P12=FeaturePlot(`2-Data_Merged_DimR`, features = c("Jchain", "Mzb1"))

#Dendritic cells
P13=FeaturePlot(`2-Data_Merged_DimR`, features = c("Lilra4", "Cd74", "Hla-dpa1", "Hla-dpb1", "Hla-dqa1", "Ccdc88a", "Hla-dra"))

#Other cells (Plateles/RBCs)
#P14=FeaturePlot(`2-Data_Merged_DimR`, features = c("Gng11", "Cavin2", "Tubb1"))


#Proliferation
P15=FeaturePlot(`2-Data_Merged_DimR`, features = c("MkI67", "Stmn1"))

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
print(P11)
print(P12)
print(P13)
#print(P14)
print(P15)

dev.off()

###Task 2: Cell Type Annotation Using Single Cell Reference

library(Matrix)
library(Seurat)
library(stringr)

# Read matrix, need to transverse the matrix, since the matrix coming from h5ad has cell_names as rows and gene_names as columns
expr <- readMM("/home/nasim/augusta.abio/expr_matrix.mtx")
expr <- t(expr)

#gene_names and cell_names is a dataframe with one column, here get the contents of the only column, it will be saved as vector of characters
gene_names <- read.csv("/home/nasim/augusta.abio/gene_names.csv", header = TRUE)[[1]]
length(gene_names)  # Should be 17984
rownames(expr) <- gene_names


cell_names <- read.csv("/home/nasim/augusta.abio/cell_names.csv", header = TRUE)[[1]]
length(cell_names)  # Should be 7294
colnames(expr) <- cell_names

#these didn't work, since there is a header
#rownames(expr) <- read.csv("/home/nasim/augusta.abio/gene_names.csv", header = FALSE)[[1]]
#colnames(expr) <- read.csv("/home/nasim/augusta.abio/cell_names.csv", header = FALSE)[[1]]


seurat_obj <- CreateSeuratObject(counts = expr)

# Add metadata
meta <- read.csv("/home/nasim/augusta.abio/obs_metadata.csv", row.names = 1)
seurat_obj <- AddMetaData(seurat_obj, metadata = meta)
print(dim(seurat_obj))

#remove "Unknown Cells" 
#seurat_obj <- subset(seurat_obj, subset=ann_level_2!="Unknown")
#print(dim(seurat_obj))

#the var meta has gene names in different format, for now add the whole file to the seurat object
var_meta <- read.csv("/home/nasim/augusta.abio/var_metadata.csv", row.names = 1)
var_meta$gene_symbol <- ""
var_meta$gene_symbol <- apply(var_meta ,1,function(x) str_split(x[7], "\\_"))
var_meta$gene_symbol <- apply(var_meta,1,function(x) x[[12]][[1]][1])

# Make sure rownames match Seurat object features
var_meta <- var_meta[rownames(seurat_obj), , drop = FALSE]
seurat_obj[["RNA"]]@misc$features <- var_meta


# Make sure rownames of var_meta match current Seurat object gene IDs
gene_symbols <- var_meta[rownames(seurat_obj), "gene_symbol"]

# Optional: check for duplicated gene symbols and make unique if needed
gene_symbols_unique <- make.unique(gene_symbols)

# Set new rownames
rownames(seurat_obj) <- gene_symbols_unique




#Task 2: Cell Type Annotation Using Single Cell Reference

out_dir <- "/home/nasim/augusta.abio/MapQuery_outs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

common.genes <- intersect(rownames(seurat_obj), rownames(`2-Data_Merged_DimR`))
length(common.genes)


#standard on reference and query. query standard is already done
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, dims = 1:30)
#seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, umap.method = "uwot",return.model = TRUE  # <-- this is the critical part
)
pdf(file.path(out_dir, 'reference_umap.pdf'), width = 14,height = 10)
DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("Reference Annotaion")
dev.off() 

library(openxlsx)
counts <- table(seurat_obj@meta.data$cell_type)
wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, sheet = "Sheet1", x = counts)
saveWorkbook(wb, file.path(out_dir, "reference_cell_types_counts.xlsx"))


options(future.globals.maxSize = 100 * 1024^3)  # 100 GB
anchors <- FindTransferAnchors( reference = seurat_obj,
                                query = `2-Data_Merged_DimR`,
                                features = common.genes, 
                                dims = 1:30,
                                reference.reduction = "pca"
)

options(future.globals.maxSize = 200 * 1024^3)  # 200 GB
mapped <- MapQuery(
  anchorset = anchors,
  query = `2-Data_Merged_DimR`,
  reference = seurat_obj,
  refdata = list(cell_type = "cell_type"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

plots_path <- file.path(out_dir, "predicted_celltypes_umap_spatial_plots.pdf")
pdf(plots_path, width = 14, height = 10)

DimPlot(mapped, reduction = "ref.umap", group.by = "predicted.cell_type", label = TRUE, repel = TRUE) +
  ggtitle("CosMx Cells Projected onto Reference UMAP Sapce")


DimPlot(mapped, reduction = "umap", group.by = "predicted.cell_type", label = TRUE, repel = TRUE) +
  ggtitle("Predicted Cell Types in CosMx UMAP Space")


DimPlot(mapped, group.by = 'predicted.cell_type', split.by= 'Group',reduction = 'umap',label=T,repel=T)+
  ggtitle("Predicted Cell Types in CosMx Split by Groups")


DimPlot(mapped, group.by = 'predicted.cell_type', split.by= 'Run_Tissue_name',reduction = 'umap',label=T,repel=T) +
  ggtitle("Predicted Cell Types in CosMx Split by Tissue")

DimPlot(mapped, reduction = "spatial", group.by = "predicted.cell_type", split.by = "Group", pt.size = 2) +
  ggtitle("Spatial Distribution of Predicted Cell Types")

DimPlot(mapped, reduction = "spatial", group.by = "predicted.cell_type",  split.by = "Run_Tissue_name", pt.size = 1.2
) +
  ggtitle("Spatial Distribution of Predicted Cell Types Split by Experimental Group")


dev.off()

pdf(file.path(out_dir, 'predicted.cell_type.score_umap.pdf'), width = 12,height = 10)
FeaturePlot(mapped, features = "predicted.cell_type.score", reduction = "ref.umap") +
  ggtitle("Prediction Confidence (Score)")

FeaturePlot(mapped, features = "predicted.cell_type.score", reduction = "umap") +
  ggtitle("Prediction Confidence (Score)")

#VlnPlot(mapped, features = "predicted.cell_type.score", group.by = "seurat_clusters")
dev.off()


mapped$predicted.cell_type.filtered <- ifelse(
  mapped$prediction.score.cell_type < 0.5,
  "Unknown",
  mapped$predicted.cell_type
)

library(openxlsx)
clusters <- table(mapped@meta.data$predicted.cell_type)
group <- table(mapped@meta.data$predicted.cell_type, mapped@meta.data$Group)
tissue <- table(mapped@meta.data$predicted.cell_type, mapped@meta.data$Run_Tissue_name)
wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
addWorksheet(wb, "Sheet2")
addWorksheet(wb, "Sheet3")
writeData(wb, sheet = "Sheet1", x = clusters )
writeData(wb, sheet = "Sheet2", x = group)
writeData(wb, sheet = "Sheet3", x = tissue)
saveWorkbook(wb, file.path(out_dir, "predicted_cell_types_counts.xlsx"))

saveRDS(mapped, file.path(out_dir, "MappedQuery.RDS"))
