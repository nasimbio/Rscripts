
# ------------------------------------------------------------------------------
# title: "10x Multiomic RNA + ATAC Profiling — Pig Kidney (Single-Nuclei)"
# Author: Nasim Rahmatpour
# Date: 2023-12-13
# ------------------------------------------------------------------------------

# Project Overview
#Multiomic single-nuclei data from two pig kidney samples (P0906 and P0909) were processed to explore gene expression and chromatin accessibility using the WNN (Weighted Nearest Neighbor) approach. This project aims to:
  
# - Perform modality-specific QC and filtering on RNA and ATAC data
# - Apply WNN-based multimodal integration to define biologically meaningful clusters
# - Annotate cell types based on known kidney-specific markers from literature

---


###Task 1: Data Processing and Quality Control
  
#libraries
library(Seurat) 
suppressMessages(require(cowplot))
suppressMessages(require(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(ggplot2))
#suppressMessages(library(Nebulosa))
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(patchwork))
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(biovizBase)

#create seurat object for each sample and add RNA and ATAC as different assay
dir.create("./Dinaqor2/dinaqor_RNA_ATAC.out")
ids <- c("P0906","P0909")
for (file in ids){
  sample <- Read10X_h5(filename=paste0("/home/nasim/Dinaqor/Dinaqor2/", file, "/filtered_feature_bc_matrix.h5"))
  rna_counts <- sample$`Gene Expression`
  atac_counts <- sample$Peaks
  metadata_cellranger<-read.csv(paste0("/home/nasim/Dinaqor/Dinaqor2/", file, "/per_barcode_metrics.csv")) 
  row.names(metadata_cellranger) <- metadata_cellranger$barcode
  sampleseurat <- CreateSeuratObject(counts = rna_counts, assay = "RNA", project = file)
  frag.file <- paste0("/home/nasim/Dinaqor/Dinaqor2/", file, "/atac_fragments.tsv.gz")
  sampleseurat[["ATAC"]] <- CreateChromatinAssay(counts=atac_counts, sep = c(":", "-"), fragments=frag.file, min.cells = 10)
  DefaultAssay(sampleseurat) <- "ATAC"
  sampleseurat <- NucleosomeSignal(sampleseurat)
  sampleseurat <- TSSEnrichment(sampleseurat)
  sampleseurat<-AddMetaData(sampleseurat,metadata=metadata_cellranger)
  assign(file, sampleseurat)
}

#add columns ro metadata
add_meta <- function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(object=x, pattern = "^MT-")
  x[["percent.hb"]] <- PercentageFeatureSet(object=x, pattern = "^HB")
  x[["percent.ribo"]] <- PercentageFeatureSet(object=x, pattern = "^RP[SL]")
  return(x)
}

`P0906`<- add_meta(`P0906`)
`P0909` <- add_meta(`P0909`)

#violin plot before filtering
pdf(paste0("/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/","QCViolin.pdf"), width=20, height=10)
violinplot<- function(x){
  VlnPlot(x, features =c('nFeature_RNA','nCount_RNA','nCount_ATAC','nFeature_ATAC','percent.mt', 'percent.ribo', 'percent.hb'), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()
} 
P1 <- violinplot(`P0906`)
P2 <- violinplot(`P0909`)
print(P1)
print(P2)

#get the median before filtering
get_median <- function(x){
  metax <- x@meta.data
  print(dim(metax))
  median.value <- apply(metax[,c('nFeature_RNA','nCount_RNA','nCount_ATAC','nFeature_ATAC','percent.mt', 'percent.ribo', 'percent.hb')], 2, median)
  return(median.value)
}

#save
write.csv(get_median(`P0906`), paste0("/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/", "P0906", "_QCmedian.csv"))
write.csv(get_median(`P0909`), paste0("/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/", "P0909", "_QCmedian.csv"))

#filter
filterseurat <- function(x){
  x <- subset(x = x, subset = nCount_ATAC < 7e4 & nCount_ATAC > 5e3 & nCount_RNA < 25000 & nCount_RNA > 1000 & percent.mt < 20
  )
}

`P0906`<- filterseurat(`P0906`)
`P0909` <- filterseurat(`P0909`)
print(dim(P0906))
print(dim(P0909))

#violin plot after filtering
pdf(paste0("./Dinaqor2/dinaqor_RNA_ATAC.out/","AfterQCViolin.pdf"), width=20, height=10)
violinplot<- function(x){
  VlnPlot(x, features =c('nFeature_RNA','nCount_RNA','nCount_ATAC','nFeature_ATAC','percent.mt', 'percent.ribo', 'percent.hb'), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()
}
P1 <- violinplot(`P0906`)
P2 <- violinplot(`P0909`)
print(P1)
print(P2)

#get the median after filtering
get_median <- function(x){
  metax <- x@meta.data
  print(dim(metax))
  median.value <- apply(metax[,c('nFeature_RNA','nCount_RNA','nCount_ATAC','nFeature_ATAC','percent.mt', 'percent.ribo', 'percent.hb')], 2, median)
  return(median.value)
}

#save
write.csv(get_median(`P0906`), paste0("./Dinaqor2/dinaqor_RNA_ATAC.out/", "P0906", "_AfterQCmedian.csv"))
write.csv(get_median(`P0909`), paste0("./Dinaqor2/dinaqor_RNA_ATAC.out/", "P0909", "_AfterQCmedian.csv"))



###Task 2: WNN Integration and Cell Type Annotation

#standard analysis for each RNA and ATAC data
standard_analysis <- function(x){
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  # ATAC analysis
  # exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(x) <- "ATAC"
  x <- RunTFIDF(x)
  x <- FindTopFeatures(x, min.cutoff = 'q0')
  x <- RunSVD(x)
  x <- RunUMAP(x, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  
  x <- FindMultiModalNeighbors(x, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  x <- RunUMAP(x, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  x <- FindClusters(x, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  
  return(x)
}

`P0906`<- standard_analysis(`P0906`)
`P0909` <- standard_analysis(`P0909`)

#visualize umap
pdf(paste0("/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/","P0906_umap.pdf"), width=20, height=10)
p1 <- DimPlot(`P0906`, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 7, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(`P0906`, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 7, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(`P0906`, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 7, repel = TRUE) + ggtitle("WNN")
print(p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()

pdf(paste0("/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/","P0909_umap.pdf"), width=20, height=10)
p1 <- DimPlot(`P0909`, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 7, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(`P0909`, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 7, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(`P0909`, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 7, repel = TRUE) + ggtitle("WNN")
print(p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
dev.off()

#finding the pg kidney markers from the pig paper
pdf(paste0("/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/", "P0906_canonical_pig_markers.pdf"))
DefaultAssay(`P0906`) <- 'SCT'
library(ggplot2)
library(RColorBrewer)
library(angrycell)
angrycell::DotPlot2(`P0906`
                    ,features = c( "AQP3", "GATA2", "AQP2", "PAX2", "TMEM52B", "CUBN", "LRP2", "SLC13A3", "SLC34A1",  "NPHS1", "NPHS2", "WT1", "CLIC5",  "SLC12A1", "COL6A1" , "DCN1", "CALD1", "COL1A1", "COL1A2", "KRT8", "KRT18", "PECAM1", "EPCAM", " NRP1", "CD3G", "CD3E", "CD3D", "CD2","CD4", "CD14", "CCR7", "TCF7", "IL2RA","FOXP3","IL7R", "CD8A", "CD8B", "NCAM1","KIT", "CSF3R","CSF1R","LILRA4","CD163", "CD68","FLT3", "JCHAIN","MZB1","CD79A","MS4A1","BLK",  "NOTCH1"), dot.scale = 8) +
  coord_flip() +
  
  theme(text = element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y = element_text(size =10),
        axis.title.x = element_text(size=10))
dev.off() 

pdf(paste0("/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/", "P0909_canonical_pig_markers.pdf"))
DefaultAssay(`P0909`) <- 'SCT'
library(ggplot2)
library(RColorBrewer)
library(angrycell)
angrycell::DotPlot2(`P0909`
                    ,features = c( "AQP3", "GATA2", "AQP2", "PAX2", "TMEM52B", "CUBN", "LRP2", "SLC13A3", "SLC34A1",  "NPHS1", "NPHS2", "WT1", "CLIC5",  "SLC12A1", "COL6A1" , "DCN1", "CALD1", "COL1A1", "COL1A2", "KRT8", "KRT18", "PECAM1", "EPCAM", " NRP1", "CD3G", "CD3E", "CD3D", "CD2","CD4", "CD14", "CCR7", "TCF7", "IL2RA","FOXP3","IL7R", "CD8A", "CD8B", "NCAM1","KIT", "CSF3R","CSF1R","LILRA4","CD163", "CD68","FLT3", "JCHAIN","MZB1","CD79A","MS4A1","BLK",  "NOTCH1"), dot.scale = 8) +
  coord_flip() +
  
  theme(text = element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y = element_text(size =10),
        axis.title.x = element_text(size=10))
dev.off() 

#add annotation
Idents(`P0906`) <- `P0906`@meta.data$seurat_clusters
`P0906`<- RenameIdents(`P0906`, `0` = "Proximal tubule cells", `1` = "Proximal tubule cells", `2` = "Proximal tubule cells", `3` = "Endothelial", `4` = "Loope of Henle cells", `5` = "Proximal tubule cells", `6` = "Distal convoluted tubule cells", `7` = "Proximal tubule cells",`8` = "Distal convoluted tubule cells", `9` = "Monocytes", `10`= "Immune cells T cells", `11` = "Fibroblast", `12`= "Epithelial", `13`= "Podocytes and Loope of Henle", `14`= "Collecting duct cells", `15` = "Endothelial", `16`= "Monocytes", `17` ="Proximal tubule cells", `18` = "Proximal tubule cells", `19` = "Endothelial", `20` = "Collecting duct cells", `21`=  "Fibroblast", `22`="Epithelial+collecting duct cells", `23`= "?", `24`="Loope of Henle cells", `25`="Loope of Henle cells", `26`= "Collecting duct cells", `27`="?" )

pdf(paste0('/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/P0906_annotation.pdf'),width=13, height=10)
DimPlot(`P0906`, reduction = "wnn.umap", label = TRUE, label.size = 5, repel = TRUE)
DimPlot(`P0906`, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")
dev.off()
`P0906`@meta.data$Annotation<- Idents(`P0906`)





Idents(`P0909`) <- `P0909`@meta.data$seurat_clusters
`P0909`<- RenameIdents(`P0909`, `0` = "Endothelial", `1` = "Loop of Henle cells", `2` = "Proximal tubule cells", `3` = "Immune cells T cells", `4` = "Epithelial", `5` = "Endothelial", `6` = "Collecting duct cells", `7` = "Collecting duct cells",`8` = "Epithelial", `9` = "Monocytes", `10`= "Collecting duct cells", `11` = "Collecting duct cells", `12`= "Fibroblast", `13`= "Fibroblast", `14`= "Collecting duct cells", `15` = "Proximal tubule cells", `16`= "Fibroblast", `17` ="Fibroblast", `18` = "Fibroblast", `19` = "Loop of Henle cells", `20` = "Loope of Henle cells+Collecting duct cells", `21`=  "Monocytes", `22`="Epithelial+ Collecting duct cells", `23`="Endothelial", `24`= "Distal covoluted tubule cells", `25`="Monocytes", `26`="Endothelial", `27`="Podocytes and Loope of Henle" )

pdf(paste0('/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/P0909_annotation.pdf'),width=13, height=10)
DimPlot(`P0909`, reduction = "wnn.umap", label = TRUE, label.size = 5, repel = TRUE)
DimPlot(`P0909`, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")
dev.off()
`P0909`@meta.data$Annotation<- Idents(`P0909`)

#save the object
saveRDS(`P0906`, paste0("/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/", "P0906_annotation.rds"))
saveRDS(`P0909`,paste0("/home/nasim/Dinaqor/Dinaqor2/dinaqor_RNA_ATAC.out/", "P0909_annotation.rds"))

