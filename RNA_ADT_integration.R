# ------------------------------------------------------------------------------
# title: "PBMC Single-Cell CITE-seq, RNA + ADT"
# Author: Nasim Rahmatpour
# Date: "6/13/2023"
# ------------------------------------------------------------------------------


# Project Overview

#This project analyzes peripheral blood mononuclear cells (PBMCs) from a clinical study using single-cell CITE-seq (RNA + surface protein) data. The goal is to characterize common and rare immune cell types by integrating transcriptomic and surface marker information, enabling a comprehensive view of immune cell diversity within human blood samples.


---

###Task 1: Data Loading, QC, and Visualization
library(Seurat) 
suppressMessages(require(cowplot))
suppressMessages(require(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(ggplot2))
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(patchwork))


dir.create("./batch7_outputs_filtered")

ids <- c("NGS008_7A_GEX_ADT","NGS008_7B_GEX_ADT",
         "NGS008_7C_GEX_ADT","NGS008_7D_GEX_ADT",
         "NGS008_7E_GEX_ADT","NGS008_7F_GEX_ADT")

#create a seurat object for each
for (file in ids){
  pbmcall <- Read10X(data.dir = paste0("./", file))
  pbmcalldata <- CreateSeuratObject(counts = pbmcall$`Gene Expression`, project = file) 
  pbmcalldata[['Protein']] = CreateAssayObject(counts = pbmcall$`Antibody Capture`)               
  assign(file, pbmcalldata)
}

#combine/merge the seurat objects
pbmc.combined <- merge(`NGS008_7A_GEX_ADT`, 
                       c(`NGS008_7B_GEX_ADT`,`NGS008_7C_GEX_ADT`,`NGS008_7D_GEX_ADT`,`NGS008_7E_GEX_ADT`,`NGS008_7F_GEX_ADT`), add.cell.ids = c("PBMC7A", "PBMC7B","PBMC7C","PBMC7D", "PBMC7E","PBMC7F"), project = "PBMCCite")

#add new column to the metadata of combined seurat
rb.genes <- "^RP[SL]"
pbmc.combined[["percent.ribo"]] <- PercentageFeatureSet(object = pbmc.combined, pattern = rb.genes)
hb.genes <- "^HB"
pbmc.combined[["percent.hb"]] <- PercentageFeatureSet(object = pbmc.combined, pattern = hb.genes)
mitochondrial_pattern <- "^MT-"
pbmc.combined[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.combined, pattern = mitochondrial_pattern)

#do filtering
pbmc.filter <- subset(x = pbmc.combined, subset = percent.mt < 20 & percent.hb < 50)
dim(pbmc.combined)
dim(pbmc.filter)

#violin plots after filtering
p1 <- VlnPlot(pbmc.filter,features=c("nFeature_RNA"))+ theme(legend.position = 'none')
p2 <- VlnPlot(pbmc.filter,features=c("nCount_RNA"))+ theme(legend.position = 'none')
p4 <- VlnPlot(pbmc.filter,features=c("nCount_Protein"))+ theme(legend.position = 'none')
p3 <- VlnPlot(pbmc.filter,features=c("nFeature_Protein"))+ theme(legend.position = 'none')
p5 <- VlnPlot(pbmc.filter,features=c("percent.hb"))+ theme(legend.position = 'none')
p6 <- VlnPlot(pbmc.filter,features=c("percent.mt"))+ theme(legend.position = 'none')
p7 <- VlnPlot(pbmc.filter,features=c("percent.ribo"))+ theme(legend.position = 'none')

#save the files
pdf(paste0("./batch7_outputs_filtered/","separate_feature_violin.pdf"),width = 20,height = 10)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
dev.off()

#defining a function to get median
tf1 <- function(x){tapply(x, Idents(pbmc.filter), median)}
pbmc.filter$Sample <- pbmc.filter$orig.ident

# median after QC
pdf(paste0("./batch7_outputs_filtered/","all_features_violin.pdf"),width = 20,height = 10)
Idents(pbmc.filter) <- pbmc.filter$Sample
median.value <- apply(pbmc.filter@meta.data[,c('nFeature_RNA','nCount_RNA','nCount_Protein','nFeature_Protein', 'percent.hb','percent.mt', 'percent.ribo')],2,tf1)
write.csv(median.value, paste0("./batch7_outputs_filtered/",'AfterQC.median.csv'))
saveRDS(pbmc.filter, paste0("./batch7_outputs_filtered/",'pbmc.filter.rds'))
VlnPlot(pbmc.filter, features = c('nFeature_RNA','nCount_RNA','nCount_Protein','nFeature_Protein', 'percent.hb','percent.mt', 'percent.ribo'),group.by='Sample', pt.size = 0,ncol = 3)
dev.off()

####Task 2: Independent RNA and Protein Integration
DefaultAssay(pbmc.filter) <- 'RNA'
pbmc.list <- SplitObject(pbmc.filter, split.by = "orig.ident")
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE, selection.method = "vst", nfeatures = 3000)
})

#Find the features 
features <- SelectIntegrationFeatures(object.list = pbmc.list)

#Find the anchors
RNA.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features)

#integrate the RNA data and save it as a new variable or rna.combined
rna.combined <- IntegrateData(anchorset = RNA.anchors)

#now rna.combined has 3 assays: RNA, Protein and Integrated which is the common data of all samples from RNA
#run the standard workflow 
DefaultAssay(rna.combined) <- "integrated"
rna.combined <- ScaleData(rna.combined, verbose = FALSE)
rna.combined <- RunPCA(rna.combined, npcs = 30, verbose = FALSE)

#now go back to pbmc.filter and change the default assay of filtered data to protein or ADT data and do the same and performing the protein integration
DefaultAssay(pbmc.filter) <- 'Protein'
pbmc.list <- SplitObject(pbmc.filter, split.by = "orig.ident")
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
  VariableFeatures(x) <- rownames(x[["Protein"]])
  x <- NormalizeData(x, verbose = FALSE, normalization.method = 'CLR', margin = 2)
})

#select the common features o
features <- SelectIntegrationFeatures(object.list = pbmc.list)

#find the anchors
ADT.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features)

#do the integration
adt.combined <- IntegrateData(anchorset = ADT.anchors)

#now the adt.combined has 3 assays too:RNA, Protein and Integrated but integrated assay here has the common protein info of the samples NOT the RNA
# Run the standard workflow 
DefaultAssay(adt.combined) <- "integrated"
adt.combined <- ScaleData(adt.combined, verbose = FALSE)
adt.combined <- RunPCA(adt.combined,reduction.name = 'apca')

DefaultAssay(rna.combined) <- "integrated"
DefaultAssay(adt.combined) <- "integrated"
#making matrix based on integrated adt
adt.data <- GetAssayData(object =  adt.combined[['integrated']], slot = 'data')
dim(adt.data)
adt.scaled <- GetAssayData(object =  adt.combined[['integrated']], slot = 'scale.data')
dim(adt.scaled)

#add the integrated adt data as an another assay in rna.combined
rna.combined[["integrated.adt"]] <- CreateAssayObject(data = adt.data )
rna.combined[['adt.pca']] <- adt.combined[['apca']]

#set the used assay in rna.combined for integrated.adt part as adt.pca 
rna.combined[['adt.pca']]@assay.used <-  "integrated.adt"
rna.combined[['adt.pca']]@assay.used

#set the default assay of rna.combined to integrated.adt to centre and scale the matrix
DefaultAssay(rna.combined) <- "integrated.adt"
rna.combined <- ScaleData(rna.combined)

#check the data and reduction method match
rna.combined[['pca']]@assay.used
rna.combined[['adt.pca']]@assay.used

#save the rna.combined 
saveRDS(rna.combined, paste0("./batch7_outputs_filtered/",'rna.combined.rds'))

#save the adt.combined 
saveRDS(adt.combined, paste0("./batch7_outputs_filtered/",'adt.combined.rds'))

###Task 2: Independent RNA and Protein Integration
DefaultAssay(pbmc.filter) <- 'RNA'
pbmc.list <- SplitObject(pbmc.filter, split.by = "orig.ident")
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE, selection.method = "vst", nfeatures = 3000)
})

#Find the features 
features <- SelectIntegrationFeatures(object.list = pbmc.list)

#Find the anchors
RNA.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features)

#integrate the RNA data and save it as a new variable or rna.combined
rna.combined <- IntegrateData(anchorset = RNA.anchors)

#now rna.combined has 3 assays: RNA, Protein and Integrated which is the common data of all samples from RNA
#run the standard workflow 
DefaultAssay(rna.combined) <- "integrated"
rna.combined <- ScaleData(rna.combined, verbose = FALSE)
rna.combined <- RunPCA(rna.combined, npcs = 30, verbose = FALSE)

#now go back to pbmc.filter and change the default assay of filtered data to protein or ADT data and do the same and performing the protein integration
DefaultAssay(pbmc.filter) <- 'Protein'
pbmc.list <- SplitObject(pbmc.filter, split.by = "orig.ident")
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
  VariableFeatures(x) <- rownames(x[["Protein"]])
  x <- NormalizeData(x, verbose = FALSE, normalization.method = 'CLR', margin = 2)
})

#select the common features o
features <- SelectIntegrationFeatures(object.list = pbmc.list)

#find the anchors
ADT.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features)

#do the integration
adt.combined <- IntegrateData(anchorset = ADT.anchors)

#now the adt.combined has 3 assays too:RNA, Protein and Integrated but integrated assay here has the common protein info of the samples NOT the RNA
# Run the standard workflow 
DefaultAssay(adt.combined) <- "integrated"
adt.combined <- ScaleData(adt.combined, verbose = FALSE)
adt.combined <- RunPCA(adt.combined,reduction.name = 'apca')

DefaultAssay(rna.combined) <- "integrated"
DefaultAssay(adt.combined) <- "integrated"
#making matrix based on integrated adt
adt.data <- GetAssayData(object =  adt.combined[['integrated']], slot = 'data')
dim(adt.data)
adt.scaled <- GetAssayData(object =  adt.combined[['integrated']], slot = 'scale.data')
dim(adt.scaled)

#add the integrated adt data as an another assay in rna.combined
rna.combined[["integrated.adt"]] <- CreateAssayObject(data = adt.data )
rna.combined[['adt.pca']] <- adt.combined[['apca']]

#set the used assay in rna.combined for integrated.adt part as adt.pca 
rna.combined[['adt.pca']]@assay.used <-  "integrated.adt"
rna.combined[['adt.pca']]@assay.used

#set the default assay of rna.combined to integrated.adt to centre and scale the matrix
DefaultAssay(rna.combined) <- "integrated.adt"
rna.combined <- ScaleData(rna.combined)

#check the data and reduction method match
rna.combined[['pca']]@assay.used
rna.combined[['adt.pca']]@assay.used

#save the rna.combined 
saveRDS(rna.combined, paste0("./batch7_outputs_filtered/",'rna.combined.rds'))

#save the adt.combined 
saveRDS(adt.combined, paste0("./batch7_outputs_filtered/",'adt.combined.rds'))


###Task 3: WNN Integration, Clustering, and Annotation*
#run wnn on rna.combined and save it as new variable or immune.combined
immune.combined  <- FindMultiModalNeighbors(
  rna.combined , reduction.list = list("pca", "adt.pca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)
#Run umap on the wnn-performed data
immune.combined  <- RunUMAP(immune.combined , nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

#Run clustring
immune.combined  <- FindClusters(immune.combined , graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)

#Plot the clusters
pdf(paste0("./batch7_outputs_filtered/","wnn_umap.pdf"),width = 20,height = 10)
p1d <- DimPlot(immune.combined , reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 8) 
p2d <- DimPlot(immune.combined , reduction = 'wnn.umap', group.by = 'orig.ident',label = TRUE, repel = TRUE, label.size = 8) 
p1d 
p2d
dev.off()

DefaultAssay(immune.combined) <- 'integrated'
markers_RNA <- FindAllMarkers(immune.combined , only.pos = TRUE)
write.xlsx(markers_RNA,file = paste0("./batch7_outputs_filtered/IntegratedRNA_markers.xlsx"),row.names = TRUE)

# get the 20 first markers of RNA integrated
top20.markers_RNA <- markers_RNA %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top20 <- c()
for(c in unique(top20.markers_RNA$cluster)){
  tmp <- top20.markers_RNA[top20.markers_RNA$cluster==c,]$gene
  top20 <- cbind(top20,tmp)
}
colnames(top20) <- paste0('C_',unique(top20.markers_RNA$cluster))
write.csv(top20,"./batch7_outputs_filtered/top20.IntegratedRNA_markers.csv")

#now set the default set of immune.combined to integrated protein data or integrated.adt to find the gene markers of those

DefaultAssay(immune.combined) <- 'integrated.adt'
markers_ADT <- FindAllMarkers(immune.combined, only.pos = TRUE,verbose = F)
write.xlsx(markers_ADT,file = paste0("./batch7_outputs_filtered/IntegratedADT_markers.xlsx"),row.names = TRUE)

#get the 20 first markers of protein integrated
top20.markers_ADT <- markers_ADT %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top20 <- c()
for(c in unique(top20.markers_ADT$cluster)){
  tmp <- top20.markers_ADT[top20.markers_ADT$cluster==c,]$gene
  top20 <- cbind(top20,tmp)
}
colnames(top20) <- paste0('C_',unique(top20.markers_ADT$cluster))
write.csv(top20,"./batch7_outputs_filtered/top20.IntegratedADT_markers.csv")

#run umap for pca
immune.combined  <- RunUMAP(immune.combined , reduction = 'pca', dims = 1:30, assay = 'integrated', 
                            reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

#run umap for apca
immune.combined  <- RunUMAP(immune.combined , reduction = 'adt.pca', dims = 1:10, assay = 'integrated.adt', 
                            reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')


pdf(paste0("./batch7_outputs_filtered/","RNA_ADT_umap.pdf"),width = 20,height = 10)
p1d <- DimPlot(immune.combined , reduction = 'rna.umap', label = TRUE, repel = TRUE, label.size = 8)
p2d <- DimPlot(immune.combined , reduction = 'adt.umap', label = TRUE, repel = TRUE, label.size = 8)
p1d 
p2d
dev.off()


pdf(paste0('./batch7_outputs_filtered/canonical_dotplot.pdf'),width = 20,height = 10)
library(ggplot2)
library(RColorBrewer)
library(angrycell)
DefaultAssay(immune.combined) <- "integrated"
angrycell::DotPlot2(immune.combined, features = c("CD4","CD3D", "FOXP3","IL7R", "CD8A", "CD8B","S100A9","S100A8","CSF3R","CSF1R","LILRA4", "CD68", "KLRG1", "CD19", "JCHAIN", "CD79A", "MS4A1","BLK","KIT", "HBB","HBA1","HBA2","MKI67"), dot.scale = 8) +
  coord_flip() +
  theme(text = element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y = element_text(size =10),
        axis.title.x = element_text(size=10)) 
dev.off()

#save the immune.combined 
saveRDS(immune.combined, paste0("./batch7_outputs_filtered/",'immune.combined.rds'))

#add annotation
Idents(immune.combined) <- immune.combined@meta.data$seurat_clusters
immune.combined <- RenameIdents(immune.combined, `0` = "Myeloid", `1` = "CD4", `2` = "B",
                                `3` = "NK", `4` = "CD8", `5` = "Treg", `6` = "NK", 
                                `7` = "NK",`8` = "Myeloid", `9` = "Myeloid", `10`= "CD8", `11` = "CD8", `12`= "NK", `13`= "T", `14`= "NK", `15` = "CD8", `16`= "CD8", `17` ="NK", `18` = "CD4", `19` = "Myeloid", `20` = "Myeloid", `21`=  "Myeloid", `22`= "Myeloid", `23`= "Plasma", `24` = "RBC", `25`="Proliferation_NK", `26`= "Proliferation_T", `27`= "pDC", `28`="B", `29`="Mast", `30`= "Neutrophil")


immune.combined@meta.data$Annotation<- Idents(immune.combined)

pdf(paste0('./batch7_outputs_filtered/annotation.pdf'))
DimPlot(immune.combined, reduction = "wnn.umap", label = TRUE, repel = TRUE)
dev.off()

#save the immune.combined 
saveRDS(immune.combined, paste0("./batch7_outputs_filtered/",'immune.combined.rds'))




###Task 4: Gene Feature Extraction and Visualization
#function to get the gene expression from the matrix and add it to the metdata of seurat
rna_T = Matrix::t(tumor_v2@assays$RNA@counts)
meta_T <- tumor_v2@meta.data

gitr = function( mylist, d = rna_T, m = meta_T ) {
  cbind(
    m,
    as.matrix(d[, colnames(d) %in% mylist, drop = FALSE])
  )
}

genelist <- c('CXCl13', 'CCl8', 'KRT1', 'Krt5')
new_meta = gitr(genelist)


new_meta$gene.pos.Ccl8[new_meta[,'CCL8'] > 0] <- 1
new_meta$gene.pos.Ccl8[new_meta[,'CCL8'] ==0] <- 0

#plot1    
library(tidyverse)
gene.sum = function(gene, gene.positive,  data=new_meta) {
  data%>%
    #both of the 2 following commands work, calling the column directly or call it through.data$
    group_by(CellMajorType, Group) %>%
    #group_by(.data$CellMajorType, .data$Group) %>%
    summarize(gene.pos = round( 100 *  sum({{gene.positive}}) / n(), digits = 2), gene.counts = sum({{gene}})) %>%  
    ggplot(aes( x =Group, y = gene.counts, size=gene.pos)) +
    #geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = Group, y = gene.counts)) +
    geom_point(color = 'darkblue', aes(alpha=gene.pos)) + coord_flip() + facet_wrap(scales = "free_x", "CellMajorType") +
    theme(axis.text.x = element_text(size=8, angle=45))+
    theme_hc()  
  
  
}


#plot2
library(tidyverse)

gene.sum = function(gene, gene.positive,  data=new_meta) {
  data %>%
    #both of the 2 following commands work, calling the column directory or call it through.data$
    group_by(CellMajorType, hash.ID, Group) %>%
    #group_by(.data$CellMajorType, .data$hash.ID, .data$Group) %>%
    summarize(gene.pos = round( 100 *  sum({{gene.positive}}) / n(), digits = 2), gene.counts = sum({{gene}})) %>%  
    
    ggplot(aes( x =Group, y = log2(1 + gene.counts) )) +
    geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = Group, y = log2(1 +gene.counts) ) ) +
    geom_point(position = position_jitterdodge(), color = 'darkblue', aes( size = gene.pos) )+ coord_flip() + facet_wrap(scales = "fixed", "CellMajorType") +
    theme(axis.text.x = element_text(size=8, angle=45))+
    theme_hc()  
  
  
}


gene.sum(CCL8, gene.pos.CCL8)
