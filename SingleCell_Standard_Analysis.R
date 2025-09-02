
#libraries
library(Seurat)
suppressMessages(require(cowplot))
suppressMessages(require(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(ggplot2))
library(ggplot2)
library(RColorBrewer)
library(angrycell)

#make seurat object for each, add hashtag data as an assay object
dir.create("/home/nasim/amgen.abio/amgen_outputs_with_negatives")
ids <- c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6", "Group7", "Group8")

for (file in ids){
  data <- Read10X(data.dir = paste0("/home/nasim/amgen.abio/", file))
  seurat_obj <- CreateSeuratObject(counts = data$`Gene Expression`,project = file)
  seurat_obj[['HTO']] = CreateAssayObject(counts = data$`Antibody Capture`)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = '^mt-')
  seurat_obj <- NormalizeData(seurat_obj, assay="HTO",normalization.method = "CLR")
  #change the lower threshod quantile=70 for demultiplexing
  seurat_obj <- MULTIseqDemux(seurat_obj,assay="HTO",quantile=0.99)
  write.csv(table(seurat_obj@meta.data$MULTI_ID, seurat_obj@meta.data$orig.ident),file = paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/', file,'_','Demultiplex.csv'),row.names = TRUE)
  write.csv(table(seurat_obj@meta.data$MULTI_classification,seurat_obj@meta.data$orig.ident),file=paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/',file,'_','Demultiplex_classification.csv'),row.names = TRUE)
  assign(file, seurat_obj)
}

#merge
combined <- merge(Group1, c(Group2, Group3, Group4, Group5, Group6, Group7, Group8), add.cell.ids = c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6", "Group7", "Group8"), project = "amgen")
combined$Group<- combined$orig.ident

#stats before filtering
pdf(paste0("/home/nasim/amgen.abio/amgen_outputs_with_negatives/Beforefilter_violin.pdf"),width = 20,height = 10)
tf <- function(x){tapply(x, Idents(combined), median)}
Idents(combined) <- combined$orig.ident
median.value <- apply(combined@meta.data[,c('nCount_RNA','nFeature_RNA','nCount_HTO','percent.mt')],2,tf)
write.csv(median.value,'/home/nasim/amgen.abio/amgen_outputs_with_negatives/Beforefilter.median.csv')
saveRDS(combined, '/home/nasim/amgen.abio/amgen_outputs_with_negatives/0-Data_RawObject.rds')
VlnPlot(combined, features = c('nCount_RNA','nFeature_RNA','nCount_HTO','percent.mt'),group.by='orig.ident', pt.size = 0,ncol = 3)
write.csv(table(combined@meta.data$Group),'/home/nasim/amgen.abio/amgen_outputs_with_negatives/Beforefilter.cellnumber_per_Group.csv')
dev.off()

#stats after filtering, keep the negatives
combined.filter <- subset(combined, subset = nFeature_RNA >= 200 & percent.mt < 10)
Idents(combined.filter) <- combined.filter@meta.data$MULTI_ID
#singlet <- subset(combined.filter, idents = "Negative", invert = TRUE)
singlet <- subset(combined.filter, idents = "Doublet", invert = TRUE)

Idents(singlet) <- singlet$orig.ident
pdf(paste0("/home/nasim/amgen.abio/amgen_outputs_with_negatives/Afterfilter_violin.pdf"),width = 20,height = 10)
tf1 <- function(x){tapply(x, Idents(singlet), median)}

median.value <- apply(singlet@meta.data[,c('nCount_RNA','nFeature_RNA','nCount_HTO','percent.mt')],2,tf1)
write.csv(median.value,'/home/nasim/amgen.abio/amgen_outputs_with_negatives/Afterfilter.median.csv')
saveRDS(singlet, '/home/nasim/amgen.abio/amgen_outputs_with_negatives/1-Data_QualityControl.rds')
VlnPlot(singlet, features = c('nCount_RNA','nFeature_RNA','nCount_HTO','percent.mt'),group.by='orig.ident', pt.size = 0,ncol = 3)
write.csv(table(singlet@meta.data$Group),'/home/nasim/amgen.abio/amgen_outputs_with_negatives/Afterfilter.cellnumber_per_Group.csv')
dev.off()


#add metadata
singlet$Timepoint<- ""
singlet$Timepoint[singlet@meta.data$orig.ident=="Group1"] <- "Day14"
singlet$Timepoint[singlet@meta.data$orig.ident=="Group2"] <- "Day14"
singlet$Timepoint[singlet@meta.data$orig.ident=="Group3"] <- "Day14"
singlet$Timepoint[singlet@meta.data$orig.ident=="Group4"] <- "Day14"

singlet$Timepoint[singlet@meta.data$orig.ident=="Group5"] <- "Day21"
singlet$Timepoint[singlet@meta.data$orig.ident=="Group6"] <- "Day21"
singlet$Timepoint[singlet@meta.data$orig.ident=="Group7"] <- "Day21"
singlet$Timepoint[singlet@meta.data$orig.ident=="Group8"] <- "Day21"

singlet$Treatment<- ""
singlet$Treatment[singlet@meta.data$orig.ident=="Group1"] <- "C"
singlet$Treatment[singlet@meta.data$orig.ident=="Group2"] <- "B"
singlet$Treatment[singlet@meta.data$orig.ident=="Group3"] <- "B+T1"
singlet$Treatment[singlet@meta.data$orig.ident=="Group4"] <- "B+T2"

singlet$Treatment[singlet@meta.data$orig.ident=="Group5"] <- "C"
singlet$Treatment[singlet@meta.data$orig.ident=="Group6"] <- "B"
singlet$Treatment[singlet@meta.data$orig.ident=="Group7"] <- "B+T1"
singlet$Treatment[singlet@meta.data$orig.ident=="Group8"] <- "B+T2"

#preprocessing
DefaultAssay(singlet) <- "RNA"
singlet <-NormalizeData(singlet, normalization.method = "LogNormalize", scale.factor = 10000)
singlet <- FindVariableFeatures(singlet, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(singlet), 10)
pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/variable_genes.pdf')
plot1 <- VariableFeaturePlot(singlet)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#print(plot1)
print(plot2)
dev.off()

#Scale and RunPCA, PCA realated plots
pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/PCA_plots.pdf')
singlet <- ScaleData(singlet)
singlet <- RunPCA(singlet, features = VariableFeatures(object = singlet))
DimPlot(singlet,group.by = 'Group',label=T,repel=T)
DimHeatmap(singlet, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# #Elbow plot
pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Elbow_plot.pdf')
P1 <- ElbowPlot(singlet)
print(P1)
dev.off()

#make pca plot for each group

pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/PCA_plots_per_Group.pdf')
for (i in unique(x = singlet@meta.data$Group)) {
  subset <- subset(singlet, subset=Group==i)
  print(DimPlot(subset, reduction= "pca", group.by = 'MULTI_ID',label=F,repel=T)+ ggtitle(i))
}
dev.off()

#Elbow plot with percentage variance shown
# Extract the standard deviations
stdev <- singlet[["pca"]]@stdev

# Calculate variance explained for each PC
variance_explained <- (stdev^2) / sum(stdev^2) * 100

# Base R plot
pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Elbow_plot_percentage_variance.pdf',width = 20,height = 10)

plot(variance_explained, type = "b", xlab = "Principal Component", ylab = "Percentage of Variance Explained",
     main = "Elbow Plot", pch = 19, col = "blue", cex.lab = 1.5)

# Add labels to each point
text(x = 1:length(variance_explained), y = variance_explained,
     labels = paste0(round(variance_explained, 1), "%"), pos = 3, cex = 0.8, col = "darkred", size=2)
dev.off()

#run clustering 
singlet <- FindNeighbors(singlet, reduction = "pca",dims = 1:20)
singlet <- FindClusters(singlet, resolution = 0.5)
singlet <- RunTSNE(singlet, dims = 1:20)
singlet <- RunUMAP(singlet, dims = 1:20)

#plots
pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/umap_tsne_plots.pdf', width = 12,height = 10)
DimPlot(singlet,group.by = 'seurat_clusters',reduction = 'umap',label=T,repel=T)
DimPlot(singlet,group.by = 'seurat_clusters',reduction = 'tsne',label=T,repel=T)
DimPlot(singlet,group.by = 'Group',reduction = 'umap')
DimPlot(singlet,group.by = 'Group',reduction = 'tsne')
DimPlot(singlet,group.by = 'Timepoint',reduction = 'umap')
DimPlot(singlet,group.by = 'Timepoint',reduction = 'tsne')
DimPlot(singlet,group.by = 'Treatment',reduction = 'umap')
DimPlot(singlet,group.by = 'Treatment',reduction = 'tsne')

DimPlot(singlet,split.by ="Timepoint", group.by = 'Treatment',reduction = 'umap')
DimPlot(singlet,split.by ="Timepoint", group.by = 'Treatment',reduction = 'tsne')
DimPlot(singlet,split.by ="Treatment", group.by = 'Timepoint',reduction = 'umap')
DimPlot(singlet,split.by ="Treatment", group.by = 'Timepoint',reduction = 'tsne')

dev.off()


#Find gene markers for each cluster
markers<- FindAllMarkers(singlet , only.pos = TRUE)
write.xlsx(markers,file = paste0("/home/nasim/amgen.abio/amgen_outputs_with_negatives/RNA_markers.xlsx"),row.names = TRUE)

#top20
top20.markers_RNA <- markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top20 <- c()
for(c in unique(top20.markers_RNA$cluster)){
  tmp <- top20.markers_RNA[top20.markers_RNA$cluster==c,]$gene
  top20 <- cbind(top20,tmp)
}
colnames(top20) <- paste0('C_',unique(top20.markers_RNA$cluster))
write.csv(top20,"/home/nasim/amgen.abio/amgen_outputs_with_negatives/top20.RNA_markers.csv")


#save
saveRDS(singlet, '/home/nasim/amgen.abio/amgen_outputs_with_negatives/2-Data_DimR.rds')


#feature and dot plot
pdf(paste0("/home/nasim/amgen.abio/amgen_outputs_with_negatives/","canonical_markers_featureplot.pdf"),width = 10,height = 10)
DefaultAssay(singlet) <- "RNA"
#Naive/ memory T
P1=FeaturePlot(singlet , features = c("Tcf7", "Lef1", "Ccr7", "Sell", "Mal", "Cd3", "Cd4", "Cd8", "Cd11a", "Cd14", "Cd19", "Cd25", "Cd27", "Cd28", "Cd44", "Cd45ra", "Cd45ro", "Cd57", "Cd621", "Cd122", "Cd127", "Tfn-y", "Il-2", "Tnf", "Il7r"),  ncol = 3)

P2=FeaturePlot(singlet , features = c("Cd3e","Cd3d", "Cd2", "Ccr1", "Ccr5", "Cxcr3", "Stat4", "T-bet", "Ccr3", "Ccr4", "Ccr8"," Cxcr4", "Stat5", "Stat6", "Gata-3","Il9", "Ccr6", "roryT", "Ccr10", "Fgf", "Ahr9", "Btla", "Bcl-6", "Stat","Il2ra", "Foxp3"),  ncol = 3)

P3=FeaturePlot(singlet, features = c("Cd8b1","Cd8a", "Ca8b", "Cx3cr1", "Gzmh", "Tbx21", "Eomes","Zfp36l2", "Gzmk", "Gzma", "Ccl5", "Infg", "Runx3",  "Cd69", "Itgae"),  ncol = 3)

P4=FeaturePlot(singlet, features = c("Klrd1", "Tyrobp", "Nkg7", "Ncam1", "Klra3", "Klra7", "Klrg1", "Kir2dl1", "Kir2dl3", "Kir2dl4", "Kir3dl1", "Kir3dl2", "Kir3dl3", "Kir2ds4"),  ncol = 3)

P5=FeaturePlot(singlet, features = c( "Ifit1", "Ifit2", "Ifit3"),  ncol = 3)

P6=FeaturePlot(singlet, features = c( "Cd79a", "Ms4a1", "Blk", "BankK1", "Jchain", "Hbb","Hba1","Hba2"),  ncol = 3)

P7=FeaturePlot(singlet, features = c( "Cd86","Cd80","Cd68", "Cd163", "Lyz", "Ms4a4a", "Ly6c2", "S100a9", "S100a8", "Kit","Hbb-bs", "Hba-a1","Slamf7", "Mzb1", "Irf7",  "Csf3r", "Ms4a2", "Csf1r", "C3ar1","Lilra4",  "H2-aa", "H2-ab1"),  ncol = 3)

P8=FeaturePlot(singlet, features = c("Ceacam8", "Arg1", "Flt1", "Icam1", "Il17ra", "Cxcr4"),  ncol = 3)

P9=FeaturePlot(singlet, features = c("Gng11", "Cavin2", "Tubb1"))

P10=FeaturePlot(singlet, features = c( "Mki67", "Stmn1"),  ncol = 2)

P11=FeaturePlot(singlet, features = c( "Pecam1", "Epcam", "Dcn", "Col1a1", "Krt8", "Krt18"),  ncol = 3)
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
dev.off()

pdf(paste0("/home/nasim/amgen.abio/amgen_outputs_with_negatives/","canonical_markers_dotplot.pdf"),width = 10,height = 10)
P12=angrycell::DotPlot2(singlet, features = c("Tcf7", "Lef1", "Ccr7", "Sell", "Mal", "Cd3", "Cd4", "Cd8", "Cd11a", "Cd14", "Cd19", "Cd25", "Cd27", "Cd28", "Cd44", "Cd45ra", "Cd45ro", "Cd57", "Cd621", "Cd122", "Cd127", "Tfn-y", "Il-2", "Tnf", "Il7r", "Cd3e","Cd3d", "Cd2", "Ccr1", "Ccr5", "Cxcr3", "Stat4", "T-bet", "Ccr3", "Ccr4", "Ccr8"," Cxcr4", "Stat5", "Stat6", "Gata-3","Il9", "Ccr6", "roryT", "Ccr10", "Fgf", "Ahr9", "Btla", "Bcl-6", "Stat","Il2ra", "Foxp3","Cd8b1","Cd8a", "Ca8b", "Cx3cr1", "Gzmh", "Tbx21", "Eomes","Zfp36l2", "Gzmk", "Gzma", "Ccl5", "Infg", "Runx3",  "Cd69", "Itgae", "Klrd1", "Tyrobp", "Nkg7", "Ncam1", "Klra3", "Klra7", "Klrg1", "Kir2dl1", "Kir2dl3", "Kir2dl4", "Kir3dl1", "Kir3dl2", "Kir3dl3", "Kir2ds4","Ifit1", "Ifit2", "Ifit3", "Cd79a", "Ms4a1", "Blk", "BankK1", "Jchain", "Hbb","Hba1","Hba2", "Cd86","Cd80","Cd68", "Cd163", "Lyz", "Ms4a4a", "Ly6c2", "S100a9", "S100a8", "Kit","Hbb-bs", "Hba-a1","Slamf7", "Mzb1", "Irf7",  "Csf3r", "Ms4a2", "Csf1r", "C3ar1","Lilra4",  "H2-aa", "H2-ab1","Ceacam8", "Arg1", "Flt1", "Icam1", "Il17ra", "Cxcr4", "Gng11", "Cavin2", "Tubb1", "Mki67", "Stmn1", "Pecam1", "Epcam", "Dcn", "Col1a1", "Krt8", "Krt18" ), dot.scale = 8) +
  coord_flip() +
  theme(text = element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=8),
        axis.title.y = element_text(size =10),
        axis.title.x = element_text(size=10))

print(P12)
dev.off()



#saving differnt dataframes
library(openxlsx)
counts <- table(singlet@meta.data$seurat_clusters)

groups <-table(singlet@meta.data$seurat_clusters, singlet@meta.data$Group)

timepoint <-table(singlet@meta.data$seurat_clusters, singlet@meta.data$Timepoint)

treatment  <-table(singlet@meta.data$seurat_clusters, singlet@meta.data$Treatment)

metadata <- singlet@meta.data
cell_counts <- metadata %>%
  group_by(seurat_clusters, Timepoint, Treatment) %>%
  summarise(num_cells = n(), .groups = "drop")

cell_counts_df <- as.data.frame(cell_counts)

# Create a new workbook
wb <- createWorkbook()

# Add worksheets to the workbook
addWorksheet(wb, "Sheet1")
addWorksheet(wb, "Sheet2")
addWorksheet(wb, "Sheet3")
addWorksheet(wb, "Sheet4")
addWorksheet(wb, "Sheet5")


# Write the dataframes to the respective sheets
writeData(wb, sheet = "Sheet1", x = counts )
writeData(wb, sheet = "Sheet2", x = groups)
writeData(wb, sheet = "Sheet3", x = timepoint)
writeData(wb, sheet = "Sheet4", x = treatment)
writeData(wb, sheet = "Sheet5", x = cell_counts_df)

# Save the workbook to a file
saveWorkbook(wb, file = "/home/nasim/amgen.abio/amgen_outputs_with_negatives/cellcount_per_cluster.xlsx", overwrite = TRUE)


#add anno
Idents(singlet) <- singlet@meta.data$seurat_clusters
singlet <- RenameIdents(singlet, `0` = "NK", `1` = "Naïve T", `2` = "CD8", `3` = "Proliferative cells", `4` = "CD4", `5` ="CD8", `6` = "Treg", `7` = "Proliferative cells",`8` = "CD8", `9` = "CD4", `10`= "CD4", `11` = "CD8", `12`= "Monocytes", `13`= "CD4", `14`= "CD4", `15`= "Neutrophils" , `16`= "Mast", `17`="Monocytes", `18`="CD8", `19`="Fibroblast", `20`="Proliferative cells" )

singlet@meta.data$Annotation<- Idents(singlet)

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/annotation.pdf'),width = 10,height = 10)
DimPlot(singlet, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

#save
saveRDS(singlet, '/home/nasim/amgen.abio/amgen_outputs_with_negatives/3-Data_Annotation.rds')


#remove some clusters
singlet_filter <- subset(singlet, subset=seurat_clusters==0 |seurat_clusters==1| seurat_clusters==2 |seurat_clusters==3 | seurat_clusters==4 |seurat_clusters==5| seurat_clusters==6 | seurat_clusters==7 | seurat_clusters==8 | seurat_clusters==9 | seurat_clusters==10 | seurat_clusters==11 | seurat_clusters==13 | seurat_clusters==14 | seurat_clusters==18 | seurat_clusters==20 )


#make the plots for remaining filetered data
singlet_filter <- FindVariableFeatures(singlet_filter, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(singlet_filter), 10)
pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/variable_genes_filetered_data.pdf')
plot1 <- VariableFeaturePlot(singlet_filter)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#print(plot1)
print(plot2)
dev.off()

pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/PCA_plots_filetered_data.pdf')
singlet_filter  <- RunPCA(singlet_filter, features = VariableFeatures(object = singlet_filter))
DimPlot(singlet_filter,reduction= "pca",group.by = 'Group',label=F,repel=T)
DimHeatmap(singlet_filter, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/PCA_plots_per_Group_filetered_data.pdf')
for (i in unique(x = singlet_filter@meta.data$Group)) {
  subset <- subset(singlet_filter, subset=Group==i)
  print(DimPlot(subset, reduction= "pca", group.by = 'MULTI_ID',label=F,repel=T)+ ggtitle(i))
}
dev.off()

# Elbow plot
pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Elbow_plot_filetered_data.pdf')
P1 <- ElbowPlot(singlet_filter)
print(P1)
dev.off()

#Elbow plot with percentage variance shown

# Extract the standard deviations
stdev <- singlet_filter[["pca"]]@stdev

# Calculate variance explained for each PC
variance_explained <- (stdev^2) / sum(stdev^2) * 100

# Base R plot
pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Elbow_plot_percentage_variance_filetered_data.pdf',width = 20,height = 10)

plot(variance_explained, type = "b", xlab = "Principal Component", ylab = "Percentage of Variance Explained",
     main = "Elbow Plot", pch = 19, col = "blue", cex.lab = 1.5)

# Add labels to each point
text(x = 1:length(variance_explained), y = variance_explained,
     labels = paste0(round(variance_explained, 1), "%"), pos = 3, cex = 0.8, col = "darkred", size=2)
dev.off()

#run clustering 
singlet_filter <- FindNeighbors(singlet_filter, reduction = "pca",dims = 1:20)
singlet_filter <- FindClusters(singlet_filter, resolution = 0.5)
singlet_filter <- RunTSNE(singlet_filter, dims = 1:20)
singlet_filter <- RunUMAP(singlet_filter, dims = 1:20)


pdf('/home/nasim/amgen.abio/amgen_outputs_with_negatives/umap_tsne_plots_filetered_data.pdf', width = 12,height = 10)
DimPlot(singlet_filter,group.by = 'seurat_clusters',reduction = 'umap',label=T,repel=T)
DimPlot(singlet_filter,group.by = 'seurat_clusters',reduction = 'tsne',label=T,repel=T)
DimPlot(singlet_filter,group.by = 'Group',reduction = 'umap')
DimPlot(singlet_filter,group.by = 'Group',reduction = 'tsne')
DimPlot(singlet_filter,group.by = 'Timepoint',reduction = 'umap')
DimPlot(singlet_filter,group.by = 'Timepoint',reduction = 'tsne')
DimPlot(singlet_filter,group.by = 'Treatment',reduction = 'umap')
DimPlot(singlet_filter,group.by = 'Treatment',reduction = 'tsne')

DimPlot(singlet_filter,split.by ="Timepoint", group.by = 'Treatment',reduction = 'umap')
DimPlot(singlet_filter,split.by ="Timepoint", group.by = 'Treatment',reduction = 'tsne')
DimPlot(singlet_filter,split.by ="Treatment", group.by = 'Timepoint',reduction = 'umap')
DimPlot(singlet_filter,split.by ="Treatment", group.by = 'Timepoint',reduction = 'tsne')
dev.off()


#Find gene markers for each cluster

markers<- FindAllMarkers(singlet_filter, only.pos = TRUE)
write.xlsx(markers,file = paste0("/home/nasim/amgen.abio/amgen_outputs_with_negatives/RNA_markers_filetered_data.xlsx"),row.names = TRUE)
#top20
top20.markers_RNA <- markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top20 <- c()
for(c in unique(top20.markers_RNA$cluster)){
  tmp <- top20.markers_RNA[top20.markers_RNA$cluster==c,]$gene
  top20 <- cbind(top20,tmp)
}
colnames(top20) <- paste0('C_',unique(top20.markers_RNA$cluster))
write.csv(top20,"/home/nasim/amgen.abio/amgen_outputs_with_negatives/top20.RNA_markers_filetered_data.csv")

#gene markers
pdf(paste0("/home/nasim/amgen.abio/amgen_outputs_with_negatives/","canonical_markers_featureplot_filtered_data.pdf"),width = 10,height = 10)
P1=FeaturePlot(singlet_filter , features = c("Tcf7", "Lef1", "Ccr7", "Sell", "Mal", "Cd3", "Cd4", "Cd8", "Cd11a", "Cd14", "Cd19", "Cd25", "Cd27", "Cd28", "Cd44", "Cd45ra", "Cd45ro", "Cd57", "Cd621", "Cd122", "Cd127", "Tfn-y", "Il-2", "Tnf", "Il7r"),  ncol = 3)

P2=FeaturePlot(singlet_filter , features = c("Cd3e","Cd3d", "Cd2", "Ccr1", "Ccr5", "Cxcr3", "Stat4", "T-bet", "Ccr3", "Ccr4", "Ccr8"," Cxcr4", "Stat5", "Stat6", "Gata-3", "Ccr6", "roryT", "Ccr10", "Fgf", "Ahr9", "Btla", "Bcl-6", "Stat","Il2ra", "Foxp3","Blk"),  ncol = 3)

P3=FeaturePlot(singlet_filter , features = c("Cd8b1","Cd8a", "Ca8b", "Cx3cr1", "Gzmh", "Tbx21", "Eomes","Zfp36l2", "Gzmk", "Gzma", "Ccl5", "Infg", "Runx3",  "Cd69", "Itgae"),  ncol = 3)

P4=FeaturePlot(singlet_filter , features = c("Klrd1", "Tyrobp", "Nkg7", "Ncam1", "Klra3", "Klra7", "Klrg1", "Kir2dl1", "Kir2dl3", "Kir2dl4", "Kir3dl1", "Kir3dl2", "Kir3dl3", "Kir2ds4"),  ncol = 3)

P5=FeaturePlot(singlet_filter , features = c("Ifit1", "Ifit2", "Ifit3"),  ncol = 3)

P6=FeaturePlot(singlet_filter , features = c("Ctla4", "Icos", "Tigit", "Pdcd1"),  ncol = 3)

P7=FeaturePlot(singlet_filter, features = c("Mki67", "Stmn1"),  ncol = 2)
print(P1)
print(P2)
print(P3)
print(P4)
print(P5)
print(P6)
print(P7)
dev.off()

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/canonical_markers_dotplot_filtered_data.pdf'), width = 15,height = 20)
library(ggplot2)
library(RColorBrewer)
library(angrycell)
Idents(singlet_filter) <- singlet_filter@meta.data$seurat_clusters
angrycell::DotPlot2(singlet_filter, features = c("Tcf7", "Lef1", "Ccr7", "Sell", "Mal", "Cd3", "Cd4", "Cd8", "Cd11a", "Cd14", "Cd19", "Cd25", "Cd27", "Cd28", "Cd44", "Cd45ra", "Cd45ro", "Cd57", "Cd621", "Cd122", "Cd127", "Tfn-y", "Il-2", "Tnf", "Il7r", "Cd3e","Cd3d", "Cd2", "Ccr1", "Ccr5", "Cxcr3", "Stat4", "T-bet", "Ccr3", "Ccr4", "Ccr8"," Cxcr4", "Stat5", "Stat6", "Gata-3", "Ccr6", "roryT", "Ccr10", "Fgf", "Ahr9", "Btla", "Bcl-6", "Stat","Il2ra", "Foxp3","Blk","Cd8b1","Cd8a", "Ca8b", "Cx3cr1", "Gzmh", "Tbx21", "Eomes","Zfp36l2", "Gzmk", "Gzma", "Ccl5", "Infg", "Runx3",  "Cd69", "Itgae", "Klrd1", "Tyrobp", "Nkg7", "Ncam1", "Klra3", "Klra7", "Klrg1", "Kir2dl1", "Kir2dl3", "Kir2dl4", "Kir3dl1", "Kir3dl2", "Kir3dl3", "Kir2ds4","Ifit1", "Ifit2", "Ifit3","Ctla4", "Icos", "Tigit", "Pdcd1","Mki67", "Stmn1"), dot.scale = 8) +
  coord_flip() +
  theme(text = element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.title.y = element_text(size =11),
        axis.title.x = element_text(size=12)) 
dev.off()

#add anno
Idents(singlet_filter) <- singlet_filter@meta.data$seurat_clusters
singlet_filter <- RenameIdents(singlet_filter, `0` = "Exhausted CD8", `1` = "Stem-like cells", `2` = "Naïve T", `3` = "Activated CD8", `4` = "CD4_Ccr8", `5` ="Stem-like cells", `6` = "Activated CD8", `7` = "CD4_Treg",`8` = "Activated CD8", `9` = "Proliferative cells", `10`= "CD8", `11` = "NK", `12`= "CD4_Treg", `13`= "B-like cells", `14`= "CD8_Ifit", `15`= "B-like cells" , `16`= "CD8", `17`="Activated CD8" )

singlet_filter@meta.data$Annotation<- Idents(singlet_filter)

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/annotation_filtered_data.pdf'),width = 10,height = 10)
DimPlot(singlet_filter, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

#change the annotation
Idents(singlet_filter) <- singlet_filter@meta.data$seurat_clusters
singlet_filter <- RenameIdents(singlet_filter, `0` = "C0_CD8_Fcrl6_Texh", `1` = "C1_Sox5_Stem-like", `2` = "C2_Id3_Sell_Naïve", `3` = "C3_CD8_Ki67_Stmn1", `4` = "C4_CD4_Ccr8_Th2", `5` ="C5_CD4_Pparg_CD30l", `6` = "C6_CD8_Batf3_mcm", `7` = "C7_Foxp3_Ox40_Il10_Treg",`8` = "C8_CD8_Ki67", `9` = "C9_Mi67_Runx3_Proliferating", `10`= "C10_CD8_Klra10", `11` = "C11_NK", `12`= "C12_Foxp3_Itgb8_Treg", `13`= "C13_Rorc_Il17a_Th17", `14`= "C14_CD8_Ifit", `15`= "C15_Il23r_Th17" , `16`= "C16_CD8_Xcl1_Irf8", `17`="C17_Foxp3_Ki67_Treg" )

singlet_filter@meta.data$Annotation<- Idents(singlet_filter)

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/annotation_filtered_data.pdf'),width = 10,height = 10)
DimPlot(singlet_filter, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(singlet_filter,split.by ="Timepoint", group.by = 'Annotation',reduction = 'umap')
DimPlot(singlet_filter,split.by ="Treatment", group.by = 'Annotation',reduction = 'umap')
dev.off()


#cell numbers per timepoint per treatment

metadata <- singlet_filter@meta.data
cell_counts <- metadata %>%
  group_by(Annotation, Group, Timepoint, Treatment, MULTI_classification) %>%
  summarise(num_cells = n(), .groups = "drop")

cell_counts_df <- as.data.frame(cell_counts)
write.csv(cell_counts_df,'/home/nasim/amgen.abio/amgen_outputs_with_negatives/cellnumbers_split_categories.csv', row.names = FALSE )


#change the annotation, braoder annotation, aggregated cluster
Idents(singlet_filter) <- singlet_filter@meta.data$seurat_clusters
singlet_filter <- RenameIdents(singlet_filter, `0` = "CD8_T-cells", `1` = "C1_Sox5_Stem-like", `2` = "C2_Id3_Sell_Naive", `3` = "CD8_T-cells", `4` = "CD4_Th2", `5` ="CD4_Th2", `6` = "CD8_T-cells", `7`="T-regs_Foxp3",`8` = "CD8_T-cells", `9` = "Proliferating_Mki67", `10`= "CD8_T-cells", `11` = "C11_NK", `12`= "T-regs_Foxp3", `13`= "Th17", `14`= "CD8_T-cells", `15`= "Th17" , `16`= "CD8_T-cells", `17`="T-regs_Foxp3" )

singlet_filter@meta.data$New_Annotation<- Idents(singlet_filter)

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/annotation_aggregated_cluster_filtered_data.pdf'),width = 10,height = 10)
DimPlot(singlet_filter, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(singlet_filter,split.by ="Group", group.by = 'New_Annotation',reduction = 'umap')
DimPlot(singlet_filter,split.by ="Timepoint", group.by = 'New_Annotation',reduction = 'umap')
DimPlot(singlet_filter,split.by ="Treatment", group.by = 'New_Annotation',reduction = 'umap')
dev.off()




#Bubble plot of cell proportion
# Extract metadata
metadata <- singlet_filter @meta.data

# Step 1: Create a summary table with counts
cell_counts <- metadata %>%
  group_by(Annotation, Timepoint, Treatment) %>%
  summarise(num_cells = n(), .groups = "drop")

# Step 2: Calculate total cells for each combination
totals <- metadata %>%
  group_by(Timepoint, Treatment) %>%
  summarise(total_cells = n(), .groups = "drop")

# Step 3: Join the total cell counts to the cell_counts table
cell_percentages <- cell_counts %>%
  left_join(totals, by = c("Timepoint", "Treatment")) %>%
  mutate(percentage = (num_cells / total_cells) * 100)

# Step 4: Convert to a dataframe
cell_percentages_df <- as.data.frame(cell_percentages)

# Load the necessary library
library(dplyr)

# Step 1: Create a new column that combines Timepoint and Treatment
cell_percentages_df$Condition <- paste(cell_percentages_df$Timepoint, cell_percentages_df$Treatment, sep = "_")

# Step 2: Pivot the data - make Annotation the row, Condition the columns, and percentage as the value
df_pivoted <- cell_percentages_df %>%
  pivot_wider(names_from = Condition, values_from = percentage)

# reset the index for easier handling in some cases
df_pivoted_reset <- df_pivoted %>%
  rownames_to_column("RowNumber")

# View the transformed data
head(df_pivoted_reset)


# Load necessary libraries
pdf(paste0("/home/nasim/amgen.abio/amgen_outputs_with_negatives/","CellProportion_BubblePlot.pdf"),width = 10,height = 10)
library(ggplot2)
library(reshape2)

df_pivoted_reset[is.na(df_pivoted_reset)] <- 0
select <- c("Annotation", "Day14_B","Day14_B+T1","Day14_B+T2", "Day14_C","Day21_B","Day21_B+T1","Day21_B+T2", "Day21_C")

df_pivoted_reset <- df_pivoted_reset[,select]

# Melt the data
df_melted <- melt(df_pivoted_reset, id.vars = "Annotation", variable.name = "Condition", value.name = "Proportion")

# Convert Proportion to numeric (if needed)
df_melted$Proportion <- as.numeric(df_melted$Proportion)

df_melted_new <- subset(df_melted, Proportion > 0)


# Custom 18-color palette using HEX values
color_palette <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", 
  "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#8b8b00", "#6a4dff", 
  "#00ff9d", "#f7b3da", "#b42b29", "#55d1c7", "#a11a6b", "#f2b802"
)

# Plot using ggplot2 with x-axis at the top
ggplot(df_melted, aes(x = Condition, y = Annotation, size = Proportion, color = Annotation)) +
  geom_point() +
  scale_size_continuous(range = c(3, 15)) + # Adjust bubble size range
  scale_color_manual(values = color_palette) + # Use the custom color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14), 
    axis.text.y = element_text(size = 14), 
    axis.title.x = element_text(size = 16, face = "bold"), 
    axis.title.y = element_text(size = 16, face = "bold"), 
    plot.title = element_text(size = 20, face = "bold"),
    axis.position = "top", 
    legend.position = "right"
  ) + 
  labs(x = "Condition", y = "Cell Type", title = "Bubble Plot of Cell Proportions") +
  guides(color = "none", size = guide_legend(title = "Proportion")) # Remove color legend and keep size lege
dev.off()



#split the UMAP for each timepoint to treatments, for 18 clusters and aggregated clusters
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_UMAP_Annotation_Split.pdf'),width = 12,height = 10)
Day14 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14")
DimPlot(Day14,split.by ="Treatment", group.by = 'Annotation',reduction = 'umap')
dev.off()
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_UMAP_NewAnnotation_Split.pdf'),width = 12,height = 10)
DimPlot(Day14,split.by ="Treatment", group.by = 'New_Annotation',reduction = 'umap')
dev.off()

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_UMAP_Annotation_Split.pdf'),width = 12,height = 10)
Day21 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21")
DimPlot(Day21, split.by="Treatment", group.by = "Annotation", reduction="umap")
dev.off()
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_UMAP_NewAnnotation_Split.pdf'),width = 12,height = 10)
DimPlot(Day21, split.by = "Treatment", group.by = "New_Annotation", reduction = "umap")
dev.off()

#split the UMAP for each timepoint to treatments, for 18 clusters and aggregated clusters
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_UMAP_B_Annotation.pdf'),width = 12,height = 10)
Day14_B <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="B")
DimPlot(Day14_B,split.by ="Treatment", group.by = 'Annotation',reduction = 'umap')
dev.off()

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_UMAP_B_NewAnnotation.pdf'),width = 12,height = 10)
#Day14_B <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="B")
DimPlot(Day14_B,split.by ="Treatment", group.by = 'New_Annotation',reduction = 'umap')
dev.off()


pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_UMAP_C_Annotation.pdf'),width = 12,height = 10)
Day14_C <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="C")
DimPlot(Day14_C,split.by ="Treatment", group.by = 'Annotation',reduction = 'umap')
dev.off()


pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_UMAP_C_NewAnnotation.pdf'),width = 12,height = 10)
#Day14_C <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="C")
DimPlot(Day14_C,split.by ="Treatment", group.by = 'New_Annotation',reduction = 'umap')
dev.off()

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_UMAP_BT1_Annotation.pdf'),width = 12,height = 10)
Day14_BT1 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="B+T1")
DimPlot(Day14_BT1,split.by ="Treatment", group.by = 'Annotation',reduction = 'umap')
dev.off()


pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_UMAP_BT1_NewAnnotation.pdf'),width = 12,height = 10)
#Day14_BT1 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="B+T1")
DimPlot(Day14_BT1,split.by ="Treatment", group.by = 'New_Annotation',reduction = 'umap')
dev.off()


pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_UMAP_BT2_Annotation.pdf'),width = 12,height = 10)
Day14_BT2 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="B+T2")
DimPlot(Day14_BT2,split.by ="Treatment", group.by = 'Annotation',reduction = 'umap')
dev.off()


pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_UMAP_BT2_NewAnnotation.pdf'),width = 12,height = 10)
#Day14_BT2 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="B+T2")
DimPlot(Day14_BT2,split.by ="Treatment", group.by = 'New_Annotation',reduction = 'umap')
dev.off()


####Day21
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_UMAP_B_Annotation.pdf'),width = 12,height = 10)
Day21_B <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="B")
DimPlot(Day21_B,split.by ="Treatment", group.by = 'Annotation',reduction = 'umap')
dev.off()

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_UMAP_B_NewAnnotation.pdf'),width = 12,height = 10)
#Day21_B <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="B")
DimPlot(Day21_B,split.by ="Treatment", group.by = 'New_Annotation',reduction = 'umap')
dev.off()


pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_UMAP_C_Annotation.pdf'),width = 12,height = 10)
Day21_C <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="C")
DimPlot(Day21_C,split.by ="Treatment", group.by = 'Annotation',reduction = 'umap')
dev.off()

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_UMAP_C_NewAnnotation.pdf'),width = 12,height = 10)
#Day21_C <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="C")
DimPlot(Day21_C,split.by ="Treatment", group.by = 'New_Annotation',reduction = 'umap')
dev.off()




pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_UMAP_BT1_Annotation.pdf'),width = 12,height = 10)
Day21_BT1 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="B+T1")
DimPlot(Day21_BT1,split.by ="Treatment", group.by = 'Annotation',reduction = 'umap')
dev.off()

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_UMAP_BT1_NewAnnotation.pdf'),width = 12,height = 10)
#Day21_BT1 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="B+T1")
DimPlot(Day21_BT1,split.by ="Treatment", group.by = 'New_Annotation',reduction = 'umap')
dev.off()




pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_UMAP_BT2_Annotation.pdf'),width = 12,height = 10)
Day21_BT2 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="B+T2")
DimPlot(Day21_BT2,split.by ="Treatment", group.by = 'Annotation',reduction = 'umap')
dev.off()

pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_UMAP_BT2_NewAnnotation.pdf'),width = 12,height = 10)
#Day21_BT2 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="B+T2")
DimPlot(Day21_BT2,split.by ="Treatment", group.by = 'New_Annotation',reduction = 'umap')
dev.off()



#split the UMAP for each timepoint or treatments for aggregated clusters for feature plot and using different colors
#pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_B_NewAnnotation_FeaturePlot.pdf'),width = 12,height = 10)
Day14_B <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="B")
P1=FeaturePlot(Day14_B, features = c("Cd8a", "Cd4", "Foxp3", "Mki67", "Ifit1", "Pdcd1", "Lag3", "Havcr2", "Tcf7", "Sell", "Il2ra", "Cx3cr1", "Gzma", "Prf1", "Ifng"),ncol = 5) & #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(P1)
#dev.off()


#Day14_B
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_B_FeaturePlot.pdf'),width = 12,height = 10)
Day14_B <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="B")
#custom_palette <- c("yellow", "orange", "purple"))
custom_palette <- colorRampPalette(c("yellow", "orange", "purple"))(100)


P1 <- FeaturePlot(Day14_B, features = c("Cd8a", "Cd4", "Foxp3", "Mki67", "Ifit1", "Pdcd1", "Lag3", "Havcr2", "Tcf7", "Sell","Il2ra", "Cx3cr1", "Gzma", "Prf1", "Ifng", "Tox", "Batf3"), ncol = 5) & scale_colour_gradientn(colours = custom_palette)

# P2 <- FeaturePlot(Day14_B, features = c("Cd8", "Tox", "Batf3", "Il-2", "4-1bb"), ncol = 2) & scale_colour_gradientn(colours = custom_palette)
print(P1)
#print(P2)
dev.off()



#Day14_C
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_C_FeaturePlot.pdf'),width = 12,height = 10)
Day14_C <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="C")
#custom_palette <- c("yellow", "orange", "purple"))
custom_palette <- colorRampPalette(c("yellow", "orange", "purple"))(100)


P1 <- FeaturePlot(Day14_C, features = c("Cd8a", "Cd4", "Foxp3", "Mki67", "Ifit1", "Pdcd1", "Lag3", "Havcr2", "Tcf7", "Sell","Il2ra", "Cx3cr1", "Gzma", "Prf1", "Ifng", "Tox", "Batf3"), ncol = 5) & scale_colour_gradientn(colours = custom_palette)

# P2 <- FeaturePlot(Day14_C, features = c("Cd8", "Tox", "Batf3", "Il-2", "4-1bb"), ncol = 2) & scale_colour_gradientn(colours = custom_palette)
print(P1)
#print(P2)
dev.off()





#Day14_BT1
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_BT1_FeaturePlot.pdf'),width = 12,height = 10)
Day14_BT1 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="B+T1")
#custom_palette <- c("yellow", "orange", "purple"))
custom_palette <- colorRampPalette(c("yellow", "orange", "purple"))(100)


P1 <- FeaturePlot(Day14_BT1, features = c("Cd8a", "Cd4", "Foxp3", "Mki67", "Ifit1", "Pdcd1", "Lag3", "Havcr2", "Tcf7", "Sell","Il2ra", "Cx3cr1", "Gzma", "Prf1", "Ifng", "Tox", "Batf3"), ncol = 5) & scale_colour_gradientn(colours = custom_palette)

# P2 <- FeaturePlot(Day14_BT1, features = c("Cd8", "Tox", "Batf3", "Il-2", "4-1bb"), ncol = 2) & scale_colour_gradientn(colours = custom_palette)
print(P1)
#print(P2)
dev.off()



#Day14_BT2
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day14_BT2_FeaturePlot.pdf'),width = 12,height = 10)
Day14_BT2 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day14" & Treatment=="B+T2")
#custom_palette <- c("yellow", "orange", "purple"))
custom_palette <- colorRampPalette(c("yellow", "orange", "purple"))(100)


P1 <- FeaturePlot(Day14_BT2, features = c("Cd8a", "Cd4", "Foxp3", "Mki67", "Ifit1", "Pdcd1", "Lag3", "Havcr2", "Tcf7", "Sell","Il2ra", "Cx3cr1", "Gzma", "Prf1", "Ifng", "Tox", "Batf3"), ncol = 5) & scale_colour_gradientn(colours = custom_palette)

# P2 <- FeaturePlot(Day14_BT2, features = c("Cd8", "Tox", "Batf3", "Il-2", "4-1bb"), ncol = 2) & scale_colour_gradientn(colours = custom_palette)
print(P1)
#print(P2)
dev.off()



#Day21_B
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_B_FeaturePlot.pdf'),width = 12,height = 10)
Day21_B <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="B")
#custom_palette <- c("yellow", "orange", "purple"))
custom_palette <- colorRampPalette(c("yellow", "orange", "purple"))(100)


P1 <- FeaturePlot(Day21_B, features = c("Cd8a", "Cd4", "Foxp3", "Mki67", "Ifit1", "Pdcd1", "Lag3", "Havcr2", "Tcf7", "Sell","Il2ra", "Cx3cr1", "Gzma", "Prf1", "Ifng", "Tox", "Batf3"), ncol = 5) & scale_colour_gradientn(colours = custom_palette)

# P2 <- FeaturePlot(Day21_B, features = c("Cd8", "Tox", "Batf3", "Il-2", "4-1bb"), ncol = 2) & scale_colour_gradientn(colours = custom_palette)
print(P1)
#print(P2)
dev.off()



#Day21_C
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_C_FeaturePlot.pdf'),width = 12,height = 10)
Day21_C <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="C")
#custom_palette <- c("yellow", "orange", "purple"))
custom_palette <- colorRampPalette(c("yellow", "orange", "purple"))(100)


P1 <- FeaturePlot(Day21_C, features = c("Cd8a", "Cd4", "Foxp3", "Mki67", "Ifit1", "Pdcd1", "Lag3", "Havcr2", "Tcf7", "Sell","Il2ra", "Cx3cr1", "Gzma", "Prf1", "Ifng", "Tox", "Batf3"), ncol = 5) & scale_colour_gradientn(colours = custom_palette)
# 
# P2 <- FeaturePlot(Day21_C, features = c("Cd8", "Tox", "Batf3", "Il-2", "4-1bb"), ncol = 2) & scale_colour_gradientn(colours = custom_palette)
print(P1)
#print(P2)
dev.off()



#Day21_BT1
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_BT1_FeaturePlot.pdf'),width = 12,height = 10)
Day21_BT1 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="B+T1")
#custom_palette <- c("yellow", "orange", "purple"))
custom_palette <- colorRampPalette(c("yellow", "orange", "purple"))(100)


P1 <- FeaturePlot(Day21_BT1, features = c("Cd8a", "Cd4", "Foxp3", "Mki67", "Ifit1", "Pdcd1", "Lag3", "Havcr2", "Tcf7", "Sell","Il2ra", "Cx3cr1", "Gzma", "Prf1", "Ifng", "Tox", "Batf3"), ncol = 5) & scale_colour_gradientn(colours = custom_palette)

# P2 <- FeaturePlot(Day21_BT1, features = c("Cd8", "Tox", "Batf3", "Il-2", "4-1bb"), ncol = 2) & scale_colour_gradientn(colours = custom_palette)
print(P1)
#print(P2)
dev.off()




#Day21_BT2
pdf(paste0('/home/nasim/amgen.abio/amgen_outputs_with_negatives/Day21_BT2_FeaturePlot.pdf'),width = 12,height = 10)
Day21_BT2 <- subset(`3-Data_Annotation_filtered_data`, subset=Timepoint=="Day21" & Treatment=="B+T2")
#custom_palette <- c("yellow", "orange", "purple"))
custom_palette <- colorRampPalette(c("yellow", "orange", "purple"))(100)


P1 <- FeaturePlot(Day21_BT2, features = c("Cd8a", "Cd4", "Foxp3", "Mki67", "Ifit1", "Pdcd1", "Lag3", "Havcr2", "Tcf7", "Sell","Il2ra", "Cx3cr1", "Gzma", "Prf1", "Ifng", "Tox", "Batf3"), ncol = 5) & scale_colour_gradientn(colours = custom_palette)

# P2 <- FeaturePlot(Day21_BT2, features = c("Cd8", "Tox", "Batf3", "Il-2", "4-1bb"), ncol = 2) & scale_colour_gradientn(colours = custom_palette)
print(P1)
#print(P2)
dev.off()


#save
saveRDS(singlet_filter, '/home/nasim/amgen.abio/amgen_outputs_with_negatives/3-Data_Annotation_filtered_data.rds')

