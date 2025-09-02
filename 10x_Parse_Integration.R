# ------------------------------------------------------------------------------
#title: " Parse + 10x Single-Cell RNA-seq Integration — Same Tissue Type"
#Author: Nasim Rahmatpour
#Date: "2024-09-04"
# ---


# Project Overview

#Two single-cell datasets were generated from the same tissue using different technologies:
# - **Parse Biosciences (split-pool barcoding)**
# - **10x Genomics (droplet-based sequencing)**
  
#To enable integrated analysis, I processed and filtered each dataset independently, identified shared genes, and then performed integration to eliminate platform-driven batch effects. Unified clustering and marker discovery were carried out to characterize shared cell populations across platforms.

---
  
###Task: Cross-Platform Integration and Cell Type Characterization  
library(Seurat)
suppressMessages(require(cowplot))
suppressMessages(require(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(ggplot2))
library(ggplot2)
library(RColorBrewer)
library(angrycell)

#create seurat for Parse data
mat_path <- "/home/nasim/genentech.abio/Parse_data/alignment_P7_P8_P9/output_combined/all-sample/DGE_filtered"
mat <- ReadParseBio(mat_path)
table(rownames(mat) =="")
rownames(mat)[rownames(mat)==""] <- "unknown"
cell_meta_P789 <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)
seurat_P789 <- CreateSeuratObject(mat,names.field = 0, meta.data=cell_meta_P789)
seurat_P789$sample[seurat_P789@meta.data$sample=="NGS095_P7"] <- "P7"
seurat_P789$sample[seurat_P789@meta.data$sample=="NGS095_P8"] <- "P8"
seurat_P789$sample[seurat_P789@meta.data$sample=="NGS095_P9"] <- "P9"
seurat_P789$Data <- ""
seurat_P789@meta.data$Data <- "Parse"


mat_path <- "/home/nasim/genentech.abio/Parse_data/alignment_P19_P20_P21/output_combined/all-sample/DGE_filtered"
mat1 <- ReadParseBio(mat_path)
table(rownames(mat1) =="")
rownames(mat1)[rownames(mat1)==""] <- "unknown"
cell_meta_P192021 <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)
seurat_P192021 <- CreateSeuratObject(mat1,names.field = 0, meta.data=cell_meta_P192021)
seurat_P192021$sample[seurat_P192021@meta.data$sample=="NGS095_19"] <- "P19"
seurat_P192021$sample[seurat_P192021@meta.data$sample=="NGS095_20"] <- "P20"
seurat_P192021$sample[seurat_P192021@meta.data$sample=="NGS095_21"] <- "P21"
seurat_P192021$Data <- ""
seurat_P192021@meta.data$Data <- "Parse"

#cerate seurat for 10x data
ids <- c("NGS095_P123", "NGS095_P456")

for (file in ids){
  data <- Read10X(data.dir = paste0("/home/nasim/genentech.abio/10X_data/", file))
  seurat_obj <- CreateSeuratObject(counts =data,project = file)
  assign(file, seurat_obj)
}

#add meta
NGS095_P123<- AddMetaData(object=NGS095_P123, metadata = P123_Cells_Per_Tag)
NGS095_P123$sample <- ""
NGS095_P123$sample[NGS095_P123@meta.data$`Cells Per Tag`=="BC001"] <- "P1"
NGS095_P123$sample[NGS095_P123@meta.data$`Cells Per Tag`=="BC002"] <- "P2"
NGS095_P123$sample[NGS095_P123@meta.data$`Cells Per Tag`=="BC003"] <- "P3"
NGS095_P123$Data <- ""
NGS095_P123@meta.data$Data <- "10x"



NGS095_P456<- AddMetaData(object=NGS095_P456, metadata = P456_Cells_Per_Tag)
NGS095_P456$sample <- ""
NGS095_P456$sample[NGS095_P456@meta.data$`Cells Per Tag`=="BC001"] <- "P4"
NGS095_P456$sample[NGS095_P456@meta.data$`Cells Per Tag`=="BC002"] <- "P5"
NGS095_P456$sample[NGS095_P456@meta.data$`Cells Per Tag`=="BC003"] <- "P6"
NGS095_P456$Data <- ""
NGS095_P456@meta.data$Data <- "10x"

# seurat_10x <- `0-Data_RawObject`
# seurat_10x$sample <- ""
# seurat_10x@meta.data$sample <- apply(seurat_10x@meta.data, 1, function(x) str_split(x[1], "Data_"))
# seurat_10x@meta.data$sample <- apply(seurat_10x@meta.data,1,function(x) x[[4]][[1]][2])

#merge
combined <- merge(NGS095_P123, c(NGS095_P456,seurat_P789, seurat_P192021), add.cell.ids = c("NGS095_P123", "NGS095_P456","seurat_P789", "seurat_P192021"), project = "Genentech.10x.Parse")
combined[["percent.mt"]] <- PercentageFeatureSet(object = combined, pattern = "^MT-")

# filtering
combined.filter <- subset(combined, subset = nFeature_RNA >= 200 & percent.mt < 10)

#Preprocess
# Step 1: Set the default assay to RNA
DefaultAssay(combined.filter) <- "RNA"

# Step 2: Split the object by sample
data.list <- SplitObject(combined.filter, split.by = "sample")

# Step 3: Intersect the gene sets across all datasets to ensure compatibility
common.genes <- Reduce(intersect, lapply(data.list, rownames))
data.list <- lapply(data.list, function(x) {
  x <- subset(x, features = common.genes)  # Subset to common genes
  return(x)
})

# Step 4: Normalize and identify variable features
data.list <- lapply(data.list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE, selection.method = "vst", nfeatures = 3000)
  return(x)
})


#Integrate
features <- SelectIntegrationFeatures(object.list = data.list)
anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)

dir.create("/home/nasim/genentech.abio/integrated_outputs")
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)

pdf('/home/nasim/genentech.abio/integrated_outputs/PCA_plots.pdf')
DimPlot(integrated,group.by = 'sample',label=T,repel=T)
DimHeatmap(integrated, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#Elbow plot
pdf('/home/nasim/genentech.abio/integrated_outputs/Elbow_plot.pdf')
P1 <- ElbowPlot(integrated)
print(P1)
dev.off()

#clustering
integrated <- FindNeighbors(integrated, dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.5)
integrated <- RunTSNE(integrated, dims = 1:20)
integrated <- RunUMAP(integrated, dims = 1:20)

#umap, tsne
pdf('/home/nasim/genentech.abio/integrated_outputs/umap_tsne_plots.pdf', width = 12,height = 10)
DimPlot(integrated,group.by = 'seurat_clusters',reduction = 'umap',label=T,repel=T)
DimPlot(integrated,group.by = 'seurat_clusters',reduction = 'tsne',label=T,repel=T)
DimPlot(integrated,group.by = 'sample',reduction = 'umap')
DimPlot(integrated,group.by = 'sample',reduction = 'tsne')
DimPlot(integrated,group.by = 'Data',reduction = 'umap')
DimPlot(integrated,group.by = 'Data',reduction = 'tsne')
dev.off()

#Find gene markers for each cluster
#if using RNA assay, normalize and join the layers
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)#then Scale too
integrated <- JoinLayers(integrated)

#if using integrated assay, join layers does not work or is not needed but makes NaN
DefaultAssay(integrated) <- "integrated"

markers<- FindAllMarkers(integrated , only.pos = TRUE)
write.xlsx(markers,file = paste0("/home/nasim/genentech.abio/integrated_outputs/RNA_markers.xlsx"),row.names = TRUE)

top20.markers_RNA <- markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top20 <- c()
for(c in unique(top20.markers_RNA$cluster)){
  tmp <- top20.markers_RNA[top20.markers_RNA$cluster==c,]$gene
  top20 <- cbind(top20,tmp)
}
colnames(top20) <- paste0('C_',unique(top20.markers_RNA$cluster))
write.csv(top20,"/home/nasim/genentech.abio/integrated_outputs/top20.RNA_markers.csv")

#save
saveRDS(integrated, '/home/nasim/genentech.abio/integrated_outputs/0-Data_Integrated.rds')

#calculate the percentage of cells in each cluster
library(openxlsx)
meta.df <- table(integrated@meta.data$seurat_clusters,integrated@meta.data$sample)
meta.df
meta.totalper.sample <- apply(meta.df,2,sum)
meta.totalper.sample
meta.ratio <- sweep(meta.df, 2, meta.totalper.sample, `/`)
meta.ratio
meta.pct <- meta.ratio*100
meta.pct 
meta.pctsig <- round(meta.pct, digits = 2)
meta.pctsig

wb <- createWorkbook()

addWorksheet(wb, "Sheet1")
addWorksheet(wb, "Sheet2")

writeData(wb, sheet = "Sheet1", x = meta.df)
writeData(wb, sheet = "Sheet2", x = meta.pctsig)

saveWorkbook(wb, file = "/home/nasim/genentech.abio/integrated_outputs/cluster_count_percentage_per_Sample.xlsx", overwrite = TRUE)

#percentage of cells
pdf(paste0('/home/nasim/genentech.abio/integrated_outputs/cluster_percentage_per_Sample_per_Data.pdf'),width=15, height=10)
library(angrycell)
integrated@meta.data$seurat_clusters <- Idents(integrated)
meta <- integrated@meta.data
angrycell::plot_fraction(meta, x = 'sample', fill.bar = 'seurat_clusters')


library(angrycell)
integrated@meta.data$seurat_clusters <- Idents(integrated)
meta <- integrated@meta.data
angrycell::plot_fraction(meta, x = 'Data', fill.bar = 'seurat_clusters')
dev.off()

#save the common genes
write.csv(as.data.frame(common.genes),"/home/nasim/genentech.abio/common_genes.csv")

#feature and dot plot for annotation
#pdf(paste0("/home/nasim/genentech.abio/","gene_markers_featureplot.pdf"),width = 10,height = 10)
pdf(paste0("/home/nasim/genentech.abio/","gene_markers_featureplot_test.pdf"),width = 10,height = 10)
#DefaultAssay(`0-Data_Integrated`) <- "RNA"
DefaultAssay(`0-Data_Integrated`) <- "integrated"
#neural_progenitor_cell
P1=FeaturePlot(`0-Data_Integrated` , features = c("SOX2", "VIM", "NES", "HOPX", "BCAN", "TNC", "FOXG1", "EMX1", "GABRA2"),  ncol = 3)

P12=FeaturePlot(`0-Data_Integrated` , features = c("GABRB1", "DLX2", "TCF7L2", "SIX3", "OTX2", "PAX7", "VSX2", "CYP26A1", "HOXA2"),  ncol = 3)


P13=FeaturePlot(`0-Data_Integrated` , features = c( "HOXB2", "UNC5C", "GPC6", "HOXD3", "HOXA3"),  ncol = 3)


#neuron
P2=FeaturePlot(`0-Data_Integrated` , features = c("STMN2", "DCX", "SLC17A7", "SLC17A6", "FOXG1","EMX1", "BCL11B","SATB2","GABRA2"),  ncol = 3)

P21=FeaturePlot(`0-Data_Integrated` , features = c("GABRB1",  "TCF7L2", "PITX2", "BARHL2", "TFAP2D","TFAP2A", "HOXA2", "HOXB2", "HOXD3"),  ncol = 3)

P22=FeaturePlot(`0-Data_Integrated` , features = c( "HOXA3","UNCX", "INSM1", "UNC5C", "ROBO1", "NEGR1"),  ncol = 3)


#inhibitory_neuron
P3=FeaturePlot(`0-Data_Integrated` , features = c("GAD1", "GAD2", "SLC32A1", "FOXG1", "DLX2", "DLX5", "NKX2-1",  "ISL1", "NR2F1"),  ncol = 3)


P31=FeaturePlot(`0-Data_Integrated` , features = c("LHX5", "LHX1", "TCF7L2", "DLX1", "OTX2", "SKOR2", "HOXA2", "HOXB2", "CA8"),  ncol = 3)

P32=FeaturePlot(`0-Data_Integrated` , features = c("LAMP5", "TFAP2A", "UNC5C", "HOXD3", "HOXA3", "LAMP5", "ROBO1", "NEGR1"),  ncol = 3)


#choroid_plexus_epithelium
P4=FeaturePlot(`0-Data_Integrated` , features = c("TTR", "GFAP", "AQP4", "OLIG1", "MBP"),  ncol = 3)

#microglia and other cells
P5=FeaturePlot(`0-Data_Integrated` , features = c("AIF1", "CLDN5", "DCN", "SOX10", "PRPH"),  ncol = 3)

print(P1)
print(P12)
print(P13)
print(P2)
print(P21)
print(P22)
print(P3)
print(P31)
print(P32)
print(P4)
print(P5)
dev.off()

#pdf(paste0("/home/nasim/genentech.abio/","gene_markers_dotplot.pdf"),width = 10,height = 10)
pdf(paste0("/home/nasim/genentech.abio/","gene_markers_dotplot_test.pdf"),width = 10,height = 10)
P6=angrycell::DotPlot2(`0-Data_Integrated`, features = c('SOX2', 'VIM', 'NES', 'FOXG1', 'EMX1', 'HOPX', 'GABRA2', 'GABRB1', 'DLX2', 'TCF7L2', 'SIX3', 'OTX2', 'PAX7', 'VSX2', 'CYP26A1', 'HOXA2', 'HOXB2', 'UNC5C', 'BCAN', 'GPC6', 'HOXD3', 'HOXA3', 'STMN2', 'DCX', 'SLC17A7', 'SLC17A6', 'BCL11B', 'SATB2', 'PITX2', 'BARHL2', 'TFAP2D', 'TFAP2A', 'UNCX', 'INSM1', 'ROBO1', 'NEGR1', 'GAD1', 'GAD2', 'SLC32A1', 'DLX5', 'NKX2-1', 'ISL1', 'NR2F1', 'LHX5', 'LHX1', 'CA8', 'LAMP5', 'SKOR2', 'AQP4', 'GFAP', 'TTR', 'OLIG1', 'MBP', 'AIF1', 'CLDN5', 'DCN', 'SOX10', 'PRPH' ), dot.scale = 8) +
  coord_flip() +
  theme(text = element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=8),
        axis.title.y = element_text(size =10),
        axis.title.x = element_text(size=10))

print(P6)
dev.off()


#add anno
Idents(`0-Data_Integrated`) <-`0-Data_Integrated`@meta.data$seurat_clusters
`0-Data_Integrated` <- RenameIdents(`0-Data_Integrated`, `0` = "Excitatory neuron", `1` = "Excitatory neuron", `2` = "Neural progenitor cells", `3` = "Excitatory neuron", `4` = "Excitatory neuron", `5` ="Excitatory neuron", `6` = "Neural progenitor cells", `7` = "Excitatory neuron",`8` = "Excitatory neuron", `9` = "Neuron", `10`= "Neural progenitor cells", `11` = "Excitatory neuron", `12`= "Inhibitory neuron", `13`="Glial Cells" , `14`= "Neuron", `15`= "Neural progenitor cells" , `16`= "Inhibitory neuron", `17`= "Inhibitory neuron", `18`= "Neural progenitor cells" )

`0-Data_Integrated`@meta.data$Annotation<- Idents(`0-Data_Integrated`)

pdf(paste0("/home/nasim/genentech.abio/annotation_umap.pdf"),width = 10,height = 10)
DimPlot(`0-Data_Integrated`, reduction = "umap", label = TRUE, repel = TRUE)
# DimPlot(`0-Data_Integrated`, split.by ="sample", group.by ="Annotation" ,reduction = 'umap', label = TRUE, repel = TRUE)
DimPlot(`0-Data_Integrated`, split.by ="Data", group.by ="Annotation" ,reduction = 'umap', label = TRUE, repel = TRUE)
dev.off()

#change the Data column

#P123=10X Flex fix-dissociate
#P456=10X Flex dissociate-fix 
#P789 and P192021=Parse

`0-Data_Integrated`$Data[`0-Data_Integrated`@meta.data$orig.ident=="NGS095_P123"] <- "10X Flex fix-dissociate"
`0-Data_Integrated`$Data[`0-Data_Integrated`@meta.data$orig.ident=="NGS095_P456"] <- "10X Flex dissociate-fix"
`0-Data_Integrated`$Data[`0-Data_Integrated`@meta.data$orig.ident=="SeuratProject"] <- "Parse"

pdf(paste0('/home/nasim/genentech.abio/clusters_qualitycheck.pdf'),width = 10,height = 10)
VlnPlot(`0-Data_Integrated`, features = c('percent.mt'),group.by='seurat_clusters', pt.size = 0)+theme(legend.position = 'none')
VlnPlot(`0-Data_Integrated`, features = c('nFeature_RNA'),group.by='seurat_clusters', pt.size = 0)+theme(legend.position = 'none')
VlnPlot(`0-Data_Integrated`, features = c('nCount_RNA'),group.by='seurat_clusters', pt.size = 0)+theme(legend.position = 'none')
dev.off()

#feature plot on top genes of cluster 2, 6, 10, 15, 18
pdf(paste0("/home/nasim/genentech.abio/","top20_markers_clusters_2_6_10_15_18_featureplot.pdf"),width = 10,height = 10)
#cluster 2
P1=FeaturePlot(`0-Data_Integrated` , features = c("TNC", "INPP5D", "ADRA1B", "HERC5", "GDF15", "EFCC1", "TRABD2B", "HAS2", "SLC17A8", "LRP2"),  ncol = 3)
P2=FeaturePlot(`0-Data_Integrated` , features = c("THBS1", "HMCN1", "SGCG", "ATP6V0D2", "COL24A1", "CD96", "ADAM12", "ATP1A4", "HKDC1", "EGF"),  ncol = 3)
print(P1)
print(P2)

#cluster6
P3=FeaturePlot(`0-Data_Integrated` , features = c("SEMA5A", "IQGAP2", "VEPH1", "IFI44L", "CRB2", "S1PR1", "LIPG", "KCNJ10", "MLC1", "FZD8"),  ncol = 3)
P4=FeaturePlot(`0-Data_Integrated` , features = c("LRRC17", "TKTL1", "OAF", "SFRP1", "ITGA2", "RHOJ", "SLCO1C1", "GPX3", "OTOF", "PTH2R"),  ncol = 3)
print(P3)
print(P4)


#cluster10
P5=FeaturePlot(`0-Data_Integrated` , features = c("RRM2", "SPC25", "MYBL2", "MCM10", "RAD51AP1", "RAD51", "RAD54L", "DTL", "E2F8", "XRCC2"),  ncol = 3)
P6=FeaturePlot(`0-Data_Integrated` , features = c("CLSPN", "POLQ", "CDCA5", "ZWINT", "EXO1", "CDC45", "FAM111B", "PCLAF", "TYMS", "ESCO2"),  ncol = 3)
print(P5)
print(P6)


#cluster15
P7=FeaturePlot(`0-Data_Integrated` , features = c("KIF2C", "CDCA8", "DLGAP5", "KIF20A", "TACC3", "CDCA2", "BUB1", "CENPF", "CDC20", "HJURP"),  ncol = 3)
P8=FeaturePlot(`0-Data_Integrated` , features = c("TOP2A", "CCNB2", "TROAP", "PIF1", "CEP55", "FAM83D", "UBE2C", "CCNB1", "AURKA", "NEK2"),  ncol = 3)
print(P7)
print(P8)

#cluster18
P9=FeaturePlot(`0-Data_Integrated` , features = c("TNC", "SALL3", "F3", "VEPH1", "LIX1", "FAM107A", "RAB31", "HEPACAM", "MOXD1", "PMP2"),  ncol = 3)
P10=FeaturePlot(`0-Data_Integrated` , features = c("CXCL14", "NAALAD2", "MMD2", "FBLN5", "EFEMP1", "ACSBG1", "SLCO1C1", "MAOB", "RPE65", "ARAP2"),  ncol = 3)
print(P9)
print(P10)

#Proliferation markers
P11=FeaturePlot(`0-Data_Integrated` , features = c("MKI67","TOP2A" ),  ncol = 2)
print(P11)
dev.off()



#add new anno
Idents(`0-Data_Integrated`) <-`0-Data_Integrated`@meta.data$seurat_clusters
`0-Data_Integrated` <- RenameIdents(`0-Data_Integrated`, `0` = "Dorsal Telen Excitatory neuron (Deep layer cortical)", `1` = "Dorsal Telen Excitatory neuron (Deep layer cortical)", `2` = "Rhombencephalic Inhibitory neuron", `3` = "Dorsal Telen Excitatory neuron (Deep layer cortical)", `4` = "Dorsal Telen Excitatory neuron (Deep layer cortical)", `5` ="Dorsal Telen Excitatory neuron (Deep layer cortical)", `6` = "Rhombencephalic Inhibitory neuron", `7` = "Dorsal Telen Excitatory neuron (Deep layer cortical)",`8` = "Dorsal Telen Excitatory neuron (Upper layer cortical)", `9` = "Rhombencephalic Excitatory Neuron", `10`= "Proliferating Telencephalic NPC", `11` = "Dorsal Telen Excitatory neuron (Deep layer cortical)", `12`= "Telencephalic Inhibitory neuron", `13`="Diencephalic NPC" , `14`= "Low Quality Cluster", `15`= "Proliferating Telencephalic NPC" , `16`= "Telencephalic Inhibitory neuron", `17`= "Non-Telencephalic Inhibitory neuron", `18`= "Ventral Telen NPC, future Inhibitory neuron" )

`0-Data_Integrated`@meta.data$Annotation<- Idents(`0-Data_Integrated`)

pdf(paste0("/home/nasim/genentech.abio/annotation_new_umap.pdf"),width = 10,height = 10)
DimPlot(`0-Data_Integrated`, reduction = "umap", label = TRUE, repel = TRUE)
# DimPlot(`0-Data_Integrated`, split.by ="sample", group.by ="Annotation" ,reduction = 'umap', label = TRUE, repel = TRUE)
DimPlot(`0-Data_Integrated`, split.by ="Data", group.by ="Annotation" ,reduction = 'umap', label = FALSE, repel = TRUE)
dev.off()



#change the Data column

#P123=10X Flex fix-dissociate
#P456=10X Flex dissociate-fix 
#P789 and P192021=Parse

`0-Data_Integrated`$Data[`0-Data_Integrated`@meta.data$orig.ident=="NGS095_P123"] <- "10X Flex fix-dissociate"
`0-Data_Integrated`$Data[`0-Data_Integrated`@meta.data$orig.ident=="NGS095_P456"] <- "10X Flex dissociate-fix"

`0-Data_Integrated`$Data[`0-Data_Integrated`@meta.data$sample=="P7" | `0-Data_Integrated`@meta.data$sample=="P8" | `0-Data_Integrated`@meta.data$sample=="P9"] <- "Parse_R1"

`0-Data_Integrated`$Data[`0-Data_Integrated`@meta.data$sample=="P19" | `0-Data_Integrated`@meta.data$sample=="P20" | `0-Data_Integrated`@meta.data$sample=="P21"] <- "Parse_R2"


#save
saveRDS(`0-Data_Integrated`, '/home/nasim/genentech.abio/StandardAnalysis_10x_Parse/integrated_outputs/Rdata/1-Data_Integrated_Annotation.rds')  