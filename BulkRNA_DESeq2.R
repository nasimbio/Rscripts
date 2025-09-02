
# ------------------------------------------------------------------------------
# title: "Bulk RNA-seq Analysis — Brain Region Comparison"
# Author: Nasim Rahmatpour
# Date: "2024-06-24"
# ------------------------------------------------------------------------------


# Project Overview

#The goal of this project is to investigate transcriptional differences within and between the hippocampus (HP) and motor cortex (MC) brain regions by analyzing bulk RNA-seq data. Each brain region is subdivided into groups A and B. Differential expression and Gene Ontology enrichment analyses were conducted to explore region-specific biological variation and experimental effects.

# **Data Preprocessing**
# - Raw FASTQ files were quality checked and aligned to the reference human genome
# - Gene count matrices were generated from aligned read

---
  
###Task: Differential Gene Expression (DGE) and GO Enrichment  
  
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")
library("ggrepel")
library("EnhancedVolcano")
#library("ggpmisc")
library(ggrepel)

t <- as.data.frame(gex_matrix, check.names = FALSE)
write.csv(t, "/home/nasim/uky/gene_expression_matrix.csv")

#dir.create("./HP-Bs.vs.HP-As")
dir.create("./MC-Bs.vs.MC-As")
#dir.create("./MC-As.vs.HP-As")
#dir.create("./MC-Bs.vs.HP-Bs")

# names(gene_expression_matrix)[names(gene_expression_matrix) == "...1"] <- 'genes'
# rownames(gene_expression_matrix) <- gene_expression_matrix$genes
# HP_Bs <- gene_expression_matrix[, c("genes", "1HP","6HP","7HP","8HP","9HP", "10HP")]
# MC_Bs <- gene_expression_matrix[, c("genes","1MC","6MC","7MC", "8MC", "9MC", "10MC")]
# #rownames(HP_As) <- rownames(gene_expression_matrix)
# #rownames(HP_Bs) <- rownames(gene_expression_matrix)
# 
# sampleNames <- c("1HP","6HP","7HP","8HP","9HP", "10HP", "1MC","6MC","7MC", "8MC", "9MC", "10MC" )
# sampleCondition <- c("HP-B", "HP-B", "HP-B", "HP-B","HP-B", "HP-B", "MC-B", "MC-B", "MC-B", "MC-B", "MC-B", "MC-B")
# sampleFiles<- c("HP", "HP", "HP","HP","HP","HP", "MC", "MC", "MC", "MC", "MC", "MC")
# SampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
# rownames(SampleTable) <- SampleTable$sampleName
# #countData_new <-cbind(HP_As, HP_Bs)
# countData_new <- merge(HP_Bs, MC_Bs, by= "genes")
# 
# rownames(countData_new) <- countData_new$genes
# 
# countData_new <- countData_new[ , c("1HP","6HP","7HP","8HP","9HP", "10HP", "1MC","6MC","7MC", "8MC", "9MC", "10MC")]
# rownames(SampleTable)
# colnames(countData_new)
# all(rownames(SampleTable) %in% colnames(countData_new))
# all(rownames(SampleTable) == colnames(countData_new))



###Remove outliers
names(gene_expression_matrix)[names(gene_expression_matrix) == "...1"] <- 'genes'
rownames(gene_expression_matrix) <- gene_expression_matrix$genes
HP_As <- gene_expression_matrix[, c("genes","3HP","4HP", "5HP","11HP")]
HP_Bs <- gene_expression_matrix[, c("genes","6HP","7HP", "8HP", "9HP", "10HP")]
#rownames(HP_As) <- rownames(gene_expression_matrix)
#rownames(HP_Bs) <- rownames(gene_expression_matrix)

sampleNames <- c("3HP","4HP", "5HP","11HP", "6HP","7HP", "8HP", "9HP", "10HP")
sampleCondition <- c("HP-A","HP-A", "HP-A", "HP-A", "HP-B", "HP-B", "HP-B", "HP-B", "HP-B")
sampleFiles<- c("HP","HP","HP","HP", "HP", "HP", "HP", "HP", "HP")
SampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
rownames(SampleTable) <- SampleTable$sampleName
#countData_new <-cbind(HP_As, HP_Bs)
countData_new <- merge(HP_As, HP_Bs, by= "genes")

rownames(countData_new) <- countData_new$genes

countData_new <- countData_new[ , c("3HP","4HP", "5HP","11HP", "6HP","7HP", "8HP", "9HP", "10HP" )]
rownames(SampleTable)
colnames(countData_new)
all(rownames(SampleTable) %in% colnames(countData_new))
all(rownames(SampleTable) == colnames(countData_new))

dds <- DESeqDataSetFromMatrix(countData=countData_new, 
                              colData=SampleTable, 
                              design=~condition)

#dds$condition <- factor(dds$condition, levels = c("MC-A","MC-B"))
#or
dds$condition <- relevel(dds$condition, ref = "HP-A")
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)
#res <- results(dds,contrast = c("HP_As.vs.HP_Bs"))
res <- results(dds, name="condition_HP.B_vs_HP.A")


#res= subset(res, padj<0.1)
#res <- res[order(res$padj),]
res <- lfcShrink(dds, coef="condition_HP.B_vs_HP.A", type="apeglm")
resOrdered <- res[order(res$pvalue),]

#volcano
pdf(paste0("/home/nasim/uky/HP-Bs.vs.HP-As/","HP-Bs.vs.HP-As_no_outlier_Volcano_lower_threshold.pdf"), width=20, height=10)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)
dev.off()

#MA plot
plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")
#plotMA(resLFC, ylim=c(-2,2))

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

mcols(res)$description

write.csv(as.data.frame(resOrdered), 
          file="/home/nasim/uky/HP-Bs.vs.HP-As/HP-Bs.vs.HP-As_no_outliers_results.csv")

vsd <- vst(dds)
head(assay(vsd), 3)

#PCA plot
#pdf(paste0("/home/nasim/uky/HP-Bs.vs.HP-As/", "/HP-Bs.vs.HP-As_PCA_label.pdf"), width=10, height=10)
plotPCA(vsd, intgroup=c("condition"))+
  #geom_label(aes(label = name), vjust=0.5,hjust= 0.5, check_overlap = TRUE,size = 1)+
  geom_text(aes(label=name),vjust=2,check_overlap = TRUE,size = 2)+
  geom_point(size=0.5)+
  theme(legend.position="none")
#dev.off()

#heatmap
pdf(paste0("/home/nasim/uky/HP-Bs.vs.HP-As/", "HP-Bs.vs.HP-As_no_ouliers_Heatmap.pdf"), width=15, height=12)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("sampleName","condition")])
heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE, show_colnames =FALSE,
                    cluster_cols=TRUE, annotation_col=df)

dev.off()

# save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
#   png(filename, width = width, height = height, res = res)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
#}

#

#GO Enrichment starts from here
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

###prepare
#these 2 ones don't have up/down regulated significantly and the 2 rest don't have enough genes for GO enrichment
#markers_new<-read.csv("/home/nasim/uky/HP-Bs.vs.HP-As/HP-Bs.vs.HP-As_results.csv")
#markers_new<-read.csv("/home/nasim/uky/MC-Bs.vs.MC-As/MC-Bs.vs.MC-As_results.csv")
#markers_new<-read.csv("/home/nasim/uky/MC-As.vs.HP-As/MC-As.vs.HP-As_results.csv")
markers_new<-read.csv("/home/nasim/uky/MC-Bs.vs.HP-Bs/MC-Bs.vs.HP-Bs_results.csv")

markers_new$expression <- "not significant"
markers_new$expression[markers_new$pvalue <  0.05 & markers_new$log2FoldChange > 1] <- "upregulated"
markers_new$expression[markers_new$pvalue < 0.05 & markers_new$log2FoldChange < -1] <- "downregulated"
ups <- subset(markers_new, expression=="upregulated")
colnames(ups)[1] <- "gene"
downs <- subset(markers_new, expression=="downregulated")
colnames(downs)[1] <- "gene"

#GOenrich
markerlist_up <- as.character(ups$gene)
eg <- bitr(markerlist_up  , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Hs.eg.db)
genelist_up <- eg$SYMBOL
Markers_GO_enrich_up <- enrichGO(genelist_up , OrgDb = org.Hs.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


markerlist_down <- as.character(downs$gene)
eg1 <- bitr(markerlist_down , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Hs.eg.db)
genelist_down <- eg1$SYMBOL
Markers_GO_enrich_down <- enrichGO(genelist_down , OrgDb = org.Hs.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)

#save
pdf("/home/nasim/uky/MC-Bs.vs.HP-Bs/upregulatedDEGs_GOenrichment.pdf")
enrichplot::dotplot(Markers_GO_enrich_up, split = "ONTOLOGY",showCategory = 10) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("")+
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=6),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))
dev.off()
upregulated_GO_enrich <- Markers_GO_enrich_up@result
write.csv(upregulated_GO_enrich, "/home/nasim/uky/MC-Bs.vs.HP-Bs/upregulatedDEGs_GOenrichment.csv")


pdf("/home/nasim/uky/MC-Bs.vs.HP-Bs/downregulatedDEGs_GOenrichment.pdf")
enrichplot::dotplot(Markers_GO_enrich_down, split = "ONTOLOGY",showCategory = 10) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("")+
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=6),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))
downregulated_GO_enrich <- Markers_GO_enrich_down@result
dev.off()
write.csv(downregulated_GO_enrich, "/home/nasim/uky/MC-Bs.vs.HP-Bs/downregulatedDEGs_GOenrichment.csv")

#KEGGenrich
markerlist_up <- as.character(ups$gene)
eg <- bitr(markerlist_up  , fromType="SYMBOL", toType=c("ENTREZID"), OrgDb = org.Hs.eg.db)
genelist_up <- as.character(eg$ENTREZID)

eg1<- bitr_kegg(genelist_up, fromType = "kegg", toType = "Path", organism = "hsa")
genes_up <-as.character(eg1$kegg)
kegg_enrich_up <- enrichKEGG(gene = genes_up, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.5)
head(kegg_enrich_up)
pdf("/home/nasim/uky/MC-Bs.vs.HP-Bs/upregulatedDEGs_KEGGenrichment.pdf")
enrichplot::dotplot(kegg_enrich_up, showCategory=10)
dev.off()

upregulated_KEGG_enrich <- kegg_enrich_up@result
write.csv(upregulated_KEGG_enrich, "/home/nasim/uky/MC-Bs.vs.HP-Bs/upregulatedDEGs_KEGGenrichment.csv")



markerlist_down <- as.character(downs$gene)
eg3 <- bitr(markerlist_down  , fromType="SYMBOL", toType=c("ENTREZID"), OrgDb = org.Hs.eg.db)
genelist_down <- as.character(eg3$ENTREZID)

eg4<- bitr_kegg(genelist_down, fromType = "kegg", toType = "Path", organism = "hsa")
genes_down <-as.character(eg4$kegg)
kegg_enrich_down <- enrichKEGG(gene = genes_down, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.5)
head(kegg_enrich_down)
pdf("/home/nasim/uky/MC-Bs.vs.HP-Bs/downregulatedDEGs_KEGGenrichment.pdf")
enrichplot::dotplot(kegg_enrich_down, showCategory=10)
dev.off()

downregulated_KEGG_enrich <- kegg_enrich_down@result
write.csv(downregulated_KEGG_enrich, "/home/nasim/uky/MC-Bs.vs.HP-Bs/downregulatedDEGs_KEGGenrichment.csv")
