#This part is for DGE anlysis based on pesudo bulk DGE in each cell line between 2 treatments, making heatmap and volcano plots,and finding the overlap of up/down regulated gene among 3 cell lines
Epithelial <- subset(BI.filter_YFP_AddModuleScore, subset= MajorCellType == "Epithelial" |  MajorCellType == "Proliferating Epithelial") 

Epithelial@meta.data$celltype <- "Epithelial"

#extract the counts and metadata to make singleCellExperiment
counts <- Epithelial@assays$RNA@counts
metadata <- Epithelial@meta.data

#vector is a list of atomic values, factor is a list of vectors, set up metadata as desired for aggregation
#metadata$MajorCellType <- factor(Epithelial@active.ident)
#metadata$MajorCellType[metadata$MajorCellType == "Proliferating Epithelial"] <- "Epithelial"
metadata$celltype <- factor(Epithelial@meta.data$celltype)
metadata$Sample <- factor(Epithelial@meta.data$Sample)
metadata$Treatment <- factor(Epithelial@meta.data$Treatment)

#create SingleCellExperiment
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

#Identify groups for aggregation of counts, (metadata is the colData in sce)
groups <- colData(sce)[, c("celltype", "Sample")]


assays(sce)
dim(counts(sce))
counts(sce)[1:6, 1:6]
dim(colData(sce))
head(colData(sce))

#determine how many cluster and how many samples do we have and what is the length of them
kids <- purrr::set_names(levels(sce$celltype))

kids
# Total number of clusters
nk <- length(kids)
nk

samples <- purrr::set_names(levels(sce$Sample))
# Total number of samples 
n_samples <- length(samples)
n_samples

#generate sample level metadata
#get the number of cells for each sample
table(sce$Sample)
#only get the numbers
n_cells <- as.numeric(table(sce$Sample))
#determine how to reorder the samples or the row of metadata to match the order of sample names in samples (look above)
m <- match(samples, sce$Sample)

#make the metadata, n_cells is a column in ei metadata now, colnames(ei). In ei maybe the rownames and the n_cells column is important
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"celltype")
ei

#aggregate the matrix per sample (per sample per MajorCellType)
groups <- colData(sce)[, c("celltype", "Sample")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum")

class(pb)

dim(pb)

pb[1:6, 1:6]

#split the cell type
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

#split the matrix data by celltype and that the genes are row names and the sample names are the column
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+[:punct:]+[:alnum:]+")))

class(pb)
# Explore the different components of list
str(pb)

#the counts per sample per major cell type can be checked
options(width = 100)
table(sce$celltype, sce$Sample)

#make sample level metadata, (maybe all samples are not present in each cell type)
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}
de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()


samples_list <- map(1:length(kids), get_sample_ids)
get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

gg_df <- data.frame(celltype = de_cluster_ids,
                    Sample = de_samples)
gg_df <- left_join(gg_df, ei[, c("Sample", "Treatment")]) 

metadata <- gg_df %>%
  dplyr::select(celltype, Sample, Treatment) 

metadata

clusters <- levels(factor(metadata$celltype))
clusters

clusters[1]

cluster_metadata <- metadata[which(metadata$celltype == clusters[1]), ]
rownames(cluster_metadata) <- cluster_metadata$Sample
head(cluster_metadata)

counts <- pb[[clusters[1]]]
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
all(rownames(cluster_metadata) == colnames(cluster_counts))

#subset based on the cell line
TC3 <- cluster_counts[which(colnames(cluster_counts) ==c("scGEX_152", "scGEX_153", "scGEX_154", "scGEX_155", "scGEX_156", "scGEX_157", "scGEX_158"))]
cluster_metadata_TC3 <- cluster_metadata[which(rownames(cluster_metadata) ==c("scGEX_152", "scGEX_153", "scGEX_154", "scGEX_155", "scGEX_156", "scGEX_157", "scGEX_158")),]

TC4 <- cluster_counts[which(colnames(cluster_counts) ==c("scGEX_159", "scGEX_160", "scGEX_161", "scGEX_162", "scGEX_163", "scGEX_164", "scGEX_165"))]
cluster_metadata_TC4 <- cluster_metadata[which(rownames(cluster_metadata) ==c("scGEX_159", "scGEX_160", "scGEX_161", "scGEX_162", "scGEX_163", "scGEX_164", "scGEX_165")),]

TC5 <- cluster_counts[which(colnames(cluster_counts) ==c("scGEX_166", "scGEX_167", "scGEX_168", "scGEX_169", "scGEX_170", "scGEX_171", "scGEX_172"))]
cluster_metadata_TC5 <- cluster_metadata[which(rownames(cluster_metadata) ==c("scGEX_166", "scGEX_167", "scGEX_168", "scGEX_169", "scGEX_170", "scGEX_171", "scGEX_172")),]

#make the DEseq object
dds_TC3 <- DESeqDataSetFromMatrix(TC3, 
                                  colData = cluster_metadata_TC3, 
                                  design = ~ Treatment)

dds_TC4 <- DESeqDataSetFromMatrix(TC4, 
                                  colData = cluster_metadata_TC4, 
                                  design = ~ Treatment)

dds_TC5 <- DESeqDataSetFromMatrix(TC5, 
                                  colData = cluster_metadata_TC5, 
                                  design = ~ Treatment)

#plotting based on Treatmant in each cell line
rld_TC3 <- rlog(dds_TC3, blind=TRUE)
DESeq2::plotPCA(rld_TC3, intgroup = "Treatment")

rld_TC4 <- rlog(dds_TC4, blind=TRUE)
DESeq2::plotPCA(rld_TC4, intgroup = "Treatment")

rld_TC5 <- rlog(dds_TC5, blind=TRUE)
DESeq2::plotPCA(rld_TC5, intgroup = "Treatment")

#plot that shows sample relationship for each cell line
rld_mat_TC3 <- assay(rld_TC3)
rld_cor_TC3 <- cor(rld_mat_TC3)
pheatmap(rld_cor_TC3, annotation = cluster_metadata_TC3[ c("celltype", "Treatment"), drop=F])

rld_mat_TC4 <- assay(rld_TC4)
rld_cor_TC4 <- cor(rld_mat_TC4)
pheatmap(rld_cor_TC4, annotation = cluster_metadata_TC4[ c("celltype", "Treatment"), drop=F])

rld_mat_TC5 <- assay(rld_TC5)
rld_cor_TC5 <- cor(rld_mat_TC5)
pheatmap(rld_cor_TC5, annotation = cluster_metadata_TC5[ c("celltype", "Treatment"), drop=F])

#for DE analysis 
dds_TC3 <- DESeq(dds_TC3)
dds_TC4 <- DESeq(dds_TC4)
dds_TC5 <- DESeq(dds_TC5)

plotDispEsts(dds_TC3)
plotDispEsts(dds_TC4)
plotDispEsts(dds_TC5)
levels(factor(cluster_metadata_TC3$Treatment))[2]
levels(factor(cluster_metadata_TC3$Treatment))[1]
contrast_TC3 <- c("Treatment", levels(factor(cluster_metadata_TC3$Treatment))[2], levels(factor(cluster_metadata_TC3$Treatment))[1])

levels(factor(cluster_metadata_TC4$Treatment))[2]
levels(factor(cluster_metadata_TC4$Treatment))[1]
contrast_TC4 <- c("Treatment", levels(factor(cluster_metadata_TC4$Treatment))[2], levels(factor(cluster_metadata_TC4$Treatment))[1])

levels(factor(cluster_metadata_TC5$Treatment))[2]
levels(factor(cluster_metadata_TC5$Treatment))[1]
contrast_TC5 <- c("Treatment", levels(factor(cluster_metadata_TC5$Treatment))[2], levels(factor(cluster_metadata_TC5$Treatment))[1])


res_TC3 <- results(dds_TC3, 
                   contrast = contrast_TC3,
                   alpha = 0.05)

# res_TC3 <- lfcShrink(dds_TC3,
#                  contrast =  contrast_TC3,
#                  res=res_TC3)
res_TC4 <- results(dds_TC4, 
                   contrast = contrast_TC4,
                   alpha = 0.05)

res_TC5 <- results(dds_TC5, 
                   contrast = contrast_TC5,
                   alpha = 0.05)

res_tbl_TC3 <- res_TC3 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

res_tbl_TC4 <- res_TC4 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

res_tbl_TC5 <- res_TC5 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

#save all the genes
write.csv(res_tbl_TC3,
          paste0("./BI_out/", "TC3_DEseq2_results_", clusters[1], "_", levels(cluster_metadata_TC3$Treatment)[2], "_vs_", levels(cluster_metadata_TC3$Treatment)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)


write.csv(res_tbl_TC4,
          paste0("./BI_out/", "TC4_DEseq2_results_", clusters[1], "_", levels(cluster_metadata_TC4$Treatment)[2], "_vs_", levels(cluster_metadata_TC4$Treatment)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)

write.csv(res_tbl_TC5,
          paste0("./BI_out/", "TC5_DEseq2_results_", clusters[1], "_", levels(cluster_metadata_TC5$Treatment)[2], "_vs_", levels(cluster_metadata_TC5$Treatment)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)

#filter and save
padj_cutoff <- 0.05
sig_res_TC3 <- dplyr::filter(res_tbl_TC3, padj < padj_cutoff) %>%
  dplyr::arrange(padj)
sig_res_TC3
write.csv(sig_res_TC3,
          paste0("./BI_out/", "TC3_DEseq2_results_", clusters[1], "_" , levels(cluster_metadata_TC3$Treatment)[2], "_vs_", levels(cluster_metadata_TC3$Treatment)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)


padj_cutoff <- 0.05
sig_res_TC4 <- dplyr::filter(res_tbl_TC4, padj < padj_cutoff) %>%
  dplyr::arrange(padj)
sig_res_TC4
write.csv(sig_res_TC4,
          paste0("./BI_out/", "TC4_DEseq2_results_", clusters[1], "_" , levels(cluster_metadata_TC4$Treatment)[2], "_vs_", levels(cluster_metadata_TC4$Treatment)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)


padj_cutoff <- 0.05
sig_res_TC5 <- dplyr::filter(res_tbl_TC5, padj < padj_cutoff) %>%
  dplyr::arrange(padj)
sig_res_TC5
write.csv(sig_res_TC5,
          paste0("./BI_out/", "TC5_DEseq2_results_", clusters[1], "_" , levels(cluster_metadata_TC5$Treatment)[2], "_vs_", levels(cluster_metadata_TC5$Treatment)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)

#get the normalized counts
# ggplot of top genes
normalized_counts_TC3 <- counts(dds_TC3, 
                                normalized = TRUE)

# Order results by padj values
top20_sig_genes_TC3 <- sig_res_TC3 %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm_TC3 <- data.frame(normalized_counts_TC3) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes_TC3)

gathered_top20_sig_TC3 <- top20_sig_norm_TC3 %>%
  gather(colnames(top20_sig_norm_TC3)[2:length(colnames(top20_sig_norm_TC3))], key = "samplename", value = "normalized_counts")

gathered_top20_sig_TC3 <- inner_join(ei[, c("Sample", "Treatment" )], gathered_top20_sig_TC3, by = c("Sample" = "samplename"))

# plot using ggplot2
ggplot(gathered_top20_sig_TC3) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = Treatment), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes in TC3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))

# ggplot of top genes
normalized_counts_TC4 <- counts(dds_TC4, 
                                normalized = TRUE)

# Order results by padj values
top20_sig_genes_TC4 <- sig_res_TC4 %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm_TC4 <- data.frame(normalized_counts_TC4) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes_TC4)

gathered_top20_sig_TC4 <- top20_sig_norm_TC4 %>%
  gather(colnames(top20_sig_norm_TC4)[2:length(colnames(top20_sig_norm_TC4))], key = "samplename", value = "normalized_counts")

gathered_top20_sig_TC4 <- inner_join(ei[, c("Sample", "Treatment" )], gathered_top20_sig_TC4, by = c("Sample" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig_TC4) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = Treatment), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes in TC4") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))

## ggplot of top genes
normalized_counts_TC5 <- counts(dds_TC5, 
                                normalized = TRUE)

## Order results by padj values
top20_sig_genes_TC5 <- sig_res_TC5 %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm_TC5 <- data.frame(normalized_counts_TC5) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes_TC5)

gathered_top20_sig_TC5 <- top20_sig_norm_TC5 %>%
  gather(colnames(top20_sig_norm_TC5)[2:length(colnames(top20_sig_norm_TC5))], key = "samplename", value = "normalized_counts")

gathered_top20_sig_TC5 <- inner_join(ei[, c("Sample", "Treatment" )], gathered_top20_sig_TC5, by = c("Sample" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig_TC5) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = Treatment), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes in TC5") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))


#heat map of significant genes
pdf("./BI_out/TC3_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_sig_genes_heatmap.pdf")
# Extract normalized counts for only the significant genes
sig_norm_TC3 <- data.frame(normalized_counts_TC3) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res_TC3$gene)

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm_TC3[ , 2:length(colnames(sig_norm_TC3))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata_TC3[, c("Treatment", "celltype")],
         #make annotation file here
         #annotation = SampleTable[, "condition"],
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 
dev.off()

pdf("./BI_out/TC4_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_sig_genes_heatmap.pdf")
# Extract normalized counts for only the significant genes
sig_norm_TC4 <- data.frame(normalized_counts_TC4) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res_TC4$gene)

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm_TC4[ , 2:length(colnames(sig_norm_TC4))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata_TC4[, c("Treatment", "celltype")],
         #make annotation file here
         #annotation = SampleTable[, "condition"],
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 
dev.off()

pdf("./BI_out/TC5_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_sig_genes_heatmap.pdf")
# Extract normalized counts for only the significant genes
sig_norm_TC5 <- data.frame(normalized_counts_TC5) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res_TC5$gene)

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm_TC5[ , 2:length(colnames(sig_norm_TC5))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata_TC5[, c("Treatment", "celltype")],
         #make annotation file here
         #annotation = SampleTable[, "condition"],
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 
dev.off()

#volcano plot
res_table_thres_TC3 <- res_tbl_TC3 %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Volcano plot
ggplot(res_table_thres_TC3) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("Volcano plot of Epithelial significant DEGs from DESeq2 (vehicle vs. SOSMEK)") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

pdf(paste0("./BI_out/", "TC3_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_volcano.pdf"))
res_tbl_TC3$expression <- "not significant"
res_tbl_TC3$expression[res_tbl_TC3$pvalue < 0.05 & res_tbl_TC3$log2FoldChange > 0.5] <- "upregulated"
res_tbl_TC3$expression[res_tbl_TC3$pvalue < 0.05 &  res_tbl_TC3$log2FoldChange < -0.5] <- "downregulated"

top_genes <- res_tbl_TC3 %>% 
  filter(expression != 'not significant')
ggplot(data=res_tbl_TC3, aes(x=log2FoldChange, y=-log10(pvalue), col=expression)) + 
  geom_point(size=0.5) + 
  #theme_minimal() +
  theme_bw() +
  theme(text=element_text(size=11))+
  #geom_text_repel(data=head(markers_rcpp_TC3_subset, 50), aes(label=gene),nudge_x = .001 ,nudge_y = 0.001, size=3) +
  #geom_text_repel(data=markers_rcpp_TC3_subset, aes(label=gene),nudge_x = 0 ,nudge_y = 0, size=3) +
  #geom_text_repel(data=top_genes, aes(label=gene)) +
  geom_label_repel(data = head(top_genes, n= 20), 
                   mapping = aes(log2FoldChange, -log10(pvalue), label = gene), size = 2, max.overlaps=Inf)+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3))+
  #geom _text is about the size of dots
  #geom_text(position=position_dodge(width=0.9),hjust=1.5, size= 0.1)+
  #scale_color_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("blue", "dark grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), col="black", linetype="dashed") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf(paste0("./BI_out/", "TC4_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_volcano.pdf"))
res_tbl_TC4$expression <- "not significant"
res_tbl_TC4$expression[res_tbl_TC4$pvalue < 0.05 & res_tbl_TC4$log2FoldChange > 0.5] <- "upregulated"
res_tbl_TC4$expression[res_tbl_TC4$pvalue < 0.05 &  res_tbl_TC4$log2FoldChange < -0.5] <- "downregulated"

top_genes <- res_tbl_TC4 %>% 
  filter(expression != 'not significant')
ggplot(data=res_tbl_TC4, aes(x=log2FoldChange, y=-log10(pvalue), col=expression)) + 
  geom_point(size=0.5) + 
  #theme_minimal() +
  theme_bw() +
  theme(text=element_text(size=11))+
  #geom_text_repel(data=head(markers_rcpp_TC3_subset, 50), aes(label=gene),nudge_x = .001 ,nudge_y = 0.001, size=3) +
  #geom_text_repel(data=markers_rcpp_TC3_subset, aes(label=gene),nudge_x = 0 ,nudge_y = 0, size=3) +
  #geom_text_repel(data=top_genes, aes(label=gene)) +
  geom_label_repel(data = head(top_genes, n= 20), 
                   mapping = aes(log2FoldChange, -log10(pvalue), label = gene), size = 2, max.overlaps=Inf)+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3))+
  #geom _text is about the size of dots
  #geom_text(position=position_dodge(width=0.9),hjust=1.5, size= 0.1)+
  #scale_color_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("blue", "dark grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), col="black", linetype="dashed") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

pdf(paste0("./BI_out/", "TC5_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_volcano.pdf"))
res_tbl_TC5$expression <- "not significant"
res_tbl_TC5$expression[res_tbl_TC5$pvalue < 0.05 & res_tbl_TC5$log2FoldChange > 0.5] <- "upregulated"
res_tbl_TC5$expression[res_tbl_TC5$pvalue < 0.05 &  res_tbl_TC5$log2FoldChange < -0.5] <- "downregulated"

top_genes <- res_tbl_TC5 %>% 
  filter(expression != 'not significant')
ggplot(data=res_tbl_TC5, aes(x=log2FoldChange, y=-log10(pvalue), col=expression)) + 
  geom_point(size=0.5) + 
  #theme_minimal() +
  theme_bw() +
  theme(text=element_text(size=11))+
  #geom_text_repel(data=head(markers_rcpp_TC3_subset, 50), aes(label=gene),nudge_x = .001 ,nudge_y = 0.001, size=3) +
  #geom_text_repel(data=markers_rcpp_TC3_subset, aes(label=gene),nudge_x = 0 ,nudge_y = 0, size=3) +
  #geom_text_repel(data=top_genes, aes(label=gene)) +
  geom_label_repel(data = head(top_genes, n= 20), 
                   mapping = aes(log2FoldChange, -log10(pvalue), label = gene), size = 2, max.overlaps=Inf)+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3))+
  #geom _text is about the size of dots
  #geom_text(position=position_dodge(width=0.9),hjust=1.5, size= 0.1)+
  #scale_color_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("blue", "dark grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), col="black", linetype="dashed") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

TC3 <- read.csv("./BI_out/TC3_DEseq2_results_Epithelial_Vehicle_vs_SOSMEK_sig_genes.csv")
TC3_up <- subset(TC3, log2FoldChange > 0)
TC3_down <- subset(TC3, log2FoldChange < 0)
TC4 <- read.csv("./BI_out/TC4_DEseq2_results_Epithelial_Vehicle_vs_SOSMEK_sig_genes.csv")
TC4_up <- subset(TC4, log2FoldChange > 0)
TC4_down <- subset(TC4, log2FoldChange < 0)
TC5 <- read.csv("./BI_out/TC5_DEseq2_results_Epithelial_Vehicle_vs_SOSMEK_sig_genes.csv")
TC5_up <- subset(TC5, log2FoldChange > 0)
TC5_down <- subset(TC5, log2FoldChange < 0)

#for venn

set_TC3_up <- c(TC3_up$gene)
set_TC4_up<- c(TC4_up$gene)
set_TC5_up <- c(TC5_up$gene)
set_TC3_down <- c(TC3_down$gene)
set_TC4_down<- c(TC4_down$gene)
set_TC5_down <- c(TC5_down$gene)

#venn diagram for up
venn.diagram(
  x = list(set_TC3_up, set_TC4_up, set_TC5_up),
  category.names = c("TC3 upregulated" , "TC4 upregulated " , "TC5 upregulated"),
  filename = './BI_out/Epithelial_DEseq2_upregulated_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#924900'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#924900',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#924900'),
  rotation = 1
)

#venn diagram for down
venn.diagram(
  x = list(set_TC3_down, set_TC4_down, set_TC5_down),
  category.names = c("TC3 downregulated" , "TC4 downregulated " , "TC5 downregulated"),
  filename = './BI_out/Epithelial_DEseq2_downregulated_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#924900'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#924900',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#924900'),
  rotation = 1
)

#save the common genes
Epithelial_Desq_common_up <- intersect(intersect(set_TC3_up,set_TC4_up),set_TC5_up)
write.csv(Epithelial_Desq_common_up, "./BI_out/Epithelial_DEseq2_common_up.csv", row.names = F)
Epithelial_Desq_common_down <- intersect(intersect(set_TC3_down,set_TC4_down),set_TC5_down)
write.csv(Epithelial_Desq_common_down, "./BI_out/Epithelial_DEseq2_common_down.csv", row.names = F)

#save other common
Epithelial_Desq_only_TC3_TC4_common_up <- setdiff(intersect(set_TC3_up, set_TC4_up), set_TC5_up)
Epithelial_Desq_only_TC4_TC5_common_up <- setdiff(intersect(set_TC4_up, set_TC5_up), set_TC3_up)
Epithelial_Desq_only_TC3_TC5_common_up <- setdiff(intersect(set_TC3_up, set_TC5_up), set_TC4_up)

write.csv(Epithelial_Desq_only_TC3_TC4_common_up, "./BI_out/Epithelial_DEseq2_only_TC3_TC4_common_up.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC4_TC5_common_up, "./BI_out/Epithelial_DEseq2_only_TC4_TC5_common_up.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC3_TC5_common_up, "./BI_out/Epithelial_DEseq2_only_TC3_TC5_common_up.csv", row.names = F)


Epithelial_Desq_only_TC3_TC4_common_down <- setdiff(intersect(set_TC3_down, set_TC4_down), set_TC5_down)
Epithelial_Desq_only_TC4_TC5_common_down <- setdiff(intersect(set_TC4_down, set_TC5_down), set_TC3_down)
Epithelial_Desq_only_TC3_TC5_common_down <- setdiff(intersect(set_TC3_down, set_TC5_down), set_TC4_down)


write.csv(Epithelial_Desq_only_TC3_TC4_common_down, "./BI_out/Epithelial_DEseq2_only_TC3_TC4_common_down.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC4_TC5_common_down, "./BI_out/Epithelial_DEseq2_only_TC4_TC5_common_down.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC3_TC5_common_down, "./BI_out/Epithelial_DEseq2_only_TC3_TC5_common_down.csv", row.names = F)


#save
Epithelial_Desq_only_TC3_up <- setdiff(setdiff(set_TC3_up, set_TC4_up), set_TC5_up)
Epithelial_Desq_only_TC4_up <- setdiff(setdiff(set_TC4_up, set_TC5_up), set_TC3_up)
Epithelial_Desq_only_TC5_up <- setdiff(setdiff(set_TC5_up, set_TC3_up), set_TC4_up)

write.csv(Epithelial_Desq_only_TC3_up, "./BI_out/Epithelial_DEseq2_only_TC3_up.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC4_up, "./BI_out/Epithelial_DEseq2_only_TC4_up.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC5_up, "./BI_out/Epithelial_DEseq2_only_TC5_up.csv", row.names = F)


Epithelial_Desq_only_TC3_down <- setdiff(setdiff(set_TC3_down, set_TC4_down), set_TC5_down)
Epithelial_Desq_only_TC4_down <- setdiff(setdiff(set_TC4_down, set_TC5_down), set_TC3_down)
Epithelial_Desq_only_TC5_down <- setdiff(setdiff(set_TC5_down, set_TC3_down), set_TC4_down)


write.csv(Epithelial_Desq_only_TC3_down, "./BI_out/Epithelial_DEseq2_only_TC3_down.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC4_down, "./BI_out/Epithelial_DEseq2_only_TC4_down.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC5_down, "./BI_out/Epithelial_DEseq2_only_TC5_down.csv", row.names = F)
