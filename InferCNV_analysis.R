#this part runs InfeCNV analysis to identify the epithelial tumor cells, setting B & T cells as the reference
#save the matrix
counts.matrix=BI_filter_YFP@assays$RNA@data; dim(counts.matrix)

cellcounts <- GetAssayData(BI.filter_YFP[['RNA']],slot='counts')
saveRDS(cellcounts,'./BI_out/cellcounts.rds')

#downsampling, split cell indices by identity
idx <- split(meta, f=meta$MajorCellType)
#idx1 <- unsplit(idx, meta$MajorCellType)
cs_keep <- lapply(idx, function(i) {
  print(class(i))
  i=i[1:100,]
})  

meta_downsample=data.frame() 
for (i in cs_keep){
  meta_downsample=rbind(meta_downsample,i)
}

#downsampling, another way
meta_downsample1 <- data.frame()
for(c in unique(meta$MajorCellType)){
  print(meta[meta$MajorCellType==c,][1:100,])
  meta_downsample1 <- rbind(meta_downsample1,meta[meta$MajorCellType==c,][1:100,])
}  
write.csv(meta_downsample,file=paste0("./BI_out/meta_downsample.csv"))
write.csv(meta_downsample1,file=paste0("./BI_out/meta_downsample1.csv"))
#meta_downsample <- read.csv(file = "./BI_out/meta_downsample.csv", row.names = 1)


#make seurat object for smaller meta
BI_filter_YFP_downsample <- CreateSeuratObject(counts = counts.matrix[, rownames(meta_downsample)], 
                                               project = "BI_downsample", 
                                               assay = "RNA", 
                                               meta.data = meta_downsample)

#get the smaller matrix
counts.matrix.small=BI_filter_YFP_downsample@assays$RNA@data; dim(counts.matrix.small)

#save
saveRDS(counts.matrix.small, "./BI_out/counts.matrix.small.rds")

#create infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/home/nasim/BI/BI_out/counts.matrix.small.rds",
                                    annotations_file="/home/nasim/BI/BI_out/annotation_down.txt",
                                    delim="\t",
                                    gene_order_file="/home/nasim/BI/BI_out/gene_pos.txt",
                                    ref_group_names=c("T", "B"))

#run infercnv
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
)