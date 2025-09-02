library(Seurat)
library(CellChat)
library(circlize)
library(ComplexHeatmap)
#function to make cell chat object
process_seurat_to_cellchat <- function(seurat_object, output_path, signal, lps, timepoint) {
  # Subset the Seurat object based on dynamic arguments
  subset_obj <- subset(seurat_object, subset = Signal == signal & lps == LPS & timepoint==Timepoint)
  print(dim(subset_obj))
  
  
  # Set the identity class of the subsetted Seurat object
  Idents(subset_obj) <- "Annotation"
  
  # Prepare Seurat object for CellChat
  data.input <- subset_obj[["RNA"]]$data
  labels <- Idents(subset_obj)
  meta <- data.frame(labels = labels, row.names = names(labels))
  
  
  # Create the CellChat object from the Seurat object
  cellChat <- createCellChat(object = subset_obj, group.by = "ident", assay = "RNA")
  
  # Set the ligand-receptor interaction database (Human in this case)
  CellChatDB <- CellChatDB.human  # Use CellChatDB.mouse for mouse data
  showDatabaseCategory(CellChatDB)
  cellChat@DB <- CellChatDB
  
  # Preprocessing for expression data
  cellChat <- subsetData(cellChat)
  
  # Setup parallel processing
  options(future.globals.maxSize = 2048 * 1024^2)
  future::plan("multisession", workers = 4)
  
  # Identify over-expressed genes and interactions
  cellChat <- identifyOverExpressedGenes(cellChat)
  cellChat <- identifyOverExpressedInteractions(cellChat)
  
  # Start the computation and measure execution time
  ptm <- Sys.time()
  cellChat <- computeCommunProb(cellChat, type = "triMean")
  cellChat <- filterCommunication(cellChat, min.cells = 10)
  
  # Infer cell-cell communication at signaling pathway level
  cellChat <- computeCommunProbPathway(cellChat)
  
  # Calculate the aggregated cell-cell communication network
  cellChat <- aggregateNet(cellChat)
  
  # Measure execution time
  execution_time <- Sys.time() - ptm
  print(as.numeric(execution_time, units = "secs"))
  cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")
  
  # Save the CellChat object
  saveRDS(cellChat, file = output_path)
  
  message("CellChat object has been saved to ", output_path)
}


# call function
process_seurat_to_cellchat(
  seurat_object = `3-Data_Annotattion`,
  output_path = "/home/nasim/feinstein1/CellChat/CellChat_Dex_1_Post.rds",
  signal = "Dex" ,
  timepoint = "Post",      
  lps = 1                
  
)

#function to compare to cell chat object and make the plots
library(CellChat)
library(ggplot2)
library(grid)

perform_cellchat_comparison <- function(output_dir, query_object, reference_object ) {
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE)
  
  # Merge the two CellChat objects into one (query and reference)
  object.list <- list(Pink_Post_1 = reference_object, Dex_Post_1 = query_object)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  
  # Total Interaction Comparison
  pdf(file.path(output_dir, "Total_interaction_comparison.pdf"))
  #ptm = Sys.time()
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  print(gg1 + gg2)
  dev.off()
  
  # gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
  # gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
  # total_interaction_plot <- gg1 + gg2
  # ggsave(file.path(output_dir, "Total_interaction_comparison.pdf"), plot = total_interaction_plot)
  
  # Differential Interaction Heatmap
  pdf(file.path(output_dir, "Differential_interaction_heatmap.pdf"), width = 10, height = 10)
  gg3 <- netVisual_heatmap(cellchat, font.size = 12, width = 8, height = 12, font.size.title = 14)
  gg4 <- netVisual_heatmap(cellchat, measure = "weight", font.size = 12, width = 8, height = 12, font.size.title = 14)
  print(gg3 + gg4)
  dev.off()
  
  
  # Outgoing Signals Heatmap
  pdf(file.path(output_dir, "outgoing_incoming_signals.pdf"), width = 10, height = 10)
  i <- 1
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], font.size = 6, width = 8,height = 12, font.size.title = 14)
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], font.size =6, width = 8,height = 12, font.size.title = 14)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], font.size = 6, width = 8,height = 12, font.size.title = 14, color.heatmap = "GnBu")
  ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], font.size = 6, width = 8,height = 12, font.size.title = 14, color.heatmap = "GnBu")
  draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))
  
  ht5 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], font.size = 6, width = 8,height = 12, font.size.title = 14, color.heatmap = "OrRd")
  ht6 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], font.size =6, width = 8,height = 12, font.size.title = 14, color.heatmap = "OrRd")
  draw(ht5 + ht6, ht_gap = unit(0.5, "cm"))
  
  dev.off()
  
  message("All comparison figures have been saved to: ", output_dir)
}


# call function
perform_cellchat_comparison(
  query_object = CellChat_Dex_1_Post,
  reference_object =CellChat_Pink_1_Post,
  output_dir = "/home/nasim/feinstein1/CellChat/CellChat_Inter/d/Dex_1_Post_vs_Pink_1_Post"
)
