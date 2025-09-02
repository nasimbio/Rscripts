# ------------------------------------------------------------------------------
#title: "BCR Clonotype Analysis from Single-Cell RNA-seq"
#Author: Nasim Rahmatpour
#Date: "8/12/2022"
# ------------------------------------------------------------------------------

# Project Overview

#This analysis extracts and summarizes BCR clonotype information from 10x Genomics single-cell V(D)J sequencing data. The input includes:
# - `all_contig_annotations.csv` to identify and separate immunoglobulin heavy (IGH), kappa (IGK), and lambda (IGL) chains
# - `clonotypes.csv` to quantify and visualize clonotype frequency and chain usage

#The results help characterize the diversity and dominance of B cell clones and light chain pairing within the sample.

---
  
###Task 1: Parsing and Saving Immunoglobulin Chain Data 
  auto_basic <- function(BCR_List = NA,
                         Sample_List = NA,
                         high_confidence="true",
                         productive="true",
                         output="/home/nasim/altarion/ABIO-NGS085_BCR/analysis/"){
    
    library(dplyr,verbose = F)
    BCR_PATH = BCR_List
    Sample_List = Sample_List
    high_confidence_p = high_confidence
    productive_p = productive
    output_parameter = output
    
    
    for (i in 1:length(BCR_PATH)) {
      bcr_tmp <- read.csv(paste0(BCR_PATH[i]))
      bcr_tmp$barcode <- paste0("Data_",Sample_List[i],"_",bcr_tmp$barcode)
      assign(x = paste0("bcr_data_",Sample_List[i]),value = bcr_tmp)
    }
    
    if (length(Sample_List) == 1) {
      m <- function(x){eval(as.name(x))}
      bcr_data <- m(ls(pattern="bcr_data_"))
    }else {
      m <- function(x){eval(as.name(x))}
      the_bcr <- ls(pattern="bcr_data_")[1]
      other_bcr <- ls(pattern="bcr_data_")[-1]
      bcr_data <- m(the_bcr)
      for (l in other_bcr) {
        bcr_data <- rbind(bcr_data,m(l))
      }
      
    }
    
    
    bcr_data <- subset(bcr_data,high_confidence == "True"|high_confidence == "true")
    bcr_data <- subset(bcr_data,productive == "True"|productive == "true")
    
    bcr_data$Clonotype <- bcr_data$cdr3_nt
    
    #Subset_IGH, IGK, IGL
    for(x in unique(bcr_data$chain)){
      assign(x, subset(bcr_data,chain == x) %>% group_by(barcode) %>% arrange(desc(umis)) %>% dplyr::slice(1:1))
    }
    if(exists("IGH")==TRUE){
      write.csv(x = IGH,file = paste0(output_parameter,"BCR_RawData_IGH.csv"))
    }else{
      print("No IGH")
    }
    
    if(exists("IGK")==TRUE){
      write.csv(x = IGK,file = paste0(output_parameter,"BCR_RawData_IGK.csv"))
    }else{
      print("No IGK")
    }
    
    if(exists("IGL")==TRUE){
      write.csv(x = IGL,file = paste0(output_parameter,"BCR_RawData_IGL.csv"))
    }else{
      print("No IGL")
    }
    
  }  
  
  
  BCRList <- c("NGS008_BCR_7A","NGS008_BCR_7B","NGS008_BCR_7C","NGS008_BCR_7D","NGS008_BCR_7E","NGS008_BCR_7F")
  
  InputDirList <- c("./NGS008_BCR_7A/all_contig_annotations.csv",
                    "./NGS008_BCR_7B/all_contig_annotations.csv",
                    "./NGS008_BCR_7C/all_contig_annotations.csv",
                    "./NGS008_BCR_7D/all_contig_annotations.csv",
                    "./NGS008_BCR_7E/all_contig_annotations.csv",
                    "./NGS008_BCR_7F/all_contig_annotations.csv")
  
  #calling the fuction, makes 3 csv file
  BCR.bird_auto_basic(BCR_List = InputDirList,
                      Sample_List = BCRList,
                      output = "./batch7_outputs_BCR/")
  
  #load the seurat object
  data = immune.combined
  
  library(dplyr)
  VDJChainList <- list.files("./batch7_outputs_BCR/",pattern = "BCR_RawData_")
  nColMetadata <- ncol(data@meta.data) + 1
  for (i in VDJChainList) {
    DF_Barcode <- as.data.frame(row.names(data@meta.data))
    colnames(DF_Barcode) <- "barcode"
    
    VDJ_Tmp <- read.csv(paste0("./batch7_outputs_BCR/",i))
    VDJ_Tmp$barcode <- gsub("BCR_", "", gsub("Data_NGS008_", "PBMC", VDJ_Tmp$barcode))
    
    ChainType <- as.character(unique(VDJ_Tmp$chain))
    DF_Barcode <- left_join(DF_Barcode, VDJ_Tmp, by = "barcode")
    DF_Barcode[is.na(DF_Barcode)] <- "0"
    data$tmp <- DF_Barcode$Clonotype
    colnames(data@meta.data)[nColMetadata] <- paste0("ClonoType",ChainType)
    nColMetadata <- nColMetadata+1
  }
  
  #get Meta
  MetaData <- data@meta.data
  MetaData[is.na(MetaData)] <- "0"
  
  #add column for IGH
  MetaData$IGH_positive <- ""
  MetaData$IGH_positive[MetaData[, 'ClonoTypeIGH']==0] <- 0
  MetaData$IGH_positive[MetaData[, 'ClonoTypeIGH']!=0] <- 1
  IGH <- table(MetaData$Annotation, MetaData$IGH_positive)
  openxlsx::write.xlsx(IGH, "./batch7_outputs_BCR/CellCountIGH.xlsx",row.names = TRUE)
  
  #add column for IGK
  MetaData$IGK_positive <- ""
  MetaData$IGK_positive[MetaData[, 'ClonoTypeIGK']==0] <- 0
  MetaData$IGK_positive[MetaData[, 'ClonoTypeIGK']!=0] <- 1
  IGK <- table(MetaData$Annotation, MetaData$IGK_positive)
  openxlsx::write.xlsx(IGH, "./batch7_outputs_BCR/CellCountIGK.xlsx",row.names = TRUE)
  
  #add column for IGL
  
  MetaData$IGL_positive <- ""
  MetaData$IGL_positive[MetaData[, 'ClonoTypeIGL']==0] <- 0
  MetaData$IGL_positive[MetaData[, 'ClonoTypeIGL']!=0] <- 1
  IGL <- table(MetaData$Annotation, MetaData$IGL_positive)
  openxlsx::write.xlsx(IGH, "./batch7_outputs_BCR/CellCountIGL.xlsx",row.names = TRUE)
  
  data@meta.data <- MetaData
  
  pdf("./batch7_outputs_BCR/IGH_IGK_IGL.pdf")
  DimPlot(data,group.by = "IGH_positive",pt.size = 1, cols= c("grey", "red"))+coord_fixed() + NoLegend()
  DimPlot(data,group.by = "IGK_positive",pt.size = 1, cols= c("grey", "red"))+coord_fixed() + NoLegend()
  DimPlot(data,group.by = "IGL_positive",pt.size = 1, cols= c("grey", "red"))+coord_fixed() + NoLegend()
  dev.off()
  
###Task 2: Clonotype Frequency and Chain Usage Visualization 
  
  clonotype <- read.csv("/home/nasim/altarion/ABIO-NGS085_BCR/ABIO-NGS085_BCR/clonotypes.csv")
  clonotype$chain <- apply(clonotype,1,function(x) str_split(x[4], "\\:"))
  clonotype$chain <- apply(clonotype,1,function(x) x[[6]][[1]][1])
  
  
  #count
  clonotype$NumberIGH <- ""
  clonotype$NumberIGH <- str_count(clonotype$cdr3s_aa, "IGH")
  clonotype$NumberIGK <- ""
  clonotype$NumberIGK <- str_count(clonotype$cdr3s_aa, "IGK")
  clonotype$NumberIGL <- ""
  clonotype$NumberIGL <- str_count(clonotype$cdr3s_aa, "IGL")
  
  
  #bar plot
  clonotype_top <- clonotype[1:100, ]
  pdf(paste0("/home/nasim/altarion/ABIO-NGS085_BCR/analysis/Barplot_BCR.pdf"), width=20, height=10)
  ggplot(aes(x = reorder(clonotype_id, -frequency), y = frequency), data = clonotype_top) +
    geom_bar(stat="identity", color="blue", fill="blue") +
    theme(legend.position="none")+
    xlab("Clonotype_id")+
    theme(text = element_text(size = 8),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y = element_text(size =8),
          axis.title.x = element_text(size=8))+
    theme(axis.text.x=element_text(angle = 90))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    coord_cartesian(expand = FALSE)
  #theme_bw()
  
  dev.off()
  
  #pie chart
  pdf("/home/nasim/altarion/ABIO-NGS085_BCR/analysis/Pie_BCR.pdf")
  clonotype_chain <- as.data.frame(clonotype[,7:9])
  clonotype_chain <- apply(clonotype_chain,2,sum)
  clonotype_chain <- as.data.frame(clonotype_chain)
  colnames(clonotype_chain) <- "value"
  clonotype_chain$BCR_chain <- ""
  clonotype_chain$BCR_chain <- c("IGH", "IGK", "IGL")
  ggplot(clonotype_chain, aes(x="", y=value, fill=BCR_chain)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    #theme_hc()+
    labs("BCR_chain")+
    
    
    
    #coord_polar(theta = "y", direction = -1)+
    geom_text(aes(label = value),
              position = position_stack(vjust = 0.6), color="red") +
    scale_fill_manual(values=c("#000000","#6db6ff","#ffff6d"))+ theme_void()
  dev.off()