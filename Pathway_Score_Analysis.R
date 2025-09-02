install.packages("qusage")
library("qusage")
gsets <- qusage::read.gmt("mh.all.v0.3.symbols.gmt")
cell.object <- AddModuleScore(Epithelial, features = gsets)
saveRDS(cell.object, "./BI_out/Epithelial_AddModuleScore.rds")

#Epithelials
colnames(Epithelial_AddModuleScore@meta.data)[21:70] <- names(gsets)

#features <- c(colnames(Epithelial_AddModuleScore@meta.data))[21:70]
features <- names(gsets)


#box plot comaprison
pdf(paste0("./BI_out/","Pathway_score_boxplot_Treatment_Cellline.pdf"),width = 20,height = 10)

features <- names(gsets)

library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library("dplyr")
library(ggpubr)
compare_list = list(
  c('SOSMEK', 'Vehicle')
)
for (i in features){
  print(ggplot(aes ( x = Treatment, y = Epithelial_AddModuleScore@meta.data[,i]), data = Epithelial_AddModuleScore@meta.data) +
          geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = Treatment, y = Epithelial_AddModuleScore@meta.data[,i], fill= CellLine),outlier.color="black")+
          geom_point(color = 'darkblue')+
          facet_wrap( 'CellLine', scales = "fixed") + theme_hc() +
          stat_compare_means(comparisons = compare_list,method = 't.test')+
          theme(legend.position="none")+
          ylab(i)+
          theme(text = element_text(size = 16),
                axis.text.x=element_text(size=16),
                axis.text.y=element_text(size=16),
                axis.title.y = element_text(size =16),
                axis.title.x = element_text(size=16)))    
  
}

dev.off()

pdf(paste0("./BI_out/","Pathway_score_boxplot_Cellline_Treatment.pdf"),width = 20,height = 10)

features <- names(gsets)

library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library("dplyr")
library(ggpubr)
compare_list = list(
  c('TC3','TC4'), c('TC3', 'TC5'), c('TC4', 'TC5'))
for (i in features){
  print(ggplot(aes ( x = CellLine, y = Epithelial_AddModuleScore@meta.data[,i]), data = Epithelial_AddModuleScore@meta.data) +
          geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = CellLine, y = Epithelial_AddModuleScore@meta.data[,i], fill= Treatment),outlier.color="black")+
          geom_point(color = 'darkblue')+
          facet_wrap( 'Treatment', scales = "free") + theme_hc() +
          stat_compare_means(comparisons = compare_list,method = 't.test')+
          theme(legend.position="none")+
          ylab(i)+
          theme(text = element_text(size = 16),
                axis.text.x=element_text(size=16),
                axis.text.y=element_text(size=16),
                axis.title.y = element_text(size =16),
                axis.title.x = element_text(size=16)))    
  
}

dev.off()

pdf(paste0("./BI_out/","Pathway_score_boxplot_Cellline_Treatment.pdf"),width = 20,height = 10)

features <- names(gsets)

library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library("dplyr")
library(ggpubr)
compare_list = list(
  c('TC3','TC4'), c('TC3', 'TC5'), c('TC4', 'TC5'))
for (i in features){
  print(ggplot(aes ( x = CellLine, y = Epithelial_AddModuleScore@meta.data[,i]), data = Epithelial_AddModuleScore@meta.data) +
          geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = CellLine, y = Epithelial_AddModuleScore@meta.data[,i], fill= Treatment),outlier.color="black")+
          geom_point(color = 'darkblue')+
          facet_wrap( 'Treatment', scales = "free") + theme_hc() +
          stat_compare_means(comparisons = compare_list,method = 't.test')+
          theme(legend.position="none")+
          ylab(i)+
          theme(text = element_text(size = 16),
                axis.text.x=element_text(size=16),
                axis.text.y=element_text(size=16),
                axis.title.y = element_text(size =16),
                axis.title.x = element_text(size=16)))    
  
}

dev.off()



