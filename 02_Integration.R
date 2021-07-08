library(Seurat)
library(dplyr)
library(Matrix)
library(SoupX)
library(tidyverse)
library(DoubletFinder)
library(ggplot2)
library(cowplot)

options(future.globals.maxSize = 10000 * 1024^2)

##############################################################################################################1

DTP_01_QC_MT25_Singlet@meta.data$ori <- "DTP_01"
DTP_02_QC_MT25_Singlet@meta.data$ori <- "DTP_02"
DTP_04_QC_MT25_Singlet@meta.data$ori <- "DTP_04"
DTP_05_QC_MT25_Singlet@meta.data$ori <- "DTP_05"

DTP_01_QC_MT25_Singlet@meta.data$condition <- "sound"
DTP_02_QC_MT25_Singlet@meta.data$condition <- "deep dental caries"
DTP_04_QC_MT25_Singlet@meta.data$condition <- "enamel caries"
DTP_05_QC_MT25_Singlet@meta.data$condition <- "deep dental caries"

DTP_01_QC_MT25_Singlet@meta.data$conditionforDE <- "sound and enamel caries"
DTP_02_QC_MT25_Singlet@meta.data$conditionforDE <- "deep dental caries"
DTP_04_QC_MT25_Singlet@meta.data$conditionforDE <- "sound and enamel caries"
DTP_05_QC_MT25_Singlet@meta.data$conditionforDE <- "deep dental caries"

DTP_list <- list(DTP_01_QC_MT25_Singlet, DTP_02_QC_MT25_Singlet, DTP_04_QC_MT25_Singlet, DTP_05_QC_MT25_Singlet)

for (i in 1:length(DTP_list)) {
  DTP_list[[i]] <- NormalizeData(DTP_list[[i]], verbose = FALSE)
  DTP_list[[i]] <- FindVariableFeatures(DTP_list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}

DTP.anchors <- FindIntegrationAnchors(object.list = DTP_list, dims = 1:30)

DTP.integrated <- IntegrateData(anchorset = DTP.anchors, dims = 1:30)

DTP.integrated <- ScaleData(DTP.integrated, verbose = FALSE)
DTP.integrated <- RunPCA(DTP.integrated, npcs = 30, verbose = FALSE)
DTP.integrated <- RunUMAP(DTP.integrated, reduction = "pca", dims = 1:30)
DTP.integrated <- FindNeighbors(DTP.integrated, reduction = "pca", dims = 1:20)
DTP.integrated <- FindClusters(DTP.integrated, resolution = 0.8)
p1 <- DimPlot(DTP.integrated, reduction = "umap", label = TRUE) + labs(title = "DTP.integrated")
p2 <- DimPlot(DTP.integrated, reduction = "umap", group.by = "condition") + labs(title = "DTP.integrated group by condition")
p3 <- DimPlot(DTP.integrated, reduction = "umap", group.by = "ori") + labs(title = "DTP.integrated group by ori")
plot_grid(p1, p2)
plot_grid(p1, p3)

DTP.integrated$condition <- factor(DTP.integrated$condition, levels = c("sound","enamel caries", "deep dental caries"))
DTP.integrated$conditionforDE <- factor(DTP.integrated$conditionforDE, levels = c("sound and enamel caries", "deep dental caries"))
DTP.integrated$ori <- factor(DTP.integrated$ori, levels = c("DTP01","DTP02","DTP04","DTP05"))

DTP.integrated.markers <- FindAllMarkers(DTP.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- DTP.integrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

top10 <- DTP.integrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(DTP.integrated, features = top10$gene) + NoLegend()

DTP.integrated_clustered <- DTP.integrated

DTP.integrated_clustered <- RenameIdents(DTP.integrated_clustered, '0' = "OMD+ odontoblasts", '1' = "Early odontoblasts", '2' = "OMD+ odontoblasts", 
                                         '3' = "Early odontoblasts", '4' = "OMD+ odontoblasts", '5' = "ALPL+ odontoblasts", '6' = "HSCs", 
                                         '7' = "OMD+ odontoblasts", '8' = "Fibroblasts", '9' = "Macrophages", '10' = "Endothelial cells", 
                                         '11' = "NK cells", '12' = "T cells", '13' = "B cells", '14' = "DMP1+ odontoblasts", '15' = "Glial cells", 
                                         '16' = "Plasma cells", '17' = "CD16+ monocytes", '18' = "Schwann cells",'19' = "CD103+ DCs")
DimPlot(DTP.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(DTP.integrated_clustered, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) + NoLegend()
DimPlot(DTP.integrated_clustered, reduction = "umap", label = FALSE, pt.size = 0.5)

DTP.integrated_clustered$celltype <- Idents(DTP.integrated_clustered)

DTP.integrated_clustered$celltype <- factor(DTP.integrated_clustered$celltype, levels = c("Early odontoblasts", "ALPL+ odontoblasts", 
                                                                                          "DMP1+ odontoblasts", "OMD+ odontoblasts", 
                                                                                          "Endothelial cells" ,"HSCs", "T cells", 
                                                                                          "NK cells", "B cells", "CD103+ DCs", "Plasma cells", 
                                                                                          "CD16+ monocytes", "Macrophages", "Fibroblasts", 
                                                                                          "Glial cells", "Schwann cells"))

saveRDS(DTP.integrated_clustered, file = "/path/DTP.integrated_clustered_with_soupx_reordercelltype_changeclustername.RDS")





