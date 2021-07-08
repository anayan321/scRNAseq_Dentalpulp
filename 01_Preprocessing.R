library(Seurat)
library(dplyr)
library(Matrix)
library(SoupX)
library(tidyverse)
library(DoubletFinder)

options(future.globals.maxSize = 10000 * 1024^2)

##############################################################################################################1

DTP_01 <- Read10X("/path/filtered_feature_bc_matrix")

DTP_01 <- CreateSeuratObject(counts = DTP_01, project = "DTP_01")

DTP_01[["percent.mt"]] <- PercentageFeatureSet(DTP_01, pattern = "^MT-")

DTP_01_QC <- NormalizeData(DTP_01, normalization.method = "LogNormalize", scale.factor = 10000)
DTP_01_QC <- FindVariableFeatures(DTP_01_QC, selection.method = "vst", nfeatures = 2000)
DTP_01_QC <- ScaleData(DTP_01_QC)
DTP_01_QC <- RunPCA(DTP_01_QC, features = VariableFeatures(object = DTP_01_QC))
DTP_01_QC <- FindNeighbors(DTP_01_QC, dims = 1:30)
DTP_01_QC <- FindClusters(DTP_01_QC, resolution = 0.5)
DTP_01_QC <- RunUMAP(DTP_01_QC, dims = 1:30)

saveRDS(DTP_01_QC, file="/path/DTP_01_QC_Genename.RDS")

########################################################################################################SoupX

DTP_01_QC.markers <- FindAllMarkers(DTP_01_QC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DTP_01_QC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

top10_DTP_01 <- DTP_01_QC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(DTP_01_QC, features = top10_DTP_01$gene) + NoLegend()

new.cluster.ids <- c("Osteo/Odontoblasts", "HSC", "MSC", "ALPL+ MSC", "Monocytes", "T/NK cells", "Fibroblasts", "Glial cells", "RBC", "DMP1+ osteoblasts")
DTP_01_QC_clustered <- DTP_01_QC
names(new.cluster.ids) <- levels(DTP_01_QC_clustered)
DTP_01_QC_clustered <- RenameIdents(DTP_01_QC_clustered, new.cluster.ids)
DimPlot(DTP_01_QC_clustered, reduction = "umap", label = TRUE, pt.size = 0.5) + labs(title = "DTP_01_QC_ident") + NoLegend()

DTP_01_QC_clustered$celltype <- Idents(DTP_01_QC_clustered)

DTP_01_QC_clustered_metadata <- data.frame(DTP_01_QC_clustered[["celltype"]])
DTP_01_QC_clustered_metadata[["cluster"]] <- DTP_01_QC_clustered@meta.data$RNA_snn_res.0.5

DTP_01_QC_clustered_umap <- data.frame(DTP_01_QC_clustered[["umap"]]@cell.embeddings)
DTP_01_QC_clustered_metadata[["umap1"]] <- DTP_01_QC_clustered_umap$UMAP_1
DTP_01_QC_clustered_metadata[["umap2"]] <- DTP_01_QC_clustered_umap$UMAP_2

toc1 = Seurat::Read10X("/path/filtered_feature_bc_matrix")
tod1 = Seurat::Read10X("/path/raw_feature_bc_matrix")
sc1 = SoupChannel(tod1, toc1)
sc1

sc1 = setClusters(sc1, setNames(DTP_01_QC_clustered_metadata$cluster, rownames(DTP_01_QC_clustered_metadata)))

sc1 = setDR(sc1, DTP_01_QC_clustered_metadata[colnames(sc1$toc), c("umap1", "umap2")])

library(ggplot2)
dd = DTP_01_QC_clustered_metadata[colnames(sc1$toc),]
mids = aggregate(cbind(umap1,umap2) ~ celltype,data=dd,FUN=mean)
gg = ggplot(dd,aes(umap1,umap2)) + 
  geom_point(aes(colour=celltype),size=0.2) +
  geom_label(data=mids,aes(label=celltype)) +
  ggtitle('DTP_01') +
  guides(colour = guide_legend(override.aes = list(size=1)))
plot(gg)

dd$IGKC = sc1$toc['IGKC',]
gg = ggplot(dd,aes(umap1,umap2)) +
  geom_point(aes(colour="IGKC">0))
plot(gg)

gg = plotMarkerMap(sc1, "IGKC")
plot(gg)

head(sc1$soupProfile[order(sc1$soupProfile$est, decreasing = TRUE), ], n = 20)

plotMarkerDistribution(sc1)

igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")

HBGenes = c("HBB", "HBA2", "HBA1")

useToEst = estimateNonExpressingCells(sc1, nonExpressedGeneList = list(c(IG = igGenes, HB =HBGenes)), 
                                      clusters = FALSE)

plotMarkerMap(sc1, geneSet = c(igGenes, HBGenes), useToEst = useToEst)

useToEst = estimateNonExpressingCells(sc1, nonExpressedGeneList = list(c(IG = igGenes, HB =HBGenes)))
plotMarkerMap(sc1, geneSet = c(igGenes, HBGenes), useToEst = useToEst)

sc1 = calculateContaminationFraction(sc1, list(c(IG = igGenes, HB =HBGenes)), useToEst = useToEst)

out1 = adjustCounts(sc1)

###################################################################################################2

DTP_01 <- CreateSeuratObject(counts = out1, project = "DTP_01", min.cells = 3, min.features = 200)

DTP_01[["percent.mt"]] <- PercentageFeatureSet(DTP_01, pattern = "^MT-")

DTP_01_QC_MT25 <- subset(DTP_01, subset = percent.mt < 25)

DTP_01_QC_MT25 <- NormalizeData(DTP_01_QC_MT25, normalization.method = "LogNormalize", scale.factor = 10000)
DTP_01_QC_MT25 <- FindVariableFeatures(DTP_01_QC_MT25, selection.method = "vst", nfeatures = 2000)
DTP_01_QC_MT25 <- ScaleData(DTP_01_QC_MT25)
DTP_01_QC_MT25 <- RunPCA(DTP_01_QC_MT25, features = VariableFeatures(object = DTP_01_QC_MT25))
DTP_01_QC_MT25 <- FindNeighbors(DTP_01_QC_MT25, dims = 1:30)
DTP_01_QC_MT25 <- FindClusters(DTP_01_QC_MT25, resolution = 0.5)
DTP_01_QC_MT25 <- RunUMAP(DTP_01_QC_MT25, dims = 1:30)

#Doublet removal
sweep.res.list <- paramSweep_v3(DTP_01_QC_MT25, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
pK_value <- find.pK(sweep.stats)
pK_value<-pK_value%>%filter(MeanBC==max(MeanBC))%>%dplyr::select(pK)%>%unlist%>%as.character()%>%as.numeric()
annotations <- DTP_01_QC_MT25@meta.data$seurat_clusters
annotations <- as.character(annotations)
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*length(DTP_01_QC_MT25%>%colnames))
DTP_01_QC_MT25 <- doubletFinder_v3(DTP_01_QC_MT25 , PCs = 1:30, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(DTP_01_QC_MT25@meta.data)[length(DTP_01_QC_MT25@meta.data)] <- "Doublet"

DTP_01_QC_MT25_Singlet <- subset(DTP_01 , cells =  rownames(DTP_01_QC_MT25@meta.data[DTP_01_QC_MT25@meta.data$Doublet == "Singlet",]))

saveRDS(DTP_01_QC_MT25_Singlet, file="/path/DTP_01_QC_MT25_Singlet_norm_with_soupx.RDS")





