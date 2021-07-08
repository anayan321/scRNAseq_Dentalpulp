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

setwd("/home/shared_data_vclab_2/SCRNAseq_test/SCRNAseq_test_sarintip/scDP_2020/NATMI")

DTP.integrated_clustered <- readRDS(file = "/path/DTP.integrated_clustered_with_soupx_reordercelltype_changeclustername.RDS")

DTP.integrated_clustered@meta.data$conditionori <- DTP.integrated_clustered@meta.data$condition
DTP.integrated_clustered@meta.data$condition <- DTP.integrated_clustered@meta.data$conditionforDE


#Deep dental caries
DTP.integrated_clustered_subset_deepcaries <- subset(DTP.integrated_clustered , cells =  rownames(DTP.integrated_clustered@meta.data[DTP.integrated_clustered@meta.data$condition == "deep dental caries",]))

DefaultAssay(DTP.integrated_clustered_subset_deepcaries) <- "RNA"
#If you use SCT, you have to normalize data before as.matrix
DTP.integrated_clustered_subset_deepcaries <- NormalizeData(DTP.integrated_clustered_subset_deepcaries, normalization.method = "LogNormalize", scale.factor = 10000)

DTP.integrated_clustered_subset_deepcaries_counts <- as.matrix(GetAssayData(DTP.integrated_clustered_subset_deepcaries, slot = "data")[, WhichCells(DTP.integrated_clustered_subset_deepcaries)])
DTP.integrated_clustered_subset_deepcaries_counts <- as.data.frame(DTP.integrated_clustered_subset_deepcaries_counts)

DTP.integrated_clustered_subset_deepcaries_meta <- data.frame(DTP.integrated_clustered_subset_deepcaries[["celltype"]])
DTP.integrated_clustered_subset_deepcaries_meta <- rownames_to_column(DTP.integrated_clustered_subset_deepcaries_meta, var = "barcode")
colnames(DTP.integrated_clustered_subset_deepcaries_meta)[colnames(DTP.integrated_clustered_subset_deepcaries_meta) == 'celltype'] <- 'annotation'

write.table(DTP.integrated_clustered_subset_deepcaries_counts, "/path/NATMI-master/all.deep.sc.em.txt", row.names = TRUE, quote = F, sep = "\t")
write.table(DTP.integrated_clustered_subset_deepcaries_meta, "/path/NATMI-master/all.deep.sc.ann.txt", row.names = FALSE, quote = F, sep = "\t")


#sound and enamel caries
DTP.integrated_clustered_subset_sound <- subset(DTP.integrated_clustered , cells =  rownames(DTP.integrated_clustered@meta.data[DTP.integrated_clustered@meta.data$condition == "sound and enamel caries",]))

DefaultAssay(DTP.integrated_clustered_subset_sound) <- "RNA"
#If you use SCT, you have to normalize data before as.matrix
DTP.integrated_clustered_subset_sound <- NormalizeData(DTP.integrated_clustered_subset_sound, normalization.method = "LogNormalize", scale.factor = 10000)

DTP.integrated_clustered_subset_sound_counts <- as.matrix(GetAssayData(DTP.integrated_clustered_subset_sound, slot = "data")[, WhichCells(DTP.integrated_clustered_subset_sound)])
DTP.integrated_clustered_subset_sound_counts <- as.data.frame(DTP.integrated_clustered_subset_sound_counts)

DTP.integrated_clustered_subset_sound_meta <- data.frame(DTP.integrated_clustered_subset_sound[["celltype"]])
DTP.integrated_clustered_subset_sound_meta <- rownames_to_column(DTP.integrated_clustered_subset_sound_meta, var = "barcode")
colnames(DTP.integrated_clustered_subset_sound_meta)[colnames(DTP.integrated_clustered_subset_sound_meta) == 'celltype'] <- 'annotation'

write.table(DTP.integrated_clustered_subset_sound_counts, "/path/NATMI-master/all.sound.sc.em.txt", row.names = TRUE, quote = F, sep = "\t")
write.table(DTP.integrated_clustered_subset_sound_meta, "/path/NATMI-master/all.sound.sc.ann.txt", row.names = FALSE, quote = F, sep = "\t")

